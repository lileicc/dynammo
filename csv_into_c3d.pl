#!/usr/bin/perl -w

#
# Usage:
#
# ./csv_into_c3d.pl input.csv input.c3d output.c3d
#


if (@ARGV != 3
 ||!($ARGV[0] =~ /\.csv$/)
 ||!($ARGV[1] =~ /\.c3d$/)
 ||!($ARGV[2] =~ /\.c3d$/)
 || ($ARGV[2] eq $ARGV[1]) ) {
	print STDERR "Usage:\n$0 <input.csv> <input.c3d> <output.c3d>\n (input and output must be different.)\n";
	exit 0;
}

my $datafile = $ARGV[0];
my $infile = $ARGV[1];
my $outfile = $ARGV[2];

print STDERR "Cramming data from '$datafile' into '$infile' and saving as '$outfile'...\n";

open(DATA, "<", $datafile) || die "ERROR: Can't open $datafile.";

open(INPUT, "<", $infile) || die "ERROR: Can't open $infile.";

if ($infile =~ /'/ || $outfile =~ /'/) {
	die "Filenames $infile or $outfile contain single quotes, which I worry about.\n";
}
#TODO: better way to copy a file?
system("cp '$infile' '$outfile'");

open(OUTPUT, "+<", $outfile) || die "ERROR: Can't open $outfile for read/write.";

#---------------------------------------------------------
#Read the header:
my $header;
read(INPUT,$header,512) == 512 || die "ERROR: Can't read header.";

my (
 $param_block,
 $c3d_key,
 $num_3d_points,
 $num_analog_measurements,
 $first_frame,
 $last_frame,
 $maxiumum_interpolation_gap, #??
 $scale_factor,
 $data_block,
 $num_analog_samples,
 $frame_rate,
 #$reserved,
 $range_key,
 $range_block,
 $event_key,
 $event_count,
 $rest
) = unpack( "CCSSSSSfSSfx[270]SSSSH*", $header);
$c3d_key == 0x50 || die "ERROR: Expected c3d magic value.";
print STDERR "$num_3d_points points.\n";
print STDERR "$num_analog_measurements analog measurements.\n";
print STDERR "Frame range: $first_frame - $last_frame\n";
print STDERR "Scale factor: $scale_factor";
if ($scale_factor < 0.0) {
	print STDERR " (data stored as floats)\n";
} else {
	print STDERR " (data stored as ints)\n";
}
print STDERR "Frame rate: $frame_rate Hz.\n";
print STDERR "Analog samples/frame: $num_analog_samples.\n";
if ($range_key == 12345) {
	print STDERR "File has label and range data.\n";
} elsif ($range_key != 0) {
	print STDERR "WARNING: unexpected range key value.\n";
}
if ($event_key == 12345) {
	print STDERR "File has $event_count events.\n";
} elsif ($event_key != 0) {
	print STDERR "WARNING: unexpected range key value.\n";
}

#---------------------------------------------------------
#Read the parameter section:
print STDERR "Reading params starting at block $param_block.\n";
seek(INPUT,($param_block-1) * 512, 0) || die "ERROR: Couldn't seek.";
read(INPUT,$params,4) == 4 || die "ERROR: couldn't read params start.";
my (
$key1,
$key2,
$num_param_blocks,
$proc_type
) = unpack("CCCC", $params);
$proc_type == 83 + 1 || die "ERROR: Got processor type $proc_type, wanted 84 (intel).";

my %id_to_group;
my %groups;
my %params;

#Read group or parameter information:
while (1) {
	read(INPUT,$head,2) == 2 || die "ERROR: couldn't read param.";
	my ($len, $id) = unpack("cc", $head);
	if ($len == 0) {
		print STDERR "WARNING: Zero length group/param record; last offset should have been zero.\n";
		($id == 0) || die "ERROR: ID in zero-length group is $id. This seems like a problem!";
		last;
	}
	my $locked = 0;
	#if $len is negative, the group/parameter is 'locked':
	if ($len < 0) { $locked = 1; $len *= -1; }
	my $offset = 0; #offset to next element.
	if ($id > 0) {
		#This is a parameter.
		my $parameter;
		read(INPUT, $parameter, $len + 2 + 1 + 1) == $len + 2 + 1 + 1 || die "ERROR: couldn't read parameter defn.";
		my ($name, $data_type, $dimensions);
		($name, $offset, $data_type, $dimensions) = unpack("A[$len]scc", $parameter);
		my @dims = ();
		my $dim_prod = 1;
		for (1 .. $dimensions) {
			read(INPUT, $parameter, 1) == 1 || die "ERROR: couldn't read parameter dimensions.";
			my $d = unpack("C", $parameter);
			push @dims, $d;
			$dim_prod *= $d;
		}
		#Read the data:
		my $format = "";
		if      ($data_type == -1) {
			$format = "A";
		} elsif ($data_type ==  1) {
			$format = "c";
		} elsif ($data_type ==  2) {
			$format = "s";
		} elsif ($data_type ==  4) {
			$format = "f";
		} else {
			die "Unknown data type $data_type."
		}
		#DEBUG:print "name:$name dims:" . join(", ", @dims) . " dim_prod:$dim_prod, abs(data_type):" . abs($data_type) . "\n";
		($dim_prod >= 0) || die "Negative dim_prod ($dim_prod).";
		read(INPUT, $parameter, $dim_prod * abs($data_type)) == $dim_prod * abs($data_type) || die "ERROR: couldn't read parameter data.";
		my @data = unpack($format x $dim_prod, $parameter);
		#Read the description:
		read(INPUT, $parameter, 1) == 1 || die "ERROR: couldn't read description length.";
		my $desc_len = unpack("c", $parameter);
		my $desc;
		read(INPUT, $desc, $desc_len) == $desc_len || die "ERROR: couldn't read description.";
		#Store everything:
		$params{$id}{$name}{"type"} = $data_type;
		$params{$id}{$name}{"dimensions"} = \@dims;
		$params{$id}{$name}{"data"} = \@data;
		print STDERR "Read parameter '$name' with type '$data_type' dimension (" . join(",", @dims) . "), and data [" . join(", ", @data) . "]\n";
		#Check offset:
		($offset == 0 || $desc_len + 2 + 1 + 1 + $dimensions + $dim_prod * abs($data_type) + 1 == $offset) || die "ERROR: Unexpected offset $offset, with $dimensions dimensions, $dim_prod elements, $data_type data type, and $desc_len description characters.";
	} else {
		#This is a group.
		my $group; #group data.
		read(INPUT, $group, $len + 2 + 1) == $len + 2 + 1 || die "ERROR: couldn't read group defn.";
		my ($name, $desc_len);
		($name, $offset, $desc_len) = unpack("A[$len]sc", $group);
		my $desc;
		read(INPUT, $desc, $desc_len) == $desc_len || die "ERROR: couldn't read group description.";
		$id_to_group{-$id} = $name;
		$groups{$name}{"id"} = -$id;
		$groups{$name}{"desc"} = $desc;
		$groups{$name}{"locked"} = $locked;
		print STDERR "Read group '$name' with description '$desc'\n";
		($offset == 0 || $desc_len + 3 == $offset) || die "ERROR: Unexpected offset $offset, when description is $desc_len characters.";
	}
	#Assume everything is adjacent, so we don't have to skip characters:
	if ($offset == 0) { last; }
}

#---------------------------------------------------------
#Read the [point] data (finally!):

{ #Check the starting block versus params:
	my @d = @{ $params{$groups{"POINT"}{"id"}}{"DATA_START"}{"data"} };
	@d == 1 || die "ERROR: expecting 1 data element for POINT:DATA_START, have " . scalar(@d) . ".\n";
	my $s = shift @d;
	$s == $data_block || die "ERROR: Data block from header (== $data_block) does not match POINT:DATA_START (== $s).";
}

{ #Check the point scale:
	my @d = @{ $params{$groups{"POINT"}{"id"}}{"SCALE"}{"data"} };
	@d == 1 || die "ERROR: expecting 1 data element for POINT:SCALE, have " . scalar(@d) . ".\n";
	my $s = shift @d;
	$s == $scale_factor || die "ERROR: Scale factor from header (== $scale_factor) does not match POINT:SCALE (== $s).";
}

{ #Check the frames:
	my @d = @{ $params{$groups{"POINT"}{"id"}}{"FRAMES"}{"data"} };
	@d == 1 || die "ERROR: expecting 1 data element for POINT:FRAMES, have " . scalar(@d) . ".\n";
	my $s = shift @d;
	$s == $last_frame - $first_frame + 1 || die "ERROR: Frames from header ($first_frame to $last_frame) does not match POINT:FRAMES (== $s).";
}

{ #Check the number of points:
	my @d = @{ $params{$groups{"POINT"}{"id"}}{"USED"}{"data"} };
	@d == 1 || die "ERROR: expecting 1 data element for POINT:USED, have " . scalar(@d) . ".\n";
	my $s = shift @d;
	$s == $num_3d_points || die "ERROR: Points from header (== $num_3d_points) does not match POINT:USED (== $s).";
}


print STDERR "Dumpping data starting at block $data_block.\n";
seek(INPUT,($data_block-1) * 512, 0) || die "ERROR: Couldn't seek.";
sysseek(OUTPUT,($data_block-1) * 512, 0) || die "ERROR: Couldn't seek output.";

$scale_factor < 0 || die "ERROR: Integer point data not supported yet.";

$num_analog_measurements == 0 || die "ERROR: Analog measurements not supported yet.";

my @labels = ();
{ #construct @labels from POINT:LABELS.
	my @d = @{ $params{$groups{"POINT"}{"id"}}{"LABELS"}{"data"} };
	my @dims = @{ $params{$groups{"POINT"}{"id"}}{"LABELS"}{"dimensions"} };
	@dims == 2 || die "ERROR: Expecting 2d labels array.\n";
	for (1 .. $dims[1]) {
		my $acc = "";
		for (1 .. $dims[0]) {
			$acc .= shift(@d);
		}
		$acc =~ s/\s+$//;
		push @labels, $acc;
	}
}
@labels == $num_3d_points || die "ERROR: expecting each point to have a label; but have " . scalar(@labels) . " labels and $num_3d_points points.";

#$slots[$i] says where to put point $i. We do this so we can have a canonical marker ordering (alphabetical, in fact.)
my @slots;
{ #build slots, and re-order labels also.
	my %sources;
	for (0 .. (@labels-1)) {
		$sources{$labels[$_]} = $_;
		push @slots, -1;
	}
	my $i = 0;
	foreach (sort @labels) {
		$slots[$sources{$_}] = $i;
		++$i;
	}
	my @new_labels;
	for (0 .. @labels-1) {
		$new_labels[$slots[$_]] = $labels[$_];
	}
	@labels = @new_labels;
}
print STDERR "Reordering 0 .. " . scalar(@slots) . " to " . join(", ", @slots) . "\n";
print STDERR "Label order: " . join(", ", @labels) . "\n";

my $first = 1;
my $expected_labels = "";
foreach (@labels) {
	if ($first) {
		$first = 0;
	} else {
		$expected_labels .= ", ";
	}
	$expected_labels .= "$_-x, $_-y, $_-z";
}

for ($first_frame .. $last_frame) {
	my @frame;
	my $line = <DATA>;
	while (1) {
		defined($line) || die "Ran out of data lines on frame $_.";
		chomp($line);
		if ($line eq $expected_labels) {
			print "Skipping a header row.\n";
			$line = <DATA>;
			next;
		}
		last;
	}
	@frame = split(/\s*,\s*/, $line);
	foreach (@frame) {
		$_ = $_ + 0;
	}
	#for (0 .. @slots-1) {
	#	push @frame, -1;
	#	push @frame, -1;
	#	push @frame, -1;
	#}
	#print "------------------\n";
	my $missing_count = 0;
	for (0 .. $num_3d_points-1) {
		my $x = $frame[$slots[$_]*3+0];
		my $y = $frame[$slots[$_]*3+1];
		my $z = $frame[$slots[$_]*3+2];
		my $info = 0;
		#TODO: ?? Handle missing values ??
		if ($x == 0 && $y == 0 && $z == 0) {
			print "WARNING: Don't have missing value handing just yet.";
		}
		
		my $point = pack("ffff", $x, $y, $z, $info);
		syswrite(OUTPUT, $point, 16) == 16 || die;
	}
}

close OUTPUT;
