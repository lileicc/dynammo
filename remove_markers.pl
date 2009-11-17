#!/usr/bin/perl -w

#
# ./remove_markers.pl [name-or-column]:start:end ... < infile.csv > outfile.csv
#


$header = <STDIN>;

print $header;

chomp $header;
$header =~ s/^\s+//;
$header =~ s/\s+$//;
@names = split(/\s*,\s*/,$header);

$frame = 0;
while ($line = <STDIN>) {
	chomp $line;
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	@vals = split(/\s*,\s*/,$line);
	scalar(@vals) == scalar(@names) || die "Wrong number of values in line";
	my ($last_start, $last_end);
	for my $rspec (@ARGV) {
		($rspec =~ /(.+):(\d+|):(\d+|)/) || die "Invalid removal spec $rspec";
		my ($name, $start, $end) = ($1, $2, $3);

		if ($start eq "") { $start = $last_start; }
		if ($end eq "") { $end = $last_end; }

		$last_start = $start;
		$last_end = $end;
		if ($frame < $start || $frame > $end) {
			next;
		}
		for my $x (0 .. (@vals - 1)) {
			if ($names[$x] =~ /$name/) {
				$vals[$x] = 0;
			}
		}
	}
	print join(", ", @vals) . "\n";
	++$frame;
}
