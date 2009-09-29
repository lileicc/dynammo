%Puts data into a local coordinate frame.
%NOTE: data assumed to be in '..x,..y,..z,..x,..y,..z' interleaved format.
% (frames in rows)
function [local] = to_local(data, C, X, Y, Z)
	assert(mod(size(data,2),3) == 0);
	local = nan(size(data));
	for i = 1:3:size(data,2)
		coords = data(:,i:i+2);
		missing = (coords(:,1) == 0) | (coords(:,2) == 0) | (coords(:,3) == 0);
		coords = coords - C;
		local(:,i+0) = dot(coords',X')';
		local(:,i+1) = dot(coords',Y')';
		local(:,i+2) = dot(coords',Z')';
		local(missing,i+0) = 0;
		local(missing,i+1) = 0;
		local(missing,i+2) = 0;
	end
	assert(sum(sum(isnan(local))) == 0);
end
