%Puts data into a local coordinate frame.
%NOTE: data assumed to be in '..x,..y,..z,..x,..y,..z' interleaved format.
% (frames in rows)
function [data] = to_global(local, C, X, Y, Z)
	assert(mod(size(local,2),3) == 0);
	data = nan(size(local));
	for i = 1:3:size(local,2)
		coords = local(:,i:i+2);
		missing = (coords(:,1) == 0) | (coords(:,2) == 0) | (coords(:,3) == 0);
		data(:,i:i+2) = ...
			  X .* [coords(:,1) coords(:,1) coords(:,1)] ...
			+ Y .* [coords(:,2) coords(:,2) coords(:,2)] ...
			+ Z .* [coords(:,3) coords(:,3) coords(:,3)] ...
			+ C;
		data(missing,i+0) = 0;
		data(missing,i+1) = 0;
		data(missing,i+2) = 0;
	end
	assert(sum(sum(isnan(data))) == 0);
end
