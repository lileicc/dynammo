%Fit a global frame to some mocap loaded from a csv (so provide the entire file)
%Return is frame center, x,y,z vectors per frame (in rows)
function [C, X, Y, Z] = fit_global(csvdata)
	x_weights = zeros(size(csvdata.data));
	y_weights = zeros(size(csvdata.data));
	z_weights = zeros(size(csvdata.data));
	front_sum = zeros(size(csvdata.data,1),3);
	front_weights = zeros(size(csvdata.data,1),3);
	back_sum = zeros(size(csvdata.data,1),3);
	back_weights = zeros(size(csvdata.data,1),3);
	for c = 1:size(csvdata.colheaders,2)
		s = char(csvdata.colheaders(c));
		fw = 0; %weight for front sum
		bw = 0; %weight for back sum
		cs = s(end-4:end-2);
		if (strcmp(cs, ':T8') ...
		 || strcmp(cs, ':C7'))
			bw = 1;
        else		
            cs = s(end-5:end-2);
            if (strcmp(cs, ':T10'))
                bw = 1;
            else
                cs = s(end-6:end-2);
        		if (strcmp(cs, ':RBAC') ...
            	 || strcmp(cs, ':LBWT') ...
                 || strcmp(cs, ':RBWT'))
                    bw = 1;
                else                 
                    if (strcmp(cs, ':STRN') ...
                     || strcmp(cs, ':RFWT') ...
                     || strcmp(cs, ':LFWT') ...
                	 || strcmp(cs, ':CLAV'))
                    	fw = 1;
                    else                    
                        cs = s(end-9:end-2);
                        if (strcmp(cs, ':NEWLBAC') ...
                        || strcmp(cs, ':NEWRBAC'))
                    	bw = 1;
                        end
                    end
                end
            end
        end
    	if (fw + bw == 0)
        	continue;
        end
		cs = s(end-1:end);
		dim = 0;
		if (strcmp(cs, '-x'))
			dim = 1;
		elseif (strcmp(cs, '-y'))
			dim = 2;
		elseif (strcmp(cs, '-z'))
			dim = 3;
		else
			assert(false);
		end
		weights = csvdata.data(:,c) ~= 0;
		%Introduce a linear fade in and out on markers that are dropped:
		weights = min(bwdist(~weights, 'cityblock')/10,1);
		front_sum(:,dim) = front_sum(:,dim) + fw * (weights .* csvdata.data(:,c));
		front_weights(:,dim) = front_weights(:,dim) + fw * weights;
		back_sum(:,dim) = back_sum(:,dim) + bw * (weights .* csvdata.data(:,c));
		back_weights(:,dim) = back_weights(:,dim) + bw * weights;
	end
	frame = front_weights;
	disp(sprintf('Using, on average, %f front and %f back markers per frame', ...
		sum(sum(front_weights)) / size(front_weights, 1) / 3, ...
		sum(sum(back_weights)) / size(back_weights, 1) / 3) );
	assert(sum(sum(front_weights == 0)) == 0);
	assert(sum(sum(back_weights == 0)) == 0);
	front = front_sum ./ front_weights;
	back = back_sum ./ back_weights;
	%plot3(front(:,1), front(:,2), front(:,3), '-r');
	%hold on;
	%plot3(back(:,1), back(:,2), back(:,3), '-b');
	%hold off;
	
	%-------------------
	%Ok, make frame by normalizing back->front vector and projecting to floor.
	C = [back(:,1:2), zeros(size(back,1),1)];
	d = front(:,1:2) - back(:,1:2);
	l = sqrt(dot(d',d'))';
	d = d ./ [l l];
	X = [d(:,1) d(:,2) zeros(size(back,1),1)];
	Y = [-d(:,2) d(:,1) zeros(size(back,1),1)];
	Z = [zeros(size(back,1),2) ones(size(back,1),1)];
end
