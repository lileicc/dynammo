function play_mocap_skel_length(data, colheaders, bones, filename, bound_in)
% a function to play the motion sequence 
% by drawing body skeletons
% data is a matrix with each column corresponding to time tick 
% colheaders is the name for each marker
% filename: optional file for saving the movie

if (nargin < 4)
  filename = 'temp.avi';
end

bounds = reshape([min(reshape(min(data'), 3, size(data,1)/3)');max(reshape(max(data'), 3, size(data,1)/3)')], 6, 1);

aviobj = avifile(filename, 'fps', 120, 'compression', 'None');
N = size(data, 2);
f = figure;
set(f, 'RendererMode', 'manual');
set(f, 'Renderer', 'opengl'); %use opengl or zbuffer

%Figure out the largest bone length distortion:
bound = 0;
for f = 1:N
	markers = reshape(data(:,f), 3, prod(size(data(:,f)))/3)';
	for i = 1:size(bones,1)
		a = bones(i,1);
		b = bones(i,2);
		d = markers(a,:)-markers(b,:);
		d = d * d';
		d = sqrt(d) - bones(i,3);
		bound = max([bound,abs(d)]);
	end
end

disp(sprintf('Largest deviation from bone length is %f',bound));
if (nargin >= 5)
	bound = bound_in;
end
%bound = 0.04; %better to see length changes.
disp(sprintf('Setting scale max at %f',bound));

%Actually draw frames:

for i = 1:N
  draw_skel_length(data(:, i), colheaders, bones, bound);
  axis(bounds);
  view(53, 34); %seems better than viewing from the back
  hold on;
  title(strcat('Time = ', num2str(i), '/', num2str(N)));
  drawnow;
  hold off;
  F = getframe(gcf);
  aviobj = addframe(aviobj, F);
  %pause(1/120);
  %pause;
end
aviobj = close(aviobj);


