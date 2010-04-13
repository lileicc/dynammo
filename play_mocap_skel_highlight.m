function play_mocap_skel_highlight(data, colheaders, bones, missing, filename)
% a function to play the motion sequence 
% by drawing body skeletons
% data is a matrix with each column corresponding to time tick 
% colheaders is the name for each marker
% filename: optional file for saving the movie

if (nargin < 5)
  filename = '';
end

assert(size(missing,2) == size(data,2));
assert(size(missing,1) == size(data,1)/3);

bounds = reshape([min(reshape(min(data'), 3, size(data,1)/3)');max(reshape(max(data'), 3, size(data,1)/3)')], 6, 1);

if strcmp(filename, '')
	aviobj = nan();
	disp('Not recording.');
else
	aviobj = avifile(filename, 'fps', 120, 'compression', 'None');
	disp(strcat('Recording to ', filename));
end
N = size(data, 2);
f = figure;
set(f, 'RendererMode', 'manual');
set(f, 'Renderer', 'opengl'); %use opengl or zbuffer
set(f, 'Position', [0 0 500 700]);

%Actually draw frames:

freecolors = [ ...
  1.0 0   0  ; ...
  0   0.9 0  ; ...
  0.2 0.2 1  ; ...
  0.8 0.9 0  ; ...
  1.0 0   1.0; ...
  0   0.9 1  ; ...
  0   0.5 1.0; ...
  0.5 1.0 0  ; ...
];

evermissing = max(missing')';
bonecolors = zeros(size(bones,1), 3);
for b = 1:size(bones,1)
	if evermissing(bones(b,1)) || evermissing(bones(b,2))
		assert(size(freecolors,1) > 0);
		bonecolors(b,:) = freecolors(1,:);
		freecolors = freecolors(2:end,:);
	end
end

for i = 1:N
  framecolors = nan(size(bones));
  for b = 1:size(bones,1)
  	if missing(bones(b,1),i) || missing(bones(b,2),i)
		framecolors(b,:) = bonecolors(b,:);
	end
  end
  draw_skel_highlight(data(:, i), colheaders, bones, framecolors, missing(:,i));
  axis(bounds);
  view(53, 34); %seems better than viewing from the back
  hold on;
  title(sprintf('Frame %d / %d', i, N));
  drawnow;
  hold off;
  if ~isnan(aviobj)
	  F = getframe(gcf);
	  aviobj = addframe(aviobj, F);
  end
  %pause(1/120);
  %pause;
end
if ~isnan(aviobj)
	aviobj = close(aviobj);
end


