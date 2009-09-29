function play_mocap_skel(data, colheaders, filename)
% a function to play the motion sequence 
% by drawing body skeletons
% data is a matrix with each column corresponding to time tick 
% colheaders is the name for each marker
% filename: optional file for saving the movie

if (nargin < 3)
  filename = 'temp.avi';
end

aviobj = avifile(filename, 'fps', 120);
N = size(data, 2);
figure;
for i = 1:N
  draw_skel(data(:, i), colheaders);
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


