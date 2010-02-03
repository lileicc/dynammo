function [X, Y, model, W, bone, bone_var, time] = run_c3d_occlusion_bone(c3dcsv, missing_bone, missing_frame_start, missing_frame_end, varargin)
% use lds with bone constraints to learn the dynamics and recover the
% occlusion
% Args:
% c3dcsv: filename for original c3d data in csv format
% missing_bone: id of the bone made missing
% missing_frame_start, missing_frame_end: start and end frame for missing
%
% Returns:
% X: original data matrixx
% Y: recovered matrix
% A, ... V0: Kalman parameter
% bone: bone estimated from data (from get_bones)
% bone_var: bone variance
%
% The program will generate the following output file:
% a .mat file: contains everything after learning
% a .csv file: contains reconstructed motion in csv format
% a .avi file: the animation of recovered motion 
% a .eps file: hardcoded plot of relavent frames.

lab = importdata(c3dcsv);

% put to body local coordinate
[frame_C, frame_X, frame_Y, frame_Z] = fit_global(lab);
X_large = to_local(lab.data, frame_C, frame_X, frame_Y, frame_Z);
X_large = X_large'; % make it M * N
X = X_large / 1000;
N = size(X, 2);
M = size(X, 1);
[bone, bone_var] = get_bones(X, varargin{:});
fprintf('estimated %d bones\n', size(bone, 1));
W = true(M/3, N);
W(missing_bone, missing_frame_start:missing_frame_end) = false;
%figure;
%plot(X(:, 97:99));

tic; 
% setup the parameters
H = 16;
% maxIter = 10000;
% learning the missing value using on_the_fly_and_bone_constraints
if (isempty(find(strcmp('Newton', varargin), 1)))
  [model, Y, LL] = learn_lds_dynammop_bone(X, 'Hidden', H, 'Observed', W, 'Bone', bone, varargin{:});
  filename = sprintf('%s_%s_%d-%d_random_bone_fly', c3dcsv, num2str(missing_bone), missing_frame_start, missing_frame_end);
else
  [model, Y, LL] = learn_lds_dynammop_bone_newton(X, 'Hidden', H, 'Observed', W, 'Bone', bone, varargin{:});
  filename = sprintf('%s_%s_%d-%d_multi_bone_fly', c3dcsv, num2str(missing_bone), missing_frame_start, missing_frame_end);  
end
time = toc;

%% save the data
% play_mocap_skel(X, lab.colheaders);
Y_large = Y * 1000;
save(sprintf('%s.mat', filename));
csvwrite(sprintf('%s.csv', filename), Y_large);

%% play the recovered motion animation
play_mocap_skel(Y, lab.colheaders, sprintf('%s.avi', filename));

%% plot signal and bone lengths
i = missing_bone(1);
figure;
subplot(4,1,1);
plot(X((i*3 -2):(i*3), :)');
title('original');
subplot(4,1,2);
plot(Y((i*3 -2):(i*3), :)');
title('reconstructed');
subplot(4,1,3);
dx = sqrt(sum((X((37*3 - 2) : (37*3), :) - X((i*3 -2):(i*3), :)).^2));
dy = sqrt(sum((Y((37*3 - 2) : (37*3), :) - Y((i*3 -2):(i*3), :)).^2));
plot(dx);
%mx = max(dy) * 1.01;
%ylim([0, mx]);
title('bone length original');
subplot(4,1,4);
plot(dy);
%ylim([0, mx]);
title('bone length reconstructed');
saveas(gcf, sprintf('%s.eps', filename), 'psc2');

