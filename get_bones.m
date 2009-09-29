function [bone, bone_var, bone_length_series] = get_bones(data, varargin)
% input the motion capture data in marker positions
% extract the bones using thresholding on variance of pairwise bone distance
% assume the position in 3d space.
% Args:
%   data: M * N matrix
% Optional Args:
%   'Threshold', followed by a number denoting the threshold used to
%   determine the bone
%   'SkeletonFile', followed by a .mat file denoting the human body
%   skeleton
% Returns:
%

N = size(data, 2);
M = size(data, 1);
k = M/3;
dist = zeros(k, k);
bone_length_series = zeros(N, k, k);
variance = zeros(k, k);
for i = 1:k
  for j = 1:k
    bone_length_series(:, i, j) = sqrt(sum((data((i*3 - 2) : (i*3), :) - data((j*3 - 2):(j*3), :)).^2));
    dist(i,j) = mean(bone_length_series(:, i, j));
    variance(i, j) = var(bone_length_series(:, i, j));
  end
end

bone_var = variance;

a = find(strcmp('Threshold', varargin));
if (~isempty(a))
  THRESHOLD = varargin{a+1};
  dist(variance > THRESHOLD) = 0;
  [x, y, d] = find(dist);
  %idx = x < y;
  %bone = [x(idx), y(idx), d(idx)];
  bone = [x, y, d];
else
  skeleton_file = 'skeleton.mat';
  load(skeleton_file);
  bone = zeros(size(skeleton, 1), 3);
  for i = 1 : size(bone, 1)
    bone(i, 1:2) = skeleton(i, 1:2);
    bone(i, 3) = dist(bone(i, 1), bone(i, 2));
  end
end



