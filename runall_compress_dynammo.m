function [error_dynammo_all, ratio_dynammo_all, h_dynammo_all] = runall_compress_dynammo(X, varargin)
% try many different hidden dimensions and get the compression ratio and
% errors using dynammo compression with dynamic programming.
% the resulting ratio is sorted descently, and the error is drawn from
% skyline shape.
% see compress_dynammo.m for additional arguments

error_dynammo_all = [];
ratio_dynammo_all = [];
h_dynammo_all = [];
cands = [1 : 4, 5:5:size(X, 1)];
for HIDDEN = cands  
  [error, ratio] = compress_dynammo(X, 'Hidden', HIDDEN, varargin{:});
  apsize = length(error);
  error_dynammo_all = [error_dynammo_all, error];
  ratio_dynammo_all = [ratio_dynammo_all, ratio];
  h_dynammo_all = [h_dynammo_all, repmat(HIDDEN, 1, apsize)];
  fprintf('Hidden = %d, error = %d, ratio = %d\n', HIDDEN, error(1), ratio(1));
end
[ratio_dynammo_all, tmpIdx] = sort(ratio_dynammo_all, 'descend');
error_dynammo_all = cummin(error_dynammo_all(tmpIdx)); % get the skyline plot
