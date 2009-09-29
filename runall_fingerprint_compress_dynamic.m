function [error_finger_all, ratio_finger_all, h_finger_all] = runall_fingerprint_compress_dynamic(X, varargin)
% try many different hidden dimensions and get the compression ratio and
% errors using fingerprint compression with dynamic programming.
% the resulting ratio is sorted descently, and the error is drawn from
% skyline shape.
% see fingerprint_compress_dynamic.m for additional arguments

error_finger_all = [];
ratio_finger_all = [];
h_finger_all = [];
cands = [1 : 4, 5:5:size(X, 1)];
for HIDDEN = cands  
  [error_f, ratio_f] = fingerprint_compress_dynamic(X, 'Hidden', HIDDEN, varargin{:});
  apsize = length(error_f);
  error_finger_all = [error_finger_all, error_f];
  ratio_finger_all = [ratio_finger_all, ratio_f];
  h_finger_all = [h_finger_all, repmat(HIDDEN, 1, apsize)];
  fprintf('Hidden = %d, error = %d, ratio = %d\n', HIDDEN, error_f(1), ratio_f(1));
end
[ratio_finger_all, tmpIdx] = sort(ratio_finger_all, 'descent');
error_finger_all = cummin(error_finger_all(tmpIdx)); % get the skyline plot
