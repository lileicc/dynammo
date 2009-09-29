function [error_pca_all, ratio_pca_all, h_pca_all] = compress_pca(X, varargin)
% compress using the pca method
% run through all possible H, get all error and ratio.
% Args:
%   X: M * N, M is number of sequences, N is number of time ticks
% Optional Args:
%   'Hidden', followed by a number denoting the number of hidden dimension.
%   if nothing provided, will run through all possibles.

X_m = mean(X, 2);
[Coeff, Score] = princomp(X', 'econ');
TOTALVAR = norm(X - repmat(X_m, 1, size(X, 2)), 'fro'); 
a = find(strcmp('Hidden', varargin), 1);
if (isempty(a))
  if (size(X, 1) < 20)
    cands = 1:size(X, 1);
  else
    cands = [1:20, 25: 5 :size(X, 1)];
  end
else
  cands = [varargin{a+1}];
end
h_pca_all = [];
error_pca_all = [];
ratio_pca_all = [];
for HIDDEN = cands
    error_pca = norm(Coeff(:, 1:HIDDEN) * Score(:, 1:HIDDEN)' + repmat(X_m, 1, size(X, 2)) - X, 'fro') / TOTALVAR;
    ratio_pca = (numel(X) + 2) / (size(X, 2) * (HIDDEN) + size(X, 1) * (HIDDEN + 1) + 3);
    error_pca_all = [error_pca_all, error_pca];
    ratio_pca_all = [ratio_pca_all, ratio_pca];
    h_pca_all = [h_pca_all, HIDDEN];
end     
