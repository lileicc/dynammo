function [feature, pred, model_clds] = clds_cluster(X, varargin)
%clds_cluster extract features using CLDS model, and cluster the 
% sequences with clds features using k-means.
%
% Please see more details in 
%   Lei Li and B. Aditya Prakash (2011), 
%   "Time Series Clustering: Complex is Simpler!", 
%   In Proceedings of the 28th International Conference on Machine learning.
%
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
%
% Optional Args:
%   "NumClusters": followed by a number indicating the number of clusters
%   (default = 2)
% Other Arguments are defined in learn_lds
%
% Returns:
%   feature: extracted feature matrix
%   pred: predicted labels
%   model_clds: the learned clds model
%
% $Author: leili $@cs.cmu.edu
% $Date: 2011-10-05 16:20:17 -0400 (Wed, 05 Oct 2011) $
% $Rev: 332 $


a = find(strcmp('NumClusters', varargin), 1);
if (isempty(a))
  k = 2;
else
  k = varargin{a+1};
end

[model_clds, LL] = learn_clds(X, varargin{:});
features_clds = abs(model_clds.C);
[coeff1, feature] = princomp(features_clds(:, 1:end), 'econ');
pred = kmeans(feature, k, 'replicates', 10);
ce_clds = condentropy(trueclass, pred);

disp('CLDS conditional entropy');
disp(ce_clds);
