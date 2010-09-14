function [group, fp, entrop, P, D, mu0, model] = fingerprint_classify(X, varargin)                                             
% kalman fingerprinting clustering
%  (PLiF method)
%
%  Lei Li, B. Aditya Prakash, Christos Faloutsos. 
%  Parsimonious Linear Fingerprinting for Time Series. 
%  VLDB 2010. Singapore
%
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
%   'Class': followed by a vector indicating the true class labels.
%   'Normalize': whether to normalize the data (default=false)
%
% Out: 
%   group: label for each row of sequences
%   fp: the final fingerprints(features) 
%   entrop: the conditional entropy of prediction vs the groundtruth. only
%   if the 'Class' argument is present. otherwise will just =0.
%   P: complex projection matrix
%   D: complex transition matrix
%   mu0: initial starting state
%   
% Example:
% cls = [1 1 2 2];
% [group, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 10,
% 'MaxIter', 100, 'Class', cls);
%
% See also: fingerprint, learn_lds
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$

if (any(strcmp('Normalize', varargin)))
  mn = mean(mean(X));
  st = sum(sum((X - mn).^2)) / (numel(X));
  X = (X - mn) ./ st;
end

[fp, P, D, mu0, zhat, model] = fingerprint(X, varargin{:});

%%%%%%%%%%%%%% *** best one %%%%%%%%%%%%%%%%%% write this in the paper
% PLiF
group = sign(fp(:,1));

a = find(strcmp('Class', varargin), 1);
if (isempty(a))
  entrop = 0;
else
  class = varargin{a + 1};
  cmpca = confusionmat(class, group);
  entrop = condentropy(cmpca);
end