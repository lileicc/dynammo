function [group, entrop, P, D, mu0, fp, component] = fingerprint_classify(X, varargin)                                             
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
%
% Out: 
%   group: label for each row of sequences
%   fp: the final fingerprints(features) 
%   entrop: the conditional entropy of prediction vs the groundtruth.
%   
% Example:
% cls = [1 1 2 2];
% [group, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 10,
% 'MaxIter', 100, 'Class', cls);

[P, D, mu0] = fingerprint(X, varargin{:});

HIDDEN = length(D);

ind = find(abs(imag(D)) > 1e-10); % assume 0==1e-10
num_real = HIDDEN;
% further eliminate the conjugate ones
if (~ isempty(ind)) 
    num_real = ind(1) - 1;
    %D = D(subind);
    Q = abs(P(:, (num_real+1):2:end)); %only the magnitude
end

%%%%%%%%%%%%%% *** best one %%%%%%%%%%%%%%%%%% write this in the paper
% PLiF
[component, fp] = princomp(Q, 'econ');
group = sign(fp(:,1));

class = varargin{find(strcmp('Class', varargin), 1) + 1};
cmpca = confusionmat(class, group);
cmpcah = condentropy(cmpca);
entrop = cmpcah;