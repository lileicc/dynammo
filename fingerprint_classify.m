function [group, entrop, P, D, mu0, coordinate, component] = fingerprint_classify(X, varargin)                                             
% kalman fingerprinting clustering
%  (PLF method)
%   X: M * N matrix, M is number of sequences, N is the time duration.
% group is a vector size of M, telling the clusters
% entropy is the conditional entropy of prediction vs the groundtruth. 
% 
% example:
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
% PLF
[component, coordinate] = princomp(Q, 'econ');
group = sign(coordinate(:,1));

class = varargin{find(strcmp('Class', varargin), 1) + 1};
cmpca = confusionmat(class, group);
cmpcah = condentropy(cmpca);
entrop = cmpcah;