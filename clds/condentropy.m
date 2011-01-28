function [h] = condentropy(a, b)
%CONDENTROPY  Estimates the conditional entropy from confusion matrix
%          H(column|row)
%
% [h] = condentropy(a) computes the conditional entropy using a as
% confusion matrix
%
% [h] = condentropy(a, b) first computes the confusion matrix from two vectors
% a, b, and then computes the conditional entropy 
%   
% $Author$@cs.cmu.edu
% $Date$
% $Rev$

if (nargin > 1)
  cmat = confusionmat(a, b);
else
  cmat = a;
end
s = sum(cmat, 2);
total = sum(s);
s(s==0) = 1;
if (total == 0) total = 1; end
condp = cmat ./ repmat(s, 1, size(cmat, 2));
condp(cmat == 0) = 1;
jointp = cmat / total;
h = - sum(sum(jointp .* log(condp)));