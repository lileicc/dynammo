function [h] = condentropy(cmat)
%CONDENTROPY  Estimates the conditional entropy from confusion matrix
%          H(column|row)


s = sum(cmat, 2);
total = sum(s);
s(s==0) = 1;
if (total == 0) total = 1; end
condp = cmat ./ repmat(s, 1, size(cmat, 2));
condp(cmat == 0) = 1;
jointp = cmat / total;
h = - sum(sum(jointp .* log(condp)));