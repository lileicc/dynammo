function [ Y ] = estimate_missing( X, Ez, model, observed )
%ESTIMATE_MISSING Summary of this function goes here
%   estimate the missing values in X
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
%

N = size(X, 2);
Y = X;
for i = 1 : N
  Y(:, i) = model.C * Ez{i};
end
end

