function [V] = genPSD(M, dia)
% generate PSD matrix
R = triu(complex(randn(M, M), randn(M, M)));
V = R' * R;
if ((nargin > 1) && (dia))
  V = diag(diag(V));
end
end
