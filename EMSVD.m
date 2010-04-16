function [Y, U, S, V, delta] = EMSVD(X, W, rank, MAX_ITER)
%estimate the Y = U S V' to approxmate X
%X: N by M 
%W: N by M to indicating which of X is observed/missing
%rank: the rank of U, S, V

if (nargin < 4)
  MAX_ITER = 100;
end

N = size(X, 1);
M = size(X, 2);
if (nargin < 3)
  rank = min(M,N);
end

%for i=1:M
%  X(W(:,i)==0,i) = mean(X(W(:,i)==1,i));
%end
Y = linear_interp(X, W);
X(W==0) = Y(W==0);
%delta = norm(Y - X, 'fro');

iter = 1;
delta = [];

CONV_BOUND = 1e-6;
ratio = 1;
delta_old = 0;

while ((ratio > CONV_BOUND) && (iter < MAX_ITER))
  [U, S, V] = svd(X,0);
  Y = U(:,1:rank) * S(1:rank,1:rank) * V(:, 1:rank)';
  delta_new = norm(Y - X, 'fro');
  ratio = abs((delta_new - delta_old) / delta_new);
  delta = [delta; delta_new];
  delta_old = delta_new;
  iter = iter + 1;
  X(W==0) = Y(W==0);
end
