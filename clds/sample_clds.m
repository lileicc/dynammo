function [ X, model ] = sample_clds( arg1, arg2, arg3 )
%generate samples for a given model (or random generate a model)
% two example usage:
% sample_clds(N, model)
% OR:
% sample_clds(N, M, H)
%

if (nargin == 3)
  N = arg1;
  M = arg2;
  H = arg3;
  model.A = diag(exp(1i * randn(H, 1)));
  model.C = eye(M, H) + complex(randn(M, H), randn(M, H));
  model.Q = genPSD(H);
  model.R = genPSD(M);
  model.mu0 = complex(randn(H, 1), randn(H, 1));
  model.Q0 = genPSD(H);  
elseif (nargin == 2)
  model = arg2;
  N = arg1;
  M = size(model.C, 1);
  H = size(model.C, 2);
end

revQ0 = chol(model.Q0, 'lower');
revQ = chol(model.Q, 'lower');
revR = chol(model.R, 'lower');

X = complex(zeros(M, N), zeros(M, N));
mu = model.mu0;
for k = 1 : N
  revmu0 = revQ0 \ mu;
  z = revQ0 * cnormrnd(revmu0);
  mux = model.C * z;
  X(:, k) = revR * cnormrnd(revR \ mux);
  muz = model.A * z;
  mu = revQ * cnormrnd(revQ \ muz);
end

end

function [x] = cnormrnd(mu)
a = real(mu);
[n,m] = size(a);
z = zeros(n, m);
s = repmat(0.5, n, m);
x = mu + complex(normrnd(z, s), normrnd(z, s));
end
