function [model] = MLE_lds(X, Ez, Ezz, Ez1z, varargin)
%Maximum likelihood estimation of linear dymaical system (M-step)
%
% Args:
%   X: M * N, M is number of sequences, N is the time duration.
%   Ez, Ezz, Ez1z as returned by forward-backward.
%
% Optional Args:
%  covariance options:
%   'DiagQ0', if presented, will learn a diagonal covariance Q0
%   'DiagQ', if presented, will learn a diagonal covariance Q
%   'DiagR', if presented, will learn a diagonal covariance R
%   'IsotropicQ0', if presented, will learn a diagonal and isotropic Q0
%   'IsotropicQ', if presented, will learn a diagonal and isotropic Q
%   'IsotropicR', if presented, will learn a diagonal and isotropic R
%   'FullQ0', if presented, will learn any possible Q0
%   'FullQ', if presented, will learn any possible Q
%   'FullR', if presented, will learn any possible R
%  Note these options could not coexist for the same covariance matrix.
%  Default (no args given) the algorithm will learn with diagonal Q0, Q, R.
%
% Returns:
%   model: a struct with the following attributes:
%     A: transition matrix (also named as F in papers), H * H
%     C: transmission  matrix(also called projection matrix, named as G
% in papers), M * H
%     Q: transition covariance, H * H
%     R: transmission covariance, M * M
%     mu0: initial states (also named as z_0 sometimes), H * 1
%     Q0: initial covariance, H * H
%

N = size(X, 2);
M = size(X, 1);
H = size(Ez{1}, 1);
Sz1z = zeros(H, H);
Szz = zeros(H, H);
Sxz = zeros(M, H);

for i = 1:(N-1)
  Sz1z = Sz1z + Ez1z{i};
end

for i = 1:N
  Szz = Szz + Ezz{i};
  Sxz = Sxz + X(:, i) * Ez{i}';
end

SzzN = Szz - Ezz{N}; % sum of E[z, z] from 1 to n-1

model.mu0 = Ez{1};
model.Q0 = Ezz{1} - Ez{1} * Ez{1}';
if (any(strcmp('IsotropicQ0', varargin)))  
  model.Q0 = diag(repmat(trace(model.Q0) / H, H, 1));
elseif (any(strcmp('FullQ', varargin)))
  % do nothing
else
  model.Q0 = diag(diag(model.Q0));
end

model.A = Sz1z / SzzN;
if (any(strcmp('IsotropicQ', varargin)))
  delta = (trace(Szz) - trace(Ezz{1}) - 2 * trace(A * Sz1z') + trace(A * SzzN * A')) / (N-1) / H;
  model.Q = diag(repmat(delta, H, 1));
elseif (any(strcmp('FullQ', varargin)))
  tmp = model.A * Sz1z';
  model.Q = (Szz - Ezz{1} - tmp - tmp' + model.A * SzzN * model.A') / (N-1);
else
  model.Q = diag((diag(Szz) - diag(Ezz{1}) - 2 * diag(model.A * Sz1z') + diag(model.A * SzzN * model.A')) / (N-1));
end

model.C = Sxz / Szz;

if (any(strcmp('IsotropicR', varargin)))
  delta = (trace(X' * X) - 2 * trace(model.C * Sxz') + trace(model.C * Szz * model.C')) / N / M;
  model.R = diag(repmat(delta, M, 1));
elseif (any(strcmp('FullR', varargin)))
  tmp = model.C * Sxz';
  model.R = (X' * X - tmp - tmp' + model.C * Szz * model.C') / N;
else
  model.R = diag((diag(X * X') - 2 * diag(model.C * Sxz') + diag(model.C * Szz * model.C')) / N);
end