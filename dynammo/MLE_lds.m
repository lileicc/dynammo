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
%  Default (no args given) the algorithm will learn with diagonal 
%  and isotropic Q0, Q, R.
%   'Fixed', followed by a model with some fieds of (A, C, Q, R, mu0, Q0), 
%           do not learn, use the provided parameters 
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
% $Author$@cs.berkeley.edu
% $Date$
% $Rev$
%

a = find(strcmp('Fixed', varargin), 1);
LEARNA = true;
LEARNQ = true;
LEARNC = true;
LEARNR = true;
LEARNMU0 = true;
LEARNQ0 = true;
if (~isempty(a))
  model = varargin{a+1};
  if (isfield(model, 'A'))
    LEARNA = false;
  end
  if (isfield(model, 'C'))
    LEARNC = false;
  end
  if (isfield(model, 'Q'))
    LEARNQ = false;
  end
  if (isfield(model, 'R'))
    LEARNR = false;
  end
  if (isfield(model, 'mu0'))
    LEARNMU0 = false;
  end
  if (isfield(model, 'Q0'))
    LEARNQ0 = false;
  end
end

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

if (LEARNMU0)
  model.mu0 = Ez{1};
end

if (LEARNQ0)
  model.Q0 = Ezz{1} - Ez{1} * Ez{1}';
  if (any(strcmp('DiagQ0', varargin)))  
    model.Q0 = diag(diag(model.Q0));
  elseif (any(strcmp('FullQ', varargin)))
    % do nothing
  else
    model.Q0 = diag(repmat(trace(model.Q0) / H, H, 1));
  end
end

if (LEARNA)
  model.A = Sz1z / SzzN;
end
if (LEARNQ)
  if (any(strcmp('DiagQ', varargin)))
    model.Q = diag((diag(Szz) - diag(Ezz{1}) - 2 * diag(model.A * Sz1z') + diag(model.A * SzzN * model.A')) / (N-1));
  elseif (any(strcmp('FullQ', varargin)))
    tmp = model.A * Sz1z';
    model.Q = (Szz - Ezz{1} - tmp - tmp' + model.A * SzzN * model.A') / (N-1);
  else
    delta = (trace(Szz) - trace(Ezz{1}) - 2 * trace(model.A * Sz1z') + trace(model.A * SzzN * model.A')) / (N-1) / H;
    model.Q = diag(repmat(delta, H, 1));
  end
end

if (LEARNC)
  model.C = Sxz / Szz;
end

if (LEARNR)
  if (any(strcmp('DiagR', varargin)))
    model.R = diag((diag(X * X') - 2 * diag(model.C * Sxz') + diag(model.C * Szz * model.C')) / N);
  elseif (any(strcmp('FullR', varargin)))
    tmp = model.C * Sxz';
    model.R = (X * X' - tmp - tmp' + model.C * Szz * model.C') / N;
  else
    delta = (trace(X * X') - 2 * trace(model.C * Sxz') + trace(model.C * Szz * model.C')) / N / M;
    model.R = diag(repmat(delta, M, 1));
  end
end