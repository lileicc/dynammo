function [model] = MLE_clds(X, Ez, Ezz, Ez1z, varargin)
% Maximum likelihood estimation of linear dymaical system (M-step)
%
% Args:
%   X: M * N, M is number of sequences, N is the time duration.
%   Ez, Ezz, Ezz1 as returned by forward-backward.
%
% Optional Args:
%  transition matrix option:
%   'DiagA', diagonal A
%   'FullA', full A
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
%    'model.A', if provided, the model.A will not be the given value
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
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
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

a = find(strcmp('model.mu0', varargin), 1);
if (~isempty(a))
  model.mu0 = varargin{a+1};
else
  model.mu0 = Ez{1};
end
a = find(strcmp('model.Q0', varargin), 1);
if (~isempty(a))
  model.Q0 = varargin{a+1};
else
  model.Q0 = Ezz{1} - Ez{1} * Ez{1}';
  if (any(strcmp('DiagQ0', varargin)))
    model.Q0 = diag(real(diag(model.Q0)));
  elseif (any(strcmp('FullQ', varargin)))
    % do nothing
    %diag(model.Q0) = real(diag(model.Q0));
  else
    model.Q0 = diag(repmat(real(trace(model.Q0) / H), H, 1));
  end
end


TAG = true; 
% full rank version
a = find(strcmp('model.A', varargin), 1);
if (isempty(a))  
  model.A = Sz1z / SzzN;
  %FIXA = false;
else
  model.A = varargin{a+1};
  %FIXA = true;
  TAG = false; % no need to do iteration
end

a = find(strcmp('model.Q', varargin), 1);
if (~isempty(a))
  model.Q = varargin{a+1};
else
  if (any(strcmp('DiagQ', varargin)))
    model.Q = real(diag((diag(Szz) - diag(Ezz{1}) - 2 * diag(model.A * Sz1z') + diag(model.A * SzzN * model.A')) / (N-1)));
  elseif (any(strcmp('FullQ', varargin)))
    tmp = model.A * Sz1z';
    model.Q = (Szz - Ezz{1} - tmp - tmp' + model.A * SzzN * model.A') / (N-1);
    %diag(model.Q) = real(diag(model.Q));
  else
    delta = (trace(Szz) - trace(Ezz{1}) - 2 * trace(model.A * Sz1z') + trace(model.A * SzzN * model.A')) / (N-1) / H;
    model.Q = diag(repmat(real(delta), H, 1));
  end
end

if (TAG)
iter = 1;
while iter < 3
  invQ = inv(model.Q);
  model.A = diag( sum((invQ .* (SzzN.')) \ (invQ .* (Sz1z).'), 2) );
  %tmp = model.A * Sz1z';
  %model.Q = (Szz - Ezz{1} - tmp - tmp' + model.A * SzzN * model.A') / (N-1);
  model.Q = real(diag((diag(Szz) - diag(Ezz{1}) - 2 * diag(model.A * Sz1z') + diag(model.A * SzzN * model.A')) / (N-1)));
  iter = iter + 1;
end
end


model.C = Sxz / Szz;

if (any(strcmp('DiagR', varargin)))
  model.R = diag(real((diag(X * X') - 2 * diag(model.C * Sxz') + diag(model.C * Szz * model.C')) / N));
elseif (any(strcmp('FullR', varargin)))
  tmp = model.C * Sxz';
  model.R = (X * X' - tmp - tmp' + model.C * Szz * model.C') / N;
else
  delta = (trace(X * X') - 2 * trace(model.C * Sxz') + trace(model.C * Szz * model.C')) / N / M;
  model.R = diag(repmat(real(delta), M, 1));
end
