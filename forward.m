function [mu, V, P, logli] = forward(X, model, varargin)
% The forward message passing for LDS (Kalman filtering)
% see learn_lds.m
% It estimates the P(z_n | x_1 ... x_n) ~ N(mu_n, V_n)
% 
%
% Args:
%   X is M * N, M is number of sequences, N is the time duration.
%   model: a struct with the following attributes:
%     A: transition matrix (also named as F in papers), H * H
%     C: transmission  matrix(also called projection matrix, named as G
% in papers), M * H
%     Q: transition covariance, H * H 
%     R: transmission covariance, M * M
%     mu0: initial states (also named as z_0 sometimes), H * 1
%     Q0: initial covariance, H * H
%
% Optional Args:
%   'Fast': use woodbery lemma for inverse of the matrix
%
% Returns:
%   mu is cell 1 * N, each with a matrix H * 1
%   V is cell 1 * N, each with a matrix H * H
%   P is cell 1 * N, each with a matrix H * H
%   logli: log-likelihood of the data under the current model. will
%

N = size(X, 2);
M = size(X, 1);
H = size(model.A, 1); %dimension of hidden variable
Ih = eye(H, H);

%predicted mean for hidden variable z
mu = cell(1, N);
V = cell(1, N);
P = cell(1, N);

% initialize
mu{1} = model.mu0;
V{1} = model.Q0;
logli = 0;

FAST = false;
a = strcmp('Fast', varargin);
if (any(a))
    FAST = true;
    invR = inv(model.R);
    invRC = invR * model.C;
    invCRC = model.C' * invRC;
end

LOGLI = true;
if (nargout < 4)
  LOGLI = false;
end

for i = 1:N
  if (i == 1)
    KP = model.Q0;    
    mu{i} =  model.mu0;
  else 
    P{i-1} = model.A * V{i-1} * model.A' + model.Q;
    KP = P{i-1};
    mu{i} =  model.A * mu{i-1};
  end
  if (FAST)
    invSig = invR - invRC / (inv(KP) + invCRC) * invRC';    
%     invSig = invR - invRC * ((inv(KP) + invCRC) \ invRC');    
  else
    sigma_c = model.C * KP * model.C' + model.R;
    invSig = inv(sigma_c);
    %invSig = inv(sigma_c);
  end
  K = KP * model.C' * invSig;
  u_c = model.C * mu{i};
  delta = X(:, i) - u_c;
  mu{i} = mu{i} + K * delta;
  V{i} = (Ih - K * model.C) * KP;
  if (LOGLI)
    posDef = delta' * invSig * delta / 2;
    if (posDef < 0)
      warning('det of not positive definite < 0');
    end
    logli = logli - M/2 * log(2 * pi) + logdet(invSig, 'chol') / 2 - posDef;
  end
end