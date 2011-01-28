function [u, UU, P, logli] = forward(X, model, varargin)
% The forward message passing for CLDS (complex dynamical systems)
%
% It estimates the P(z_n | x_1 ... x_n) ~ N(mu_n, V_n)
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
% see also learn_slds.m
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
%

N = size(X, 2);
M = size(X, 1);
H = size(model.A, 1); %dimension of hidden variable
Ih = eye(H, H);
Im = eye(M, M);

%predicted mean for hidden variable z
u = cell(1, N);
UU = cell(1, N);
P = cell(1, N);

% initialize
%u{1} = model.mu0;
%UU{1} = model.Q0;
logli = 0;

% FAST = false;
% a = strcmp('Fast', varargin);
% if (any(a))
%     FAST = true;
%     invR = inv(model.R);
%     invRC = invR * model.C;
%     invCRC = model.C' * invRC;
% end

invR = inv(model.R);
invRC = invR * model.C;
invCRC = model.C' * invRC;

LOGLI = true;
if (nargout < 4)
  LOGLI = false;
end

for i = 1:N
  if (i == 1)
    %KP = model.Q0;
    P{i} = model.Q0;
    u{i} =  model.mu0;
  else 
    P{i} = model.A * UU{i-1} * model.A' + model.Q;
    %KP = P{i-1};
    u{i} =  model.A * u{i-1};
  end 
  %P{i}(Ih ~= 0) = abs(P{i}(Ih ~= 0));
  if (abs(imag(P{i}(Ih ~= 0))) > 1E-10) 
    warning('det of not positive definite < 0 @ forward propagation');
  end
  K = (inv(P{i}) + invCRC) \ invRC';
  %sigma_c = model.C * KP * model.C' + model.R;
  %invSig = inv(sigma_c);
  %K = P{i} * model.C' * invSig;
  u_c = model.C * u{i};
  delta = X(:, i) - u_c;
  u{i} = u{i} + K * delta;
  UU{i} = (Ih - K * model.C) * P{i};
  if (LOGLI)
    %sigma_c = model.C * P{i-1} * model.C' + model.Q;
    invSig = P{i} * model.C' \ K;
    invSig(Im ~= 0) = abs(invSig(Im ~= 0));
    posDef = real(delta' * invSig * delta);    
    if (posDef < 0)
      warning('det of not positive definite < 0');
    end
    logli = logli - M * log(pi) + logdet(invSig, 'chol') / 2 - posDef;
  end
end