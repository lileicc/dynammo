function [ model, LL ] = learn_clds(X, varargin)
%learn_clds learning a complex dynamical systems
% Linear Dynamical Systems are described by the following equations:
% z_1 = mu0 + w_1
% z_n = A z_{n-1} + w_n
% x_n = C x_n + v_n
% w_1 ~ CN(0, Q0)
% w_n ~ CN(0, Q)
% v_n ~ CN(0, R)
% here A can be diagonal
%
% Please see more details in 
%   Lei Li and B. Aditya Prakash (2011), 
%   "Time Series Clustering: Complex is Simpler!", 
%   In Proceedings of the 28th International Conference on Machine learning.
% 
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
%
% Optional Args:
%   'Hidden', followed by an integer indicating the number of hidden
%   variables.
%   'MaxIter', followed by an integer indicating the number of max
%   iterations.
%   'Model', followed by a struct denoting a model to start with (see
%   below).
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
%  Default (no args given) the algorithm will learn with H=M, MaxIter=10,
%  diagonal and isotropic Q0, Q, R.
%    'model.A', if provided, the model.A will remain the given value
%    'model.mu0', if provided, the model.mu0 will remain the given value
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
% Example:
% t = 1:100;
% x1 = sin(2 * pi * t / 50);
% x2 = sin(2 * pi * t / 50 + pi / 4);
% X = [x1; x2];
% model = learn_lds(X, 'Hidden', 2, 'MaxIter', 100);
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
%
% change log:
%
%

N = size(X, 2);
M = size(X, 1);

% get number of hidden variables
a = find(strcmp('MaxIter', varargin), 1);
if (isempty(a))
  maxIter = 10;
else
  maxIter = varargin{a+1};
end

% get number of hidden variables
a = find(strcmp('Hidden', varargin), 1);
if (isempty(a))
  H = M;
else
  H = varargin{a+1};
end

% get the initial model
a = find(strcmp('Model', varargin), 1);
if (isempty(a))
  model.A = diag(exp(1i * randn(H, 1)));
  %model.A = diag(complex(randn(H, 1), randn(H, 1)));
  %model.C = ones(M, H) + complex(randn(M, H), randn(M, H));
  model.C = zeros(M, H);
  model.Q = eye(H, H) * 1E-2;
  model.R = eye(M, M) * 1E-2;
  model.mu0 = ones(H, 1) + complex(randn(H, 1), randn(H, 1));
  model.Q0 = model.Q;
else
  model = varargin{a+1};
end

a = find(strcmp('model.A', varargin), 1);
if (~isempty(a))
  model.A = varargin{a+1};
end
a = find(strcmp('model.mu0', varargin), 1);
if (~isempty(a))
  model.mu0 = varargin{a+1};
end
a = find(strcmp('model.Q0', varargin), 1);
if (~isempty(a))
  model.Q0 = varargin{a+1};
end
a = find(strcmp('model.Q', varargin), 1);
if (~isempty(a))
  model.Q = varargin{a+1};
end
a = find(strcmp('model.R', varargin), 1);
if (~isempty(a))
  model.R = varargin{a+1};
end
a = find(strcmp('model.C', varargin), 1);
if (~isempty(a))
  model.C = varargin{a+1};
end

LOGLI = true;
if (nargout < 2)
  LOGLI = false;
end

CONV_BOUND = 1e-5;

ratio = 1;
diff = 1;
iter = 0;
oldLogli = -inf;

while ((ratio > CONV_BOUND || diff > CONV_BOUND) && (iter < maxIter) && (~ (isTiny(model.Q0) || isTiny(model.Q) || isTiny(model.R))))
  iter = iter + 1;
  if (LOGLI)
    [u, UU, P, logli] = forward(X, model, varargin{:});
    logli = real(logli);
    diff = (logli - oldLogli);
    if (logli < oldLogli)
      warning('Loglikelihood decreases!');
    end
    ratio = abs(diff/logli) ;
    LL(iter) = logli;
    oldLogli = logli;
    fprintf('iteration = %d, logli = %d\n', iter, logli);
  else 
    [u, UU, P] = forward(X, model, varargin{:});
    fprintf('iteration = %d\n', iter);
  end  
  [Ez, Ezz, Ez1z] = backward(u, UU, P, model);
  oldmodel = model;  
  model = MLE_clds(X, Ez, Ezz, Ez1z, varargin{:});
end
model = oldmodel;
end
  
function [t] = isTiny(sigma)
% test whether the matrix sigma is close to zero
%
t = (norm(sigma, 1) < eps) || (any(diag(sigma) < eps));

end

