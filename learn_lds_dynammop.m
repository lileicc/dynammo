function [model, Xhat, LL] = learn_lds_dynammop(X, varargin)
% learning model parameters for Linear Dynamical Systems (LDS), also known
% as Kalman Filters. 
% Recover missing values using DynaMMo+ algorithm. 
% an improved version of the algorithm in 
% Lei Li, Jim McCann, Nancy Pollard, Christos Faloutsos. DynaMMo: Mining 
% and Summarization of Coevolving Sequences with Missing Values. 
% KDD '09, Paris, France.
%
% Linear Dynamical Systems are described by the following equations:
% z_1 = mu0 + w_1
% z_n = A z_{n-1} + w_n
% x_n = C x_n + v_n
% w_1 ~ N(0, Q0)
% w_n ~ N(0, Q)
% v_n ~ N(0, R)
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
%   'Observed', followed by a matrix with the same size as X, with binary
%   values denoting whether X(i,j) is observed(1) or missing (0). If not
%   provided, will automatically treat 0's in X as missing values.
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
%  diagonal Q0, Q, R.
%   'Plotfun', followed by a function handle (should take in X), for plotting purpose.
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
% derived from old function 
% function [A, Gamma, C, Sigma, u0, V0, LL] = learn_kalman(x, H, maxIter)

X_original = X;
N = size(X, 2);
M = size(X, 1);

% get number of hidden variables
a = find(strcmp('MaxIter', varargin));
if (isempty(a))
  maxIter = 10;
else
  maxIter = varargin{a+1};
end

% get number of hidden variables
a = find(strcmp('Hidden', varargin));
if (isempty(a))
  H = M;
else
  H = varargin{a+1};
end

% get the initial model
a = find(strcmp('Model', varargin));
if (isempty(a))
  model.A = eye(H, H) + randn(H, H);
  model.C = eye(M, H) + randn(M, H);
  model.Q = eye(H, H);
  model.R = eye(M, M);
  model.mu0 = randn(H, 1);
  model.Q0 = model.Q;
else
  model = varargin{a+1};
end

% get the observed
a = find(strcmp('Observed', varargin));
if (isempty(a))
  observed = (abs(X) > eps);
else
  observed = varargin{a+1};
end

% get plot function 
a = find(strcmp('Plotfun', varargin));
if (~isempty(a))
  plotFun = varargin{a+1};
end

CONV_BOUND = 1e-5;

ratio = 1;
diff = 1;
iter = 0;
oldLogli = -inf;

while ((ratio > CONV_BOUND || diff > CONV_BOUND) && (iter < maxIter) && (~ (isTiny(model.Q0) || isTiny(model.Q) || isTiny(model.R))))
  oldmodel = model;
  iter = iter + 1;
  if (iter > 1)
    [mu, V, P, logli, X] = forward_fly(X, model, observed);
  else
    [mu, V, P, logli] = forward(X, model);
  end
  [Ez, Ezz, Ez1z] = backward(mu, V, P, model);  
  Y = estimate_missing(X, Ez, model, observed);
  X(~observed) = Y(~observed);
  model = MLE_lds(X, Ez, Ezz, Ez1z, varargin{:});  
  logli = real(logli);
  diff = (logli - oldLogli);
  if (logli < oldLogli)
    warning('Loglikelihood decreases!');
  end
  ratio = abs(diff/logli) ;
  LL(iter) = logli;
  oldLogli = logli;
  fprintf('iteration = %d, logli = %d\n', iter, logli);
  if (exist(plotFun))
    plotFun(X);
  end
end
model = oldmodel;
Xhat = X;

function [t] = isTiny(sigma)
% test whether the matrix sigma is close to zero
%
t = (norm(sigma, 1) < eps) || (any(diag(sigma) < eps));