function [model, Xhat, LL] = learn_lds_dynammop_bone(X, varargin)
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
%   'Bone', followed by a matrix indicating the bone contraints
%      .. * 3 matrix: each row with 
%      <bone 1, bone 2, bone length>
%      assume bone is symmetric. i.e. <bone 2, bone 1, bone length> will
%      also be in bone matrix.
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
% [model, Xhat, LL] = learn_lds_dynammop(X, 'Hidden', 2, 'MaxIter', 100, 'PlotFun', @(X)plot(X'));
%
% derived from old function 
% function [A, Gamma, C, Sigma, u0, V0, LL] = learn_kalman(x, H, maxIter)
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
%

X_original = X;
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

% get the bone
a = find(strcmp('Bone', varargin), 1);
if (isempty(a))
  % no bone contraints
  bone = [];
else
  bone = varargin{a+1};
end

% get the initial model
a = find(strcmp('Model', varargin), 1);
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
a = find(strcmp('Observed', varargin), 1);
if (isempty(a))
  observed = ~isnan(X);
  if (~any(~observed))
    % if there is no NAN in the matrix,
    % try to regard the 0 as missing
    observed = (abs(X) > eps);
  end
else
  observed = varargin{a+1};
  %Dim = 3; % default
  % observed can be on bone joints or on marker coordinates
  
  Dim = int32(M / size(observed, 1));
  observed = (reshape(repmat(observed', 3, 1), N, M))';
end

% get plot function 
a = find(strcmp('PlotFun', varargin), 1);
if (~isempty(a))
  plotFun = varargin{a+1};
end

% use linear interpolation as an initialization
Y = linear_interp(X, observed);
X(~observed) = Y(~observed);

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
  
  % make the bone contraints
  if (((iter > 20) && (rem(floor(iter / 4), 25) ~= 0)) || (iter > 500))
    % do bone length constraint
    %[u, V, P, logli] = forward(X, A, Gamma, C, Sigma, u0, V0);
    
    %ucap is the estimated E[z_n]
    %[ucap, Vcap, J] = backward(u, V, P, A);
    sq = randperm(size(bone, 1));
    % estimate the missing values and lagrangian multiplier
    for t = 1 : N
      % for i = 1 : size(bone, 1)
      % use random optimization order
      for i = sq
        %for i = 1:size(bone, 1)
        k = bone(i, 1);
        j = bone(i, 2);
        %                 if (~W(t, k) && W(t, j))
        if (~observed(k * Dim, t))
          %[y, eta] = estimate_eta(Sigma, k, j, bone(i, 3), X(t, :)', C, ucap{t}, Dim);
          [y, eta] = estimate_eta(model.R, k, j, bone(i, 3), X(:, i), Dim);
          X(((k-1) * Dim + 1) : k * Dim, t) = y;
          %                 elseif (W(t, k) && ~W(t, j))
          %                     [y, eta] = estimate_eta(Sigma, j, k, bone(i, 3), X(t, :)', C, ucap{t}, Dim);
          %                     X(t, ((j-1) * Dim + 1) : j * Dim) = y;
        end
      end
    end
    
    [mu, V, P, logli] = forward(X, model);
    
    %ucap is the estimated E[z_n]
    [Ez, Ezz, Ez1z] = backward(mu, V, P, model);  
  end

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
  if (exist('plotFun'))
    plotFun(X);
    drawnow;
  end
end
model = oldmodel;
Xhat = X;

function [t] = isTiny(sigma)
% test whether the matrix sigma is close to zero
%
t = (norm(sigma, 1) < 1e-10) || (any(diag(sigma) < 1e-10));



%%

function [X, A, Gamma, C, Sigma, u0, V0, LL] = learn_kalman_single_bone_fly(X, A, C, H, W, bone, maxIter, diagG, diagS, initFunc)
%learning kalman filter parameters
%X is N * M dimension
%H is the hidden dimension
%Z_n = A z_{n-1} + w_n
%X_n = C Z_n + v_n
%A: H * H, Gamma: H * H, if given, then no need to learn A; if given [],
%then to learn A
%C: M * H, Sigma: M * M, if given, then no need to learn C; if given [],
%then to learn C
%u0: H * 1, V0: H * H
%W is N * B, indicating whether x(i,j) is observed (1) or missing (0)
%         we assume the motion is in 3-D spaces, so if one marker is
%         missing there should be missing observations across three
%         dimensions. B is the number of bone joints
%         (modified by leili 2009-5-5)
%bone: .. * 3 matrix: each row with 
%         <bone 1, bone 2, bone length>
%         assume bone is symmetric. i.e. <bone 2, bone 1, bone length> will
%         also be in bone matrix.
%initFunc: 0 or not defined , initialize with linear interpolation
%          1,  initialize with SVD, also initialize C,
if (nargin < 7)
    maxIter = 100;
end

if (nargin < 8)
    diagG = 1;
    diagS = 1;
end

if (nargin < 10)
    initFunc = 0;
end

N = size(X, 1);
M = size(X, 2);
Dim = int32(M / size(W,2)); %the dimension of marker position (default 3)
% build a missing indication matrix of the same size of X
observed = reshape(repmat(W, Dim, 1), N, M);
%observed = ~missing;
missing = ~observed;
% modified by leili 2009-5-5

if (isempty(A))
    % needs to learn A
    tagA = 1;
else
    H = size(A, 1);
    tagA = 0;
end
if (isempty(C))
    % needs to learn C
    tagC = 1;
else
    H = size(C, 2);
    tagC = 0;
end
if (tagA)
%    A = eye(H, H) + randn(H, H);
    % H must be provided
    A = eye(H, H);
end
if (tagC)
    C = eye(M, H) + randn(M, H);
end

Gamma = eye(H, H);
Sigma = eye(M, M);

%% modification for the case of H > M
%by leili
%Feb 25, 2008
if (H <= M) 
    u0 = X(1, 1:H)';
else
    u0 = [X(1, :)'; zeros(H-M, 1)];
end
V0 = Gamma;


CONV_BOUND = 1e-6;



%%
%initialize mssing value
%for i=1:M
%  X(observed(:,i)==0,i) = mean(X(observed(:,i)==1,i));
%end
% if (initFunc == 1)
%     [Y, U, S, V] = EMSVD(X, observed, H, 100);
%     C = V(:,1:H) * S(1:H, 1:H);
% else
%     Y = linear_interp(X, observed);
% end
% X(missing) = Y(missing);

%Y is estimate of X
%Y = zeros(N,M);

%% initial estimation to get loglikelihood
%[u, V, P, oldLogli, Y] = forward_fly(X, A, Gamma, C, Sigma, u0, V0, observed);
%X(missing) = Y(missing);

% do not use on the fly initially, since we use linear interp
% modified by leili 2009-5-5
%[u, V, P, oldLogli] = forward(X, A, Gamma, C, Sigma, u0, V0);

%LL = oldLogli;
oldLogli = -Inf;
LL = [];

ratio = 1;
diff = 1;
iter = 0;
workingtag = 5;

%tmp = length(ek_num);
%ek = zeros(M, tmp);
%for i = 1:tmp
%    ek(tmp * (i-1) + ek_num(i), i) = 1;
%end
%tmp = length(ej_num);
%ej = zeros(M, tmp);
%for i = 1:tmp
%    ej(ej_num(i), i) = 1;
%end    


while ((workingtag~=0) && (iter < maxIter))
%while (iter < maxIter)
    iter = iter + 1;
    
    % forward belief propagation
    [u, V, P, logli, Y] = forward_fly(X, A, Gamma, C, Sigma, u0, V0, observed);
    X(missing) = Y(missing);
    
    % backward propagation
    [ucap, Vcap, J] = backward(u, V, P, A);
    for i = 1:N
        %using forward message to estimate the output layer, used in online case
        %Y(i, :) = u{i}' * C';
	
        %using backward message(posterior marginal probablility) to estimate output
        Y(i, :) = ucap{i}' * C'; 
    end
    X(missing) = Y(missing);
    
    
    %if ((rem(iter, 37) == 0) || iter>500)
    if ((iter > 100) && (rem(floor(iter / 4), 25) ~= 0))
        % do bone length constraint
        %[u, V, P, logli] = forward(X, A, Gamma, C, Sigma, u0, V0);
    
        %ucap is the estimated E[z_n]
        %[ucap, Vcap, J] = backward(u, V, P, A);
        sq = randperm(size(bone, 1));
        % estimate the missing values and lagrangian multiplier
        for t = 1 : N
            % for i = 1 : size(bone, 1)
            % use random optimization order
            for i = sq
            %for i = 1:size(bone, 1)
                k = bone(i, 1);
                j = bone(i, 2);
%                 if (~W(t, k) && W(t, j))
                if (~W(t, k))
                    %[y, eta] = estimate_eta(Sigma, k, j, bone(i, 3), X(t, :)', C, ucap{t}, Dim);
                    [y, eta] = estimate_eta(Sigma, k, j, bone(i, 3), X(t, :)', Dim);
                    X(t, ((k-1) * Dim + 1) : k * Dim) = y;
%                 elseif (W(t, k) && ~W(t, j))
%                     [y, eta] = estimate_eta(Sigma, j, k, bone(i, 3), X(t, :)', C, ucap{t}, Dim);
%                     X(t, ((j-1) * Dim + 1) : j * Dim) = y;
                end
            end
        end
        
        [u, V, P, logli] = forward(X, A, Gamma, C, Sigma, u0, V0);
    
        %ucap is the estimated E[z_n]
        [ucap, Vcap, J] = backward(u, V, P, A);
    end

%     %% debug only
%     [bb, bv, bls] = get_bones(X);
%     subplot(2, 1, 1);
%     plot([bls(:, 35, 31), bls(:, 31, 33), bls(:, 33, 21)]);
%     subplot(2, 1, 2);
%     plot(X(:, [91:93, 97:99]));
%     %subplot(3, 1, 3);
%     %plot(X(:, [91:93, 97:99]));
%     %plot(X);
%     drawnow;      
%     %pause;
    
    
    %%
    
%     if ((tagA==1) && (tagC==1))
%         [A, noG, C, Sigma, u0, V0] = MLE_kalman_bone(X, ucap, Vcap, J, observed, eta, ek, ej, diagG, diagS);
%     elseif ((tagA==1) && (tagC==0))
%         [A, Gamma, noC, noS, u0, V0] = MLE_kalman_bone(X, ucap, Vcap, J, observed, eta, ek, ej, diagG, diagS);
%     elseif ((tagA==0) && (tagC==1))
%         [noA, noG, C, Sigma, u0, V0] = MLE_kalman_bone(X, ucap, Vcap, J, observed, eta, ek, ej, diagG, diagS);
%     else 
%         [noA, noG, noC, noS, u0, V0] = MLE_kalman_bone(X, ucap, Vcap, J, observed, eta, ek, ej, diagG, diagS);
%     end
    
    %assume to learn all variables.
    [A, Gamma, C, Sigma, u0, V0] = MLE_kalman(X, ucap, Vcap, J, diagG, diagS);
    
    %%
    
    diff = abs(logli - oldLogli);
    ratio = abs(diff/logli) ;

    if ((workingtag>0) && (ratio < CONV_BOUND))
        workingtag = workingtag - 1;
    else
        workingtag = 5;
    end
    
    LL = [LL;logli];
    oldLogli = logli;
    
    
    

    
%    plot(X(:, 52:60));
%    plot(X);
%    drawnow;      
    %plot figure
    %% debug only
    [bb, bv, bls] = get_bones(X);
    subplot(2, 1, 1);
    plot([bls(:, 37, 25), bls(:, 25, 28), bls(:, 28, 38), bls(:, 28, 39), bls(:, 38, 39), bls(:, 39, 27)]);
    subplot(2, 1, 2);
    plot(X(:, [73:75, 82:84, 112:114]));
    %subplot(3, 1, 3);
    %plot(X(:, [91:93, 97:99]));
    %plot(X);
    drawnow;      
    %pause;
    
    
    
    
    fprintf('iteration = %d, logli = %d\n', iter, logli);
end


%% final smoothing
%     %for kk = 1:10
%         [u, V, P, logli] = forward(X, A, Gamma, C, Sigma, u0, V0);
%     
%         %ucap is the estimated E[z_n]
%         [ucap, Vcap, J] = backward(u, V, P, A);
%         [A, Gamma, C, Sigma, u0, V0] = MLE_kalman(X, ucap, Vcap, J, diagG, diagS);
%     %end
%     [u, V, P, logli, Y] = forward_fly(X, A, Gamma, C, Sigma, u0, V0, observed);
%     X(missing) = Y(missing);
%     
%     % backward propagation
%     [ucap, Vcap, J] = backward(u, V, P, A);
%     for i = 1:N
%         %using forward message to estimate the output layer, used in online case
%         %Y(i, :) = u{i}' * C';
% 	
%         %using backward message(posterior marginal probablility) to estimate output
%         Y(i, :) = ucap{i}' * C'; 
%     end
%     X(missing) = Y(missing);

%% 

function [x_k, eta] = estimate_eta(Sigma, k, j, bone, x, dim) 
% newton-raphson method to find the solution to
% f(eta) = 0
% where f(eta) is the difference between estimated bone and expected
% bone
% eta <- eta - f(eta) / f'(eta)
% now converted to x_bar = x(idk), i.e. x contains intial estimation in k-th bone %%%x_bar = G * z%%%%
% k is the bone# of target bone (occlusion bone),
% j is the bone# of observed one.

if nargin < 8 
    dim = 3;
end

alpha = 1; % learning rate

%x_bar = G * z;
x_bar = x;
M = size(Sigma, 1);
idk = ((k-1) * dim + 1) : (k * dim);
idj = ((j-1) * dim + 1) : (j * dim);
Sigmakj = zeros(M);
Sigmakj(:, idk) = Sigma(:, idk) - Sigma(:, idj);
Sigmakj(:, idj) = - Sigmakj(:, idk);

eta = 0;
deltaeta = 1;
MAXITER = 100;
df = 1;
iter = 0;
old_sign = 0;
f = 0;
while ((deltaeta > 1e-6) && (iter < MAXITER))
    Ieta = eye(M) + eta * Sigmakj;
    invEtax = Ieta \ x_bar;
    invEtaa = Ieta \ Sigmakj;
    deltax = invEtax(idk) - x(idj);
    df = - 2 * (invEtaa(idk, :) * invEtax)' * deltax;
    if (abs(df) > 1e-10) 
        f = sum(deltax .^ 2) - bone ^ 2;
        deltaeta = f / df;
        eta = eta - alpha * deltaeta;
        deltaeta = abs(deltaeta / eta) / alpha;
        
        if (old_sign * f < 0) 
            alpha = alpha / 2;
        end
        old_sign = sign(f);
    else 
        deltaeta = 0;
    end
    iter = iter + 1;
end
if (abs(f) > 0.1) 
    warning(sprintf('not getting the right eta when estimating %d from %d, f=%d', k, j, f));
    x_k = x_bar(idk);
else
    x_k = invEtax(idk);
end


%%
% MAXITER = 2;
% eta = NaN;
% iter = 0;
% exitflag = 0;
% while ((exitflag ~= 1) && (iter < MAXITER))
%     [eta, fval, exitflag] = fzero(@(eta)dist_diff_eta(eta, Sigma, k, j, bone, x, dim), randn);
%     iter = iter + 1;
% end
% if (exitflag ~= 1)
%     eta = 0;
%     warning(sprintf('fzero returns eta with exitflag=%d, for estimating %d from %d', exitflag, k, j));
% end
% %x_bar = G * z;
% M = size(Sigma, 1);
% idk = ((k-1) * dim + 1) : (k * dim);
% idj = ((j-1) * dim + 1) : (j * dim);
% Sigmakj = zeros(M);
% Sigmakj(:, idk) = Sigma(:, idk) - Sigma(:, idj);
% Sigmakj(:, idj) = - Sigmakj(:, idk);
% 
% Ieta = eye(M) + eta * Sigmakj;
% %invEtax = Ieta \ x_bar;
% invEtax = Ieta \ x;
% x_k = invEtax(idk);
