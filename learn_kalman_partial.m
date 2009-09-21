function [Y, A, Gamma, C, Sigma, u0, V0, LL, RMSE, Frob] = learn_kalman_partial(X, A, C, observed, H, maxIter, diagG, diagS, initFunc)
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
%observed is N * M, indicating whether x(i,j) is observed or missing
%initFunc: 0 or not defined , initialize with linear interpolation
%          1,  initialize with SVD, also initialize C,

if (nargin < 7)
    diagG = 1;
    diagS = 1;
end

if (nargin < 9)
    initFunc = 0;
end

data = X;

N = size(X, 1);
M = size(X, 2);
%A = eye(H, H) + randn(H, H);
%C = eye(M, H) + randn(M, H);
if (isempty(A))
    tagA = 1;
else
    H = size(A, 1);
    tagA = 0;
end
if (isempty(C))
    tagC = 1;
else
    H = size(C, 2);
    tagC = 0;
end
if (tagA)
%    A = eye(H, H) + randn(H, H);
    A = eye(H, H);
end
if (tagC)
    C = eye(M, H) + randn(M, H);
end

Gamma = eye(H, H);
Sigma = eye(M, M);

%%
%initialize mssing value
%for i=1:M
%  X(observed(:,i)==0,i) = mean(X(observed(:,i)==1,i));
%end
if (initFunc == 1)
    [Y, U, S, V] = EMSVD(X, observed, H, 100);
    C = V(:,1:H) * S(1:H, 1:H);
else
    Y = linear_interp(X, observed);
end
X(observed==0) = Y(observed==0);


%modification for the case of H > M
%by leili
%Feb 25, 2008
if (H <= M) 
    u0 = X(1, 1:H)';
else
    u0 = [X(1, :)'; zeros(H-M, 1)];
end
V0 = Gamma;



%Y is estimate of X
%Y = zeros(N,M);

CONV_BOUND = 1e-6;

ratio = 1;
diff = 1;
iter = 0;

%%
[u, V, P, oldLogli] = forward(X, A, Gamma, C, Sigma, u0, V0);
LL= oldLogli;
RMSE = mean_sqerror(data, X, observed);
Frob = +inf;

workingtag = 2;

% %% for plot
% ht = plot(X);
% drawnow;


%%
while ((workingtag~=0) && (iter < maxIter))
%while (iter < maxIter)
    iter = iter + 1;
    [ucap, Vcap, J] = backward(u, V, P, A);
    
    if ((tagA==1) && (tagC==1))
        [A, Gamma, C, Sigma, u0, V0] = MLE_kalman(X, ucap, Vcap, J, diagG, diagS);
        %[A, noG, C, Sigma, u0, V0] = MLE_kalman(X, ucap, Vcap, J, diagG, diagS);
    elseif ((tagA==1) && (tagC==0))
        [A, Gamma, noC, noS, u0, V0] = MLE_kalman(X, ucap, Vcap, J, diagG, diagS);
    elseif ((tagA==0) && (tagC==1))
        [noA, noG, C, Sigma, u0, V0] = MLE_kalman(X, ucap, Vcap, J, diagG, diagS);
    else 
        [noA, noG, noC, noS, u0, V0] = MLE_kalman(X, ucap, Vcap, J, diagG, diagS);
    end
    
    [u, V, P, logli] = forward(X, A, Gamma, C, Sigma, u0, V0);
    
    for i = 1:N
	%using forward message to estimate the output layer, used in online case
        %Y(i, :) = u{i}' * C';
	
	%using backward message(posterior marginal probablility) to estimate output
	Y(i, :) = ucap{i}' * C'; 
    end
    X(observed==0) = Y(observed==0);
    
    %logli = real(logli);
    diff = abs(logli - oldLogli);
    ratio = abs(diff/logli) ;
    
   
    LL = [LL;logli];
    oldLogli = logli;
    
    RMSE = [RMSE; mean_sqerror(data, X, observed)];
    thisfrob = norm(Y - X, 'fro');
    Frob = [Frob; thisfrob];
    
    %if ((workingtag>0) && (ratio < CONV_BOUND) && (diff < CONV_BOUND))
    if ((workingtag>0) && (ratio < CONV_BOUND))
        workingtag = workingtag - 1;
    else
        workingtag = 2;
    end
    
     fprintf('iteration = %d, logli = %d, frob=%d\n', iter, logli, thisfrob);
     
     %plot figure
     %
%      plot(X);
%      drawnow;
end
