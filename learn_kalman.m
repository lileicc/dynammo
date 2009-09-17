function [A, Gamma, C, Sigma, u0, V0, LL] = learn_kalman(x, H, maxIter)
% learning kalman filter parameters
% x is N * M dimension
% z_n = A z_{n-1} + w_n
% x_n = C x_n + v_n
% A: H * H, Gamma: H * H
% C: M * H, Sigma: M * M
% u0: H * 1, V0: H * H
% version 29, by leili

N = size(x, 1);
M = size(x, 2);
A = eye(H, H) + randn(H, H);
C = eye(M, H) + randn(M, H);
%A = eye(H, H);
%C = eye(M, H);
Gamma = eye(H, H);
Sigma = eye(M, M);
m = mean(x);

%modification for the case of H > M
%by leili
%Feb 25, 2008
if (H <= M) 
    %u0 = x(1, 1:H)';
    u0 = randn(H, 1);
else
    u0 = [x(1, :)'; ones(H-M, 1)];
end
V0 = Gamma;

CONV_BOUND = 1e-5;

ratio = 1;
diff = 1;
iter = 0;
[u, V, P, oldLogli] = forward(x, A, Gamma, C, Sigma, u0, V0);

while (ratio > CONV_BOUND || diff > CONV_BOUND) && (iter < maxIter)
    iter = iter + 1;
    [ucap, Vcap, J] = backward(u, V, P, A);
    [A, Gamma, C, Sigma, u0, V0] = MLE_kalman(x, ucap, Vcap, J);
    [u, V, P, logli] = forward(x, A, Gamma, C, Sigma, u0, V0);
    logli = real(logli);
    diff = (logli - oldLogli);
    ratio = abs(diff/logli) ;
    LL(iter) = logli;
    oldLogli = logli;
    fprintf('iteration = %d, logli = %d\n', iter, logli);
end