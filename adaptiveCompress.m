function [z, obs, x, error] = adaptiveCompress(data, ucap, A, Gamma, C, Sigma, u0, V0, ratio)
%
%adaptive compression
z=ucap;
N = size(data, 1);
M = size(data, 2);
x = zeros(N, M);
for i = 1:N
    x(i, :) = C * z{i};
end
baseerr = norm(x-data, 'fro')^2;

%hop = 4;
if (nargin < 9) 
    ratio = 4;
end

thresh = ratio * baseerr / N;
obs = zeros(N, 1);
obs(1) = 1;
mu = ucap{1};
curerr = norm(data(1,:)' - C * mu, 'fro')^2;
z{1} = ucap{1};
ss = 1;
for i = 2:N
    obs(i) = 0;
    mu = A * mu;
    curerr = curerr + norm(data(i,:)' - C * mu, 'fro')^2;
    if (curerr > thresh)
        obs(i) = 1;
        ss = ss + 1;
        z{ss} = ucap{i};
        mu = ucap{i};
        curerr = norm(data(i,:)' - C * mu, 'fro')^2;
    end
end

x = recover(z, A, Gamma, C, Sigma, u0, V0, obs);
error = norm(x-data, 'fro');
