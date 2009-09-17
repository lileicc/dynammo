function [A, Gamma, C, Sigma, u0, V0] = MLE_kalman(x, ucap, Vcap, J, diagG, diagS)
%Maximum likelihood estimation of kalman filter learning in one step
%x: N * M
%ucap: N * H
%A: H * H, Q: H * H
%C: M * H, R: M * M
%u0: H * 1, V0: H * H
%diagG, diagS: whether to put diagonal (1), diagonal and same (2), matrix (3)
%version 29, by leili

if (nargin < 5) 
    diagG = 1;
    diagS = 1;
end

N = size(x, 1);
M = size(x, 2);
H = size(ucap{1}, 1);
%Ezz1 = zeros(N-1, H, H);
%Ezz = zeros(N, H, H);
Ezz1 = cell(N-1, 1);
Ezz = cell(N, 1);

Szz1 = zeros(H, H);
Szz = zeros(H, H);
Sxz = zeros(M, H);

for i = 1:(N-1)
    % old version, there is on bug here
    % Ezz1{i} = J{i} * Vcap{i+1} + ucap{i+1} * ucap{i}';
    
    % new version
    Ezz1{i} = Vcap{i+1} * J{i}' + ucap{i+1} * ucap{i}';
    Szz1 = Szz1 + Ezz1{i};
end

for i = 1:N
    Ezz{i} = Vcap{i} + ucap{i} * ucap{i}';
    Szz = Szz + Ezz{i};
    Sxz = Sxz + x(i, :)' * ucap{i}';
end

SzzN = Szz - Ezz{N}; % sum of E[z, z] from 1 to n-1


u0 = ucap{1};
V0 = Vcap{1};
%A = Szz1 * inv(SzzN);
A = Szz1 / SzzN;
if (diagG==1) 
    delta = (trace(Szz) - trace(Ezz{1}) - 2 * trace(A * Szz1') + trace(A * SzzN * A')) / (N-1) / H;
    Gamma = diag(repmat(delta, H, 1));
elseif (diagG==2)
    %use diagonal covariance
    Gamma = diag((diag(Szz) - diag(Ezz{1}) - 2 * diag(A * Szz1') + diag(A * SzzN * A')) / (N-1));
else 
    Gamma = (Szz - Ezz{1} - A * Szz1' - Szz1 * A' + A * SzzN * A') / (N-1);
end

%Sxz = x' * ucap;
%C = Sxz * inv(Szz);
C = Sxz / Szz;


if (diagS == 1) 
    delta = (trace(x' * x) - 2 * trace(C * Sxz') + trace(C * Szz * C')) / N / M;
    Sigma = diag(repmat(delta, M, 1));
elseif (diagS == 2)
    Sigma = diag((diag(x' * x) - 2 * diag(C * Sxz') + diag(C * Szz * C')) / N);
else
    Sigma = (x' * x - C * Sxz' - Sxz * C' + C * Szz * C') / N;
end