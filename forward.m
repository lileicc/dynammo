function [u, V, P, logli] = forward(x, A, Gamma, C, Sigma, u0, V0, obs)
%the forward message passing for LDS
%using forward algorithm to inference x_(N+1) given X
%x is N * M dimension
%z_n = A z_{n-1} + w_n
%x_n = C x_n + v_n
%A: H * H, Gamma: H * H
%C: M * H, Sigma: M * M
%u0: H * 1, V0: H * H
%u is cell N * 1
%V is cell N * H * H
%P is cell N * H * H
%obs is optional argument for indicating whether x(obs(i)) is observed or
%missing (modified Sep 11)
%version 29, by leili


N = size(x, 1);
M = size(x, 2);
H = size(A, 1); %dimension of hidden variable
Ih = eye(H, H);

if (nargin < 8) 
    obs = ones(N,1);
end

%c_n = p(x_n|x_1, ... x_{n-1}) ~ N(x_n| u_c(n), sigma_c(n))
u_c = zeros(M, 1); 
sigma_c = zeros(M, M);

%predicted mean for hidden variable z
%u = zeros(N, H);
%V = zeros(N, H, H);
u = cell(N, 1);
V = cell(N, 1);

%P = zeros(N, H, H);
P = cell(N, 1);

%initialize
if (obs(1) == 1)
    %if x_1 is observed, this is normal forward in LDS
    u_c = C * u0;
    sigma_c = C * V0 * C' + Sigma;
%    R = chol(sigma_c);
%    detSigmaC = det(R) ^ 2;
%    R = inv(R);
%    invSig = R * R';
    %invSig = inv(sigma_c);
    %could we use /
    %K = V0 * C' * invSig;
    K = V0 * C' / sigma_c;
    u{1} = u0 + K * (x(1,:)' - u_c);
    V{1} = (Ih - K * C) * V0;
    delta = x(1, :) - u_c';
    logdetSigmaC = logDet(sigma_c);
    %if (detSigmaC < 0) 
    %    warning('det of sigma_c < 0');
    %end
    posDef = (delta / sigma_c) * delta'/2;
    if (posDef < 0)
        warning('det of not positive definite < 0');
    end 
    logli = - M/2 * log(2 * pi) - logdetSigmaC / 2 - posDef;
else 
    %if x_1 is missing, this is just initial value
    u{1} = u0;
    V{1} = V0;
    logli = 0;
end


for i = 2:N
    P{i-1} = A * V{i-1} * A' + Gamma;
    
    if (obs(i) == 1)
        %if x_i is observed, this is normal forward in LDS    
        sigma_c = C * P{i-1} * C' + Sigma;  
        %invSig = inv(sigma_c);
        %use /
        %K = P{i-1} * C' * invSig;
        K = P{i-1} * C' / sigma_c;
        u{i} =  A * u{i-1};
        u_c = C * u{i};
        u{i} = u{i} + K * (x(i, :)' - u_c);
        V{i} = (Ih - K * C) * P{i-1};
        delta = x(i, :) - u_c';
    
        logdetSigmaC = logDet(sigma_c);
        %if (detSigmaC < 0) 
        %    warning('det of sigma_c < 0');
        %end
        posDef = (delta / sigma_c) * delta'/2;
        if (posDef < 0)
            warning('det of not positive definite < 0');
        end 
        logli = logli - M/2 * log(2 * pi) - logdetSigmaC / 2 - posDef;
    else
        %if x_i is missing, this is just initial value
        u{i} = A * u{i-1};
        V{i} = P{i-1};
    end
end