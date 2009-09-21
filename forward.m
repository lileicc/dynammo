function [mu, V, P, logli] = forward(X, model)
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
% Returns:
%   mu is cell 1 * N, each with a matrix H * 1
%   V is cell 1 * N, each with a matrix H * H
%   P is cell 1 * N, each with a matrix H * H
%   logli: log-likelihood of the data under the current model.
%   
%obs is optional argument for indicating whether x(obs(i)) is observed or
%missing (modified Sep 11)
%version 29, by leili

%A, Gamma, C, Sigma, u0, V0,  obs

N = size(X, 2);
M = size(X, 1);
H = size(model.A, 1); %dimension of hidden variable
Ih = eye(H, H);

% if (nargin < 8) 
%     obs = ones(N,1);
% end

%predicted mean for hidden variable z
mu = cell(1, N);
V = cell(1, N);
P = cell(1, N);

% %initialize
% if (obs(1) == 1)
%     %if x_1 is observed, this is normal forward in LDS
%     u_c = C * u0;
%     sigma_c = C * V0 * C' + Sigma;
% %    R = chol(sigma_c);
% %    detSigmaC = det(R) ^ 2;
% %    R = inv(R);
% %    invSig = R * R';
%     %invSig = inv(sigma_c);
%     %could we use /
%     %K = V0 * C' * invSig;
%     K = V0 * C' / sigma_c;
%     u{1} = u0 + K * (x(1,:)' - u_c);
%     V{1} = (Ih - K * C) * V0;
%     delta = x(1, :) - u_c';
%     logdetSigmaC = logDet(sigma_c);
%     %if (detSigmaC < 0) 
%     %    warning('det of sigma_c < 0');
%     %end
%     posDef = (delta / sigma_c) * delta'/2;
%     if (posDef < 0)
%         warning('det of not positive definite < 0');
%     end 
%     logli = - M/2 * log(2 * pi) - logdetSigmaC / 2 - posDef;
% else 
%     %if x_1 is missing, this is just initial value
%     u{1} = u0;
%     V{1} = V0;
%     logli = 0;
% end

% initialize
mu{1} = model.mu0;
V{1} = model.Q0;
logli = 0;

for i = 2:N
    P{i-1} = model.A * V{i-1} * model.A' + model.Q;
    
%     if (obs(i) == 1)
        %if x_i is observed, this is normal forward in LDS    
        sigma_c = model.C * P{i-1} * model.C' + model.R;  
        invSig = pinv(sigma_c);
        %use /
        %K = P{i-1} * C' * invSig;        
        K = P{i-1} * model.C' * invSig;
        mu{i} =  model.A * mu{i-1};
        u_c = model.C * mu{i};
        delta = X(:, i) - u_c;
        mu{i} = mu{i} + K * delta;
        V{i} = (Ih - K * model.C) * P{i-1};
        posDef = delta' * invSig * delta / 2;
        if (posDef < 0)
            warning('det of not positive definite < 0');
        end 
        logli = logli - M/2 * log(2 * pi) - logdet(sigma_c, 'chol') / 2 - posDef;
%     else
%         %if x_i is missing, this is just initial value
%         u{i} = A * u{i-1};
%         V{i} = P{i-1};
%     end
end