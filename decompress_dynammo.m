function [X] = decompress_dynammo(data)
% decompress from compressed data using dynammo learned features,
% there are three compression methods: fixed hop, adaptive, and optimal by
% dynamic programming.
% Args:
%   data: a vectore of numbers
%     the first number identifies the compression type: 
%     1 --> fixed
%     2 --> adaptive
%     3 --> optimal
%     the second, third, firth numbers are N, M, H
%     followed by mu0, A, C
%     the remaining will be 
%     for fixed, a number for hop, else, each H numbers for H
%     hidden variables for the time tick.
%     for adaptive and optimal, each H+1 numbers will be time tick, and H
%     hidden variables for the time tick.
%
% Returns:
%   X: decompressed data
%
% $Author$
% $Date$
% $Revision$

N = size(obs, 1);
M = size(C, 1);
H = size(A, 1); %dimension of hidden variable
Ih = eye(H, H);

if (nargin < 8) 
    obs = ones(N,1);
end


%predicted mean for hidden variable z
%u = zeros(N, H);
%V = zeros(N, H, H);
u = cell(N, 1);
V = cell(N, 1);
ucap = cell(N, 1);


%P = zeros(N, H, H);
P = cell(N, 1);
J = cell(N, 1);

last = 0;
head = 0;
%initialize
if (obs(1) == 1)
    %if x_1 is observed, this is normal forward in LDS   
    
    %option 1: use 0
    u{1} = z{1};
    V{1} = zeros(H, H);
    %option 2: use 
    %sigma_c = C * V0 * C' + Sigma;
    %K = V0 * C' / sigma_c;
    %u{1} = u0;
    %V{1} = (Ih - K * C) * V0;
    
    ucap{1} = z{1};
    last = 1;
    head = 1;
else 
    %if x_1 is missing, this is just initial value
    u{1} = u0;
    V{1} = V0;
end

for i = 2:N
    P{i-1} = A * V{i-1} * A' + Gamma;
    
    %option 1: use P
    %V{i} = P{i-1}; 
    %option 2: use other
    sigma_c = C * P{i-1} * C' + Sigma;  
    K = P{i-1} * C' / sigma_c;
    V{i} = (Ih - K * C) * P{i-1}; 
    
    if (obs(i) == 0)    
        u{i} =  A * u{i-1};
    else
        head = head + 1;
        ucap{i} = z{head};
        u{i} = z{head};
        V{i} = zeros(H, H);
        
        for k = (i-1): (-1) : (last+1)
            J{k} = V{k} * A' / P{k};
            ucap{k} = u{k} + J{k} * (ucap{k+1} - A * u{k});
        end
        last = i;
    end
end

if (obs(N) == 0)
    for k =(last+1):1:N
        ucap{k} = u{k};
    end
end

x = zeros(N, M);
for i = 1:N
    x(i, :) = C * ucap{i};
end