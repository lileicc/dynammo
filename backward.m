function [ucap, Vcap, J] = backward(u, V, P, A)
%backward calculation of gamma(z) = alpha(z) * beta(z)
%ucap: cell N * 1 with each H * 1
%Vcap: cell N * 1 with each H * H
%P: cell (N-1) * 1 with each H * H
%backward algorithm is the same no matter x_i is observed or missing.
%version 29, by leili

N = size(u, 1);
%H = size(u, 2);
%ucap = zeros(N, H);
%Vcap = zeros(N, H, H);
%J = zeros(N, H, H);
ucap = cell(N,1);
Vcap = cell(N,1);
J = cell(N, 1);

ucap{N} = u{N};
Vcap{N} = V{N};

%delta = A * u(N,:)';

for i = (N-1): (-1) : 1
    %J{i} = V{i} * A' * inv(P{i});
    J{i} = V{i} * A' / P{i};
    ucap{i} = u{i} + J{i} * (ucap{i+1} - A * u{i});
    Vcap{i} = V{i} + J{i} * (Vcap{i+1} - P{i}) * J{i}';
end
