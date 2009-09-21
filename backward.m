function [Ez, Ezz, Ez1z] = backward(mu, V, P, model)
% backward calculation of gamma(z) = alpha(z) * beta(z)
% 
% Ez: cell 1 * N with each H * 1, E[ z_n | x_1...x_N ]
% Ezz: cell 1 * N with each H * H, E[ z_n * z_n' | x_1...x_N ]
% Ez1z: cell (N-1) * 1 with each H * H, E[ z_n+1 * z_n | x_1 ...x_N ]
% backward algorithm is the same no matter x_i is observed or missing.

N = length(mu);
Ez = cell(1, N);
Ezz = cell(1, N);
Ez1z = cell(1, N);
Ez{N} = mu{N};
Vhat = V{N};
Ezz{N} = Vhat + Ez{N} * Ez{N}';

for i = (N-1): (-1) : 1
    J = V{i} * model.A' / P{i};
    Ez{i} = mu{i} + J * (Ez{i+1} - model.A * mu{i});
    Ez1z{i} = Vhat * J' + Ez{i+1} * Ez{i}';
    Vhat = V{i} + J * (Vhat - P{i}) * J';
    Ezz{i} = Vhat + Ez{i} * Ez{i}';
end

