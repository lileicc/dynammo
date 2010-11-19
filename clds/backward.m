function [Ez, Ezz, Ez1z] = backward(u, UU, P, model)
% backward calculation of posteriors P(z_n | x_1...x_N)
% 
% Ez: cell 1 * N with each H * 1, E[ z_n | x_1...x_N ]
% Ezz: cell 1 * N with each H * H, E[ z_n * z_n' | x_1...x_N ]
% Ezz1: cell (N-1) * 1 with each H * H, E[ z_n * z_n+1 | x_1 ...x_N ]
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
%

N = length(u);
Ez = cell(1, N);
Ezz = cell(1, N);
Ez1z = cell(1, N);
Ez{N} = u{N};
Vhat = UU{N};
Ezz{N} = Vhat + Ez{N} * Ez{N}';

for i = (N-1): (-1) : 1
    J = UU{i} * model.A' / P{i+1}; % this is the back-propagation factor
    Ez{i} = u{i} + J * (Ez{i+1} - model.A * u{i});
    Ez1z{i} = Vhat * J' + Ez{i+1} * Ez{i}';
    Vhat = UU{i} + J * (Vhat - P{i+1}) * J';
    Ezz{i} = Vhat + Ez{i} * Ez{i}';
end


