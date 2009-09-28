function [P, D, mu0, zhat] = fingerprint(X, varargin)
% kalman fingerprinting (PLF method)
%
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
%
% the usage is like: fingerprint(X, 'Hidden', 10, 'Iteration', 100)
%
% Example: 
%
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$

model = learn_lds(X, varargin{:});

%
[VV, DD] = eig(model.A); % take the eigen value decomposition of A
[no, ind] = sort(abs(imag(diag(DD)))); % sort according the (abs) imagenary part
D = diag(DD(ind, ind)); % rearrange
V = VV(:, ind); % rearrange
P = model.C * V;
mu0 = V \ model.mu0; % mapped initial

% get the hidden states
[u_k, V_k, P_k] = forward(X, model);
[Ex] = backward(u_k, V_k, P_k, model);
zhat = cell(size(Ex));
for i = 1 : length(zhat)
    zhat{i} = V \ Ex{i};
end