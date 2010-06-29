%FINGERPRINT extract features from time series data using PLiF method

function [P, D, mu0, zhat] = fingerprint(X, varargin)
% feature extraction for time series data (PLiF method)
%
%  Lei Li, B. Aditya Prakash, Christos Faloutsos. 
%  Parsimonious Linear Fingerprinting for Time Series. 
%  VLDB 2010. Singapore
%
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
% 
% Out:
%   D: complex diagonal matrix describing the dynamics of hidden varibles.
%   P: the (complex) matrix mapping from hidden varible to observation
%   mu0: initial hidden varible
%   zhat: hidden variables under PLiF model (might be complex)   
%
% the usage is like: fingerprint(X, 'Hidden', 10, 'MaxIter', 100)
%
% Example: 
%
%   See also learn_lds
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$

model = learn_lds(X, varargin{:}); % learn parameters of LDS

%
[VV, DD] = eig(model.A); % take the eigen value decomposition of A
[no, ind] = sort(abs(imag(diag(DD)))); % sort according the (abs) imagenary part
D = diag(DD(ind, ind)); % rearrange
V = VV(:, ind); % rearrange
P = model.C * V;
mu0 = V \ model.mu0; % mapped initial

% get the hidden states
[u_k, V_k, P_k] = forward(X, model, varargin{:});
[Ex] = backward(u_k, V_k, P_k, model);
zhat = cell(size(Ex));
for i = 1 : length(zhat)
    zhat{i} = V \ Ex{i};
end