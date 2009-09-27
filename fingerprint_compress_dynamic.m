function [errors, ratios, D, mu0, Pr, Qc, Qs, Qm, Rc, Rs, Rm, P, Y, zhat] = fingerprint_compress_dynamic(X, varargin)
% kalman fingerprinting (PLF method) compressiong,
% compression using dynamic programming on the optimal error 
%
% X: N * M matrix, N is the duration, M is number of sequences
%
% the usage is like: fingerprint_compress_dynamic(X, 'Hidden', 10, 'Iteration',
% 100)
%
% Option:
% 'MaxRatio': followed by a number max allowed ratio for compression

N = size(X, 1);
M = size(X, 2);


[P, D, mu0, zhat] = fingerprint(X, varargin{:});
oldP = P;

% number of hidden dimension
% usually put an even number
HIDDEN = length(D);

ind = find(abs(imag(D)) > 1e-10); % assume 0==1e-10
num_real = HIDDEN;
% further eliminate the conjugate ones
if (~ isempty(ind)) 
    num_real = ind(1) - 1; % the eigenvalues before ind(1) are assumed to be real values.
    subind = [(1:(ind(1)-1))'; ind(1:2:end)];
    PP = P(:, subind);
    %D = D(subind);

    R = angle(PP(:, (num_real+1):end)); %only the imaginary part
    Rm = mean(R);
    [Rcoeff, Rscore, Rlatent] = princomp(R, 'econ');
    Rsum = cumsum(Rlatent);
   
    Rnum = find(Rsum > (Rsum(end) * 0.95), 1);
    %Rnum = size(R, 2);
    Rc = Rcoeff(:, 1:Rnum);
    Rs = Rscore(:, 1:Rnum);
end

Pr = PP(:, 1:num_real);
Q = abs(PP(:, (num_real+1):end));
Qm = mean(Q);
[Qcoeff, Qscore, Qlatent] = princomp(Q, 'econ');
Qsum = cumsum(Qlatent);
Qnum = find(Qsum > (Qsum(end) * 0.95), 1);
%Qnum = size(Q, 2);
Qc = Qcoeff(:, 1:Qnum);
Qs = Qscore(:, 1:Qnum);

%% dynamic compression

a = find(strcmp('MaxRatio', varargin));
if (isempty(a))
    maxRatio = 0.99;
else
    maxRatio = varargin{a+1};
end
[errs, last] = dynamicCompress(X, zhat, diag(D), [], oldP, [], mu0, [], maxRatio);
errors = sqrt(errs(N, :)) / norm(X - repmat(mean(X), N, 1), 'fro');
ratios = (N * M + 2) ./ ((1:length(errors)) .* HIDDEN + 2 + (HIDDEN + 1) + numel(Pr) + numel(Q) + numel(R));



%% reconstruction
%Q = Qs * Qc' + repmat(Qm, size(Qs, 1), 1);
%R = Rs * Rc' + repmat(Rm, size(Rs, 1), 1);
k = size(R, 2);
tmp = reshape([Q(:, (end-k+1):end) .* exp(1i .* R); Q(:, (end-k+1):end) .* exp(-1i .* R)], size(Q, 1), []);
P = [Pr, tmp];
%tmpd = reshape([D((end-k+1):end), conj(D((end-k+1):end))]', [], 1);
%D = [D(1:(end-k)); tmpd];

z = mu0;
Y = X;
for k = 1 : N
    Y(k, :) = real(P * z);
    z = D .* z;
end

error = norm(Y - X, 'fro') / norm(X - repmat(mean(X), N, 1), 'fro');
%ratio = (N * M + 2) / (2 + (2 * HIDDEN + 1) + numel(Pr) + numel(Qs) + numel(Qc) + numel(Qm) + 1 + numel(Rs) + numel(Rc) + numel(Rm) + 1);
ratio = (N * M + 2) / (2 + (2 * HIDDEN + 1) + numel(Pr) + numel(Q) + numel(R));

P = oldP;