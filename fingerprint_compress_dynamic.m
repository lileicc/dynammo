function [errors, ratios, D, mu0, Pr, Qc, Qs, Qm, Rc, Rs, Rm, P, Y, zhat] = fingerprint_compress_dynamic(X, varargin)
% kalman fingerprinting (PLF method) compressiong,
% compression using dynamic programming on the optimal error 
%
%   X: M * N matrix, M is number of sequences, N is the time duration.
%
% the usage is like: fingerprint_compress_dynamic(X, 'Hidden', 10, 'MaxIter',
% 100)
%
% Option:
% 'MaxRatio': followed by a number max allowed ratio for compression

N = size(X, 2);
M = size(X, 1);


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
   
    %Rnum = find(Rsum > (Rsum(end) * 0.95), 1);
    Rnum = size(R, 2);
    Rc = Rcoeff(:, 1:Rnum);
    Rs = Rscore(:, 1:Rnum);
else
    PP = P;
end

Pr = PP(:, 1:num_real);
Q = abs(PP(:, (num_real+1):end));
Qm = mean(Q);
[Qcoeff, Qscore, Qlatent] = princomp(Q, 'econ');
Qsum = cumsum(Qlatent);
%Qnum = find(Qsum > (Qsum(end) * 0.95), 1);
Qnum = size(Q, 2);
Qc = Qcoeff(:, 1:Qnum);
Qs = Qscore(:, 1:Qnum);

%% dynamic compression
model = struct;
model.A = diag(D);
model.C = oldP;
model.mu0 = mu0;
H = length(model.mu0);
ORIGINAL = numel(X) + 2;
Ex = zhat;
%[errs, last] = dynamicCompress(X, zhat, diag(D), [], oldP, [], mu0, [], maxRatio);
%[errors, ratios] = compress_dynammo(X, 'model', model, zhat, varargin{:});
baseerr = zeros(N + 1, 1);
for i = 1 : N
  baseerr(i) = norm(X(:, i) - real(model.C * Ex{i}), 'fro') ^ 2;
end


a = find(strcmp('MaxSpaceRatio', varargin), 1);
if (isempty(a))
  maxRatio = 1;
else
  maxRatio = varargin{a+1};
end
thresh = ceil(N*maxRatio);
errs = zeros(N + 1, thresh);
last = zeros(N + 1, thresh);
errs(1, 1) = baseerr(1);

for ss = 1 : N
  mu = Ex{ss};
  tt = ss + 1;
  curerr = 0;
  while (tt <= N + 1)
    mu = model.A * mu;
    if (ss > 1)
      i = 2;
    else
      i = 1;
    end
    while ((i <= ss) && (i < thresh))
      tmperr = curerr + baseerr(tt) + errs(ss, i);
      if ((last(tt, i+1) == 0) || (tmperr < errs(tt, i + 1)))
        errs(tt, i+1) = tmperr;
        last(tt, i+1) = ss;
      end
      i = i + 1;
    end
    if (tt <= N)
      curerr = curerr + norm(X(:, tt) - real(model.C * mu), 'fro')^2;
    end
    tt = tt+1;
  end
  if (ss <= 1)
    errs(N+1, 1) = curerr;
  end
end

errors = errs(N+1, :);
OVERHEAD = numel(model.mu0) * 2 + numel(model.C) + 4;
K = length(errors);
store_space = ((1:K) - 1) .* (H + 1) + OVERHEAD;
ratios = ORIGINAL ./ store_space;
comp_data = cell(K, 1);
TAG = 3;
for i = 1 : K
  j = i;
  tt = N + 1;
  temp = [];
  while (last(tt, j) > 1)
    tt = last(tt, j);
    temp = [temp, tt];
  end
  comp_data{i} = zeros(store_space(i), 1);
  comp_data{1}(1:OVERHEAD) = [TAG; N; M; H; model.mu0; reshape(diag(model.A), [], 1); reshape(model.C, [], 1)];
  ss = OVERHEAD + 1;
  for tt = temp(end:(-1):1)
    comp_data{i}(ss) = tt;
    comp_data{i}((ss+1) : (ss+H)) = Ex{tt};
    ss = ss + H + 1;
  end
end

if (any(strcmp('RelativeError', varargin)))
  X_m = mean(X, 2);
  TOTALVAR = norm(X - repmat(X_m, 1, size(X, 2)), 'fro') ^ 2;
  errors = errors ./ TOTALVAR;
end

% %% reconstruction
% %Q = Qs * Qc' + repmat(Qm, size(Qs, 1), 1);
% %R = Rs * Rc' + repmat(Rm, size(Rs, 1), 1);
% k = size(R, 2);
% tmp = reshape([Q(:, (end-k+1):end) .* exp(1i .* R); Q(:, (end-k+1):end) .* exp(-1i .* R)], size(Q, 1), []);
% P = [Pr, tmp];
% %tmpd = reshape([D((end-k+1):end), conj(D((end-k+1):end))]', [], 1);
% %D = [D(1:(end-k)); tmpd];
% 
% z = mu0;
% Y = X;
% for k = 1 : N
%     Y(:, k) = real(P * z);
%     z = D .* z;
% end
% 
% error = norm(Y - X, 'fro') / norm(X - repmat(mean(X), N, 1), 'fro');
% %ratio = (N * M + 2) / (2 + (2 * HIDDEN + 1) + numel(Pr) + numel(Qs) + numel(Qc) + numel(Qm) + 1 + numel(Rs) + numel(Rc) + numel(Rm) + 1);
% ratio = (N * M + 2) / (2 + (2 * HIDDEN + 1) + numel(Pr) + numel(Q) + numel(R));
% 
% P = oldP;