% test fingerprint_compress
clear;
load 'data_global_new.mat';

%data = load('../data/chlorine_level_data_cl2fullLarge.dat');

X = data{22}(:, 4:96)';
%X = data;
error_finger_all = [];
ratio_finger_all = [];
error_svd_all = [];
ratio_svd_all = [];

Gap = 5;
for iter = 1:(size(X, 1) / Gap)
    HIDDEN = iter * Gap;
%     HIDDEN = 10;

    [error_f, ratio_f, D, mu0, Pr, Qc, Qs, Qm, Rc, Rs, Rm, P] = fingerprint_compress_dynamic(X, 'Hidden', HIDDEN, 'MaxIter', 10, 'RelativeError');
    error_finger_all = [error_finger_all, error_f];
    ratio_finger_all = [ratio_finger_all, ratio_f];
    
    %fprintf('Hidden = %d, error = %d, ratio = %d', HIDDEN, error_f, ratio_f);   
    fprintf('Hidden = %d, error = %d, ratio = %d', HIDDEN, error_f(1), ratio_f(1));   
end    
[tmpvalue, tmpIdx] = sort(ratio_finger_all);
ratio_finger_all = tmpvalue;
error_finger_all = error_finger_all(tmpIdx);

X_m = mean(X, 2);
[Coeff, Score] = princomp(X', 'econ');
TOTALVAR = norm(X - repmat(X_m, 1, size(X, 2)), 'fro'); 
if (size(X, 1) < 20)
  cands = 1:size(X, 1);
else
  cands = [1:20, 25: 5 :size(X, 1)];
end
for HIDDEN = cands
    error_svd = norm(Coeff(:, 1:HIDDEN) * Score(:, 1:HIDDEN)' + repmat(X_m, 1, size(X, 2)) - X, 'fro') / TOTALVAR;
    ratio_svd = (numel(X) + 2) / (size(X, 2) * (HIDDEN) + size(X, 1) * (HIDDEN + 1) + 3);
    error_svd_all = [error_svd_all, error_svd];
    ratio_svd_all = [ratio_svd_all, ratio_svd];
end     


h = figure;

hold on;

plot(ratio_svd_all, error_svd_all, 'b', 'LineWidth', 2, 'DisplayName', 'SVD');

plot(ratio_finger_all, error_finger_all, 'r', 'LineWidth', 2, 'DisplayName', 'PLF');

legend show;


% [P, D, mu0] = fingerprint(X, 'Hidden', 10, 'Iteration', 100);
% 
% z = mu0;
% Y = X;
% for k = 1 : size(X, 1)
%     Y(k, :) = real(P * z);
%     z = D .* z;
% end
% 
% error = norm(Y - X, 'fro') / norm(X - repmat(mean(X), size(X, 1), 1), 'fro');


HIDDEN = 50;
%MAX_ITER = 100;
HIDDEN = 10;
    [ucap, A, Gamma, C, Sigma, u0, V0] = learn_kalman(X, HIDDEN, 100);
    [z, obs, x, error_k] = fixedCompress(X, ucap, A, Gamma, C, Sigma, u0, V0, size(X, 1));
    
    
    error_k = error_k / norm(X - repmat(mean(X), size(X, 1), 1), 'fro');
    
    z = u0;
    Y = X;
    for k = 1 : size(X, 1)
        Y(k, :) = C * z;
        z = A * z;
    end
    error_kalman = norm(X - Y, 'fro') / (norm(X - repmat(mean(X), size(X, 1), 1), 'fro'));
    
    figure;
    hold all;
    plot(X(:, 1));
    plot(Y(:, 1));