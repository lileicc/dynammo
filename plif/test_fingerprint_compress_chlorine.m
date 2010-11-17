% test fingerprint_compress on chlorine dataset
clear;
data = load('chlorine_level_data_cl2fullLarge.dat');
X = data';
tic;
[error_finger_all, ratio_finger_all, h_finger_all] = runall_fingerprint_compress_dynamic(X, 'Hidden', [60:(-5):5, 1:4], 'MaxIter', 50, 'RelativeError', 'Fast');
time_finger = toc;

tic;
[error_dynammo_all, ratio_dynammo_all, h_dynammo_all] = runall_compress_dynammo(X, 'Hidden', [60:(-5):5, 1:4], 'MaxIter', 50, 'RelativeError', 'Fast');
time_dynammo = toc;

[error_svd_all, ratio_svd_all, h_svd_all] = compress_pca(X);

h = figure;
hold on;
plot(ratio_svd_all, error_svd_all, 'b--', 'LineWidth', 2, 'DisplayName', 'PCA');
plot(ratio_dynammo_all, error_dynammo_all, 'g:', 'LineWidth', 2, 'DisplayName', 'DynaMMo');
plot(ratio_finger_all, error_finger_all, 'r', 'LineWidth', 2, 'DisplayName', 'PLiF');
legend show;
