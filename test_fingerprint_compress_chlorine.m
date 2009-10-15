% test fingerprint_compress on chlorine dataset
clear;
data = load('chlorine_level_data_cl2fullLarge.dat');
X = data';
tic;
[error_finger_all, ratio_finger_all, h_finger_all] = runall_fingerprint_compress_dynamic(X, 'Hidden', [60:(-5):5, 1:4] , 'MaxIter', 50, 'RelativeError');
time_finger = toc;

[error_svd_all, ratio_svd_all, h_svd_all] = compress_pca(X);

h = figure;
hold on;
plot(ratio_svd_all, error_svd_all, 'b', 'LineWidth', 2, 'DisplayName', 'PCA');
plot(ratio_finger_all, error_finger_all, 'r', 'LineWidth', 2, 'DisplayName', 'PLF');
legend show;
