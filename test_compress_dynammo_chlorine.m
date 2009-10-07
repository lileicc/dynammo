% test fingerprint_compress on chlorine dataset
clear;
data = load('chlorine_level_data_cl2fullLarge.dat');
X = data';
tic;
[error_dynammo_all, ratio_dynammo_all, h_dynammo_all] = runall_compress_dynammo(X, 'MaxIter', 50, 'RelativeError', 'method', 'optimal');
time_dynammo = toc;

tic;
[error_ad_all, ratio_ad_all, h_ad_all] = runall_compress_dynammo(X, 'MaxIter', 50, 'RelativeError', 'method', 'adaptive');
time_ad = toc;

[error_svd_all, ratio_svd_all, h_svd_all] = compress_pca(X);

h = figure;
hold on;
plot(ratio_svd_all, error_svd_all, 'b', 'LineWidth', 2, 'DisplayName', 'PCA');
plot(ratio_ad_all, error_ad_all, 'g', 'LineWidth', 2, 'DisplayName', 'DynaMMo_a');
plot(ratio_dynammo_all, error_dynammo_all, 'r', 'LineWidth', 2, 'DisplayName', 'DynaMMo_d');
legend show;
