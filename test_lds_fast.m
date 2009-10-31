
clear;
load 'data_global_new.mat';
X = data{22}(:, 4:96)';
X = data{1}(:, 4:96)';
X = data{45}(:, 4:96)';

tic;
[model_fast] = learn_lds(X, 'Hidden', 15, 'MaxIter', 100, 'Fast');
time_fast = toc;

tic;
[model] = learn_lds(X, 'Hidden', 15, 'MaxIter', 100);
time = toc;