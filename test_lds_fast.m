
clear;
load 'data_global_new.mat';

%X = data{22}(:, 4:96)';
%X = data{1}(:, 4:96)';

N = length(data);
time_fast = zeros(N, 1);
time_basic = zeros(N, 1);
for i = 1 : N
  X = data{i}(:, 4:96)';
  
  tic;
  [model_fast] = learn_lds(X, 'Hidden', 15, 'MaxIter', 100, 'Fast');
  time_fast(i) = toc;
  
  tic;
  [model] = learn_lds(X, 'Hidden', 15, 'MaxIter', 100);
  time_basic(i) = toc;
end

figure;
scatter(time_basic, time_fast);
hold on;
xx = xlim;
plot(xx, xx, '--black');