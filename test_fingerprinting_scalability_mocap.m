
clear;
load 'data_global_new.mat';

%X = data{22}(:, 4:96)';
%X = data{1}(:, 4:96)';

N = length(data);
time_fast = zeros(N, 1);
time_basic = zeros(N, 1);
for i = [15, 32, 45] %1 : N
  X = data{i}(:, 4:96)';
  
  tic;
  fingerprint(X, 'Hidden', 15, 'MaxIter', 20, 'Fast');
  time_fast(i) = toc;
  
  tic;
  fingerprint(X, 'Hidden', 15, 'MaxIter', 20);
  time_basic(i) = toc;
end

figure;
scatter(time_basic, time_fast);
axis equal;
hold on;
xx = xlim;
plot(xx, xx, '--black');
xlabel('PLiF-basic');
ylabel('PLiF-basic');

ttt = [time_basic([15, 32, 45]) time_fast([15, 32, 45])];
bar(ttt)
