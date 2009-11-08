
clear;
load 'data_global_new.mat';

%X = data{22}(:, 4:96)';
%X = data{1}(:, 4:96)';

N = length(data);
time_fast = zeros(N, 1);
time_basic = zeros(N, 1);
for i = 1 : N %[15, 32, 45] %1 : N
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
axis square;
xx = xlim;
xlim([0, xx(2)]);
ylim([0, xx(2)]);
hold on;
line([0, 0], [xx(2), xx(2)], 'Color', 'black');
xlabel('PLiF-basic');
ylabel('PLiF');

figure;
ttt = [time_basic([15, 32, 45]) time_fast([15, 32, 45])];
bar(ttt)
