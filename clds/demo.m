% demo the learning of clds
% 
clear;
disp('Learning CLDS demo');
N = 500;
t = (1:N)';
a = zeros(N, 2);
f = 1 / 100;
b = [sin(2 * pi * f * t) cos(2 * pi * f * t) sin(2 * pi * f * t + pi/6)];

f1 = 1 / 110;
f2 = 1 / 30;
t1 = 2 * pi * f1 * t;
t2 = 2 * pi * f2 * t;
c = [sin(t1) + 0.2 * sin(t2), cos(t1) + 0.2 * sin(t2 + pi/4)];

X = [b c]';
M = size(X, 1);
figure;
set(gcf, 'defaultlinelinewidth', 1);
for i = 1:M
  subplot(M, 1, i);
  hold on;
  plot(X(i, :));
  if (i > 1)
    plot(X(1, :), 'black--');
  end
end

model_train = learn_clds(X, 'Hidden', 6, 'MaxIter', 100);
score = abs(model_train.C);
figure;
imagesc(score);
title('features from CLDS');
disp('mu0:');
disp(model_train.mu0);
disp('Q0:');
disp(model_train.Q0);
disp('A:');
disp(model_train.A);
disp('Q:');
disp(model_train.Q);
disp('C:');
disp(model_train.C);
disp('R:');
disp(model_train.R);

save('demo.mat');