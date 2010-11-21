%% generate random clds
clear;
[X, model_true] = sample_clds(100, 5, 2);
Y = real(X);
figure;
plot(Y');

model_train = learn_clds(Y, 'Hidden', 2, 'MaxIter', 1000);
model_train1 = learn_clds(Y, 'Hidden', 2, 'MaxIter', 1000);

Xhat = sample_clds(model_train, 100);
Yhat = real(Xhat);
figure;
subplot(2, 1, 1);
plot(Y');
subplot(2, 1, 2);
plot(Yhat');

%% generate sample sin data
clear;
N = 500;
t = (1:N)';
a = zeros(N, 2);
f = 1 / 100;
t1 = 2 * pi * f * t;
b = [sin(2 * pi * (1/100) * t) cos(2 * pi * (1/100) * t) sin(2 * pi * (1/98) * t + pi/6)];

f1 = 1 / 110;
f2 = 1 / 30;
t1 = 2 * pi * f1 * t;
t2 = 2 * pi * f2 * t;
c = [sin(t1) + 0.2 * sin(t2), cos(t1) + 0.2 * sin(t2 + pi/4)];

X = [b c]';
M = size(X, 1);
figure;
%set(gcf, 'defaultlinelinewidth', 1);
for i = 1:M
  subplot(M, 1, i);
  hold on;
  plot(X(i, :));
  if (i > 1)
    plot(X(1, :), 'black--');
  end
end
export_fig('sin_data.pdf', '-pdf');


model_train = learn_clds(X, 'Hidden', 4, 'MaxIter', 1000);

score = abs(model_train.C);
figure;
imagesc(score);
export_fig('sin_data_clds_C.pdf', '-pdf');