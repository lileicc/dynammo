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


model_train = learn_clds(X, 'Hidden', 4, 'MaxIter', 100);

score = abs(model_train.C);
figure;
imagesc(score);
export_fig('sin_data_clds_C.pdf', '-pdf');


%% 
clear;
N = 32;
t = (0:(N-1))';
f1 = 1 / 32;
f2 = 1 / 50;
t1 = 2 * pi * f1 * t;
t2 = 2 * pi * f2 * t;
%X = [sin(t1) + 0.2 * sin(t2), cos(t1) + 0.2 * sin(t2 + pi/4)];
X = sin(t1);
xft = fft(X);




X = X';
%A = diag([exp(2i*pi*f1), exp(-2i*pi*f1), exp(2i*pi*f2), exp(-2i*pi*f2)]);
A = diag([exp(-2i*pi*(0:(N-1)) / N)]);
mu0 = ones(N, 1);
Q0 = 0.001 * eye(N);
Q = 0.001 * eye(N);
model_train = learn_clds(X, 'model.A', A, 'model.mu0', mu0, 'Hidden', N, 'MaxIter', 500);
model_train = learn_clds(X, 'model.A', A, 'Hidden', N, 'MaxIter', 500);
model_train = learn_clds(X, 'model.A', A, 'model.mu0', mu0, 'model.Q0', Q0, 'model.Q', Q, 'Hidden', N, 'MaxIter', 500);
%model_train = learn_clds(X, 'Hidden', 2, 'MaxIter', 100);

xft_sp = 2 * abs(xft) / N;
x_clds_sp = abs(model_train.C);

figure;
hold all;
plot((0:(N-1))/N, xft_sp);
plot((0:(N-1))/N, x_clds_sp);
legend('DFT', 'CLDS', 'Location', 'Best');
title('spectrum');
xlabel('frequency');

figure;
hold all;
bar((0:(N-1))/N, [xft_sp, x_clds_sp']);
%bar((0:(N-1))/N, x_clds_sp);
legend('DFT', 'CLDS', 'Location', 'Best');
title('spectrum');
xlabel('frequency');


%% test on the mocap data
clear;
load('../motion16-labeled.1.mat');
classind = [find(class==2); find(class==3)];
% 2 = walking, 3 = running
classind = classind(randperm(length(classind)));

%%%%%%%%%%%%%%

% figure;
% idx = 33:33;
% for i=idx
%     subplot(length(idx), 1, i - idx(1) + 1);
%     plot(motion_dim{i}(:, [33:38]));
% end

% build data in X
X9 = motion_dim{9}(:, classind);
X15 = motion_dim{15}(:, classind);
X18 = motion_dim{18}(:, classind);
X24 = motion_dim{24}(:, classind);
% rfoot.z
X33 = motion_dim{33}(:, classind);


X = X33';
trueclass = class(classind);
%[group, features, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 8, 'MaxIter', 100, 'Class', trueclass, 'IsotropicQ', 'IsotropicR', 'IsotropicQ0');
for i = 1 : 10
  model{i} = learn_clds(X, 'Hidden', i, 'MaxIter', 10000);
  features{i} = abs(model{i}.C);
  pred{i} = kmeans(features{i}, 2);
  cm{i} = confusionmat(trueclass, pred{i});
  ce{i} = condentropy(trueclass, pred{i});
  [coeff, score] = princomp(features{i}, 'econ');
  pred1 = kmeans(score(:,1:2), 2);
  cm1 = confusionmat(trueclass, pred1);
end
figure;
hold all;
%scatter(features(trueclass==2, 1), features(trueclass==2, 2));
%scatter(features(trueclass==3, 1), features(trueclass==3, 2));
scatter(score(trueclass==2, 1), score(trueclass==2, 2));
scatter(score(trueclass==3, 1), score(trueclass==3, 2));
axis equal;
xlabel('FP1');
ylabel('FP2');

figure; plot(D, '*');
hold on;
plot(cos([1:100]*2*pi/100), sin([1:100]*2*pi/100), 'black')
axis equal;

figure;
[nouse, tmpidx] = sort(trueclass);
imagesc(features(tmpidx, 1:2));
set(gca,'XTick', 1:2, 'XTickLabel',['FP1'; 'FP2']);
colormap(grayColor);

figure; 
imagesc(abs(P(tmpidx, 2:2:6)));

%% for ploting data
tt = [15, 22, 45, 38, 8]; 
m = length(tt);
for i = 1:m;  
  subplot(m, 1, i); 
  plot(X(:, tt(i))); 
end

%% for baseline: PCA + kmean
[coeff, score] = princomp(X);
figure;
hold all;
scatter(score(trueclass==2, 1), score(trueclass==2, 2), 'bo');
scatter(score(trueclass==3, 1), score(trueclass==3, 2), 'r*');
axis equal;
xlabel('PC1');
ylabel('PC2');

