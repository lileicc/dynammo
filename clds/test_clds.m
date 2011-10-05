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
b = [sin(2 * pi * (1/100) * t) cos(2 * pi * (1/100) * t) sin(2 * pi * (1/100) * t + pi/6)];

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


model_train = learn_clds(X, 'Hidden', 6, 'MaxIter', 100);

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
A = diag([exp(2i*pi*(0:(N-1)) / N)]);
mu0 = ones(N, 1);
Q0 = 0.001 * eye(N);
Q = 0.001 * eye(N);
model_train = learn_clds(X, 'model.A', A, 'model.mu0', mu0, 'Hidden', N, 'MaxIter', 500);
model_train = learn_clds(X, 'model.A', A, 'Hidden', N, 'MaxIter', 500);
[model_train, LL] = learn_clds(X, 'model.A', A, 'model.mu0', mu0, 'Hidden', N, 'MaxIter', 500);
[model_train, LL] = learn_clds(X, 'model.A', A, 'model.mu0', mu0, 'model.Q0', Q0, 'model.Q', Q, 'Hidden', N, 'MaxIter', 500);
%model_train = learn_clds(X, 'Hidden', 2, 'MaxIter', 100);
[model_train, LL] = learn_clds(X, 'model.A', A, 'model.mu0', mu0, 'Hidden', N, 'MaxIter', 1000, 'DiagQ0', 'DiagQ');


modelfft.mu0=mu0;
modelfft.A = A;
modelfft.C = xft.' / N;
modelfft.Q = Q;
modelfft.R = 1;
[model_train, LL] = learn_clds(X, 'Model', modelfft, 'Hidden', N, 'MaxIter', 100);

xft_sp = abs(xft) / N;
x_clds_sp = abs(model_train.C);

figure;
hold all;
plot((0:(N-1))/N, xft_sp);
plot((0:(N-1))/N, x_clds_sp);
legend('DFT', 'CLDS', 'Location', 'Best');
title('spectrum');
xlabel('frequency');

figure;
colormap colorGray;
bar((0:(N-1))/N, [xft_sp, x_clds_sp'], 0.98);
%bar((0:(N-1))/N, x_clds_sp);
legend('DFT', 'CLDS', 'Location', 'Best');
ylabel('spectrum');
xlabel('frequency');
xlim([0, 1]);

figure;
colormap colorGray;
semilogx(LL);
ylabel('log-likelihood');
xlabel('iteration');



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
X = X33(1:50, :)';
X = X33(1:100, :)';

trueclass = class(classind);
%[group, features, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 8, 'MaxIter', 100, 'Class', trueclass, 'IsotropicQ', 'IsotropicR', 'IsotropicQ0');
[model_train, LL] = learn_clds(X, 'Hidden', 4, 'MaxIter', 10000);
features = abs(model_train.C);
pred = kmeans(features, 2);
cm = confusionmat(trueclass, pred);
ce = condentropy(trueclass, pred);

% check PLiF
[Feature, P, D, mu0, zhat, model] = fingerprint(X, 'Hidden', 8, 'MaxIter', 100, 'Class', trueclass, 'IsotropicQ', 'IsotropicR', 'IsotropicQ0');
 
% use PLiF to initialize
model0.mu0 = mu0;
H = length(mu0);
model0.Q0 = eye(H) * 0.01;
model0.A = diag(D);
model0.Q = model0.Q0;
model0.C = P + rand(size(P));
model0.R = eye(size(X, 1)) * 0.01;
[model_train1, LL1] = learn_clds(X, 'Model', model0, 'Hidden', 8, 'MaxIter', 100);
features1 = abs(model_train1.C(:, 1:end));
%pred1 = kmeans(features1, 2, 'replicates', 10);
pred1 = kmeans(features1, 2, 'distance', 'correlation', 'replicates', 10, 'Options', statset('Display','final'));
cm1 = confusionmat(trueclass, pred1);
ce1 = condentropy(trueclass, pred1);

[coeff1, score1] = princomp(features1(:, 1:end), 'econ');
pred1_pca = kmeans(score1(:,1:2), 2, 'replicates', 10);
cm1_pca = confusionmat(trueclass, pred1_pca);
ce1_pca = condentropy(trueclass, pred1_pca);


figure;
hold all;
scatter(features1(trueclass==2, 1), features1(trueclass==2, 2));
scatter(features1(trueclass==3, 1), features1(trueclass==3, 2));

figure;
hold all;
scatter(score1(trueclass==2, 1), score1(trueclass==2, 2), 'd', 'DisplayName', 'walking');
scatter(score1(trueclass==3, 1), score1(trueclass==3, 2), 'r*', 'DisplayName', 'running');
legend('show', 'Location', 'best');
export_fig 'scatter-mocap-clds-100.pdf' '-pdf'

for i = 1 : 10
  model{i} = learn_clds(X, 'Hidden', i, 'MaxIter', 10000);
  features{i} = abs(model{i}.C);
  pred{i} = kmeans(features{i}, 2, 'replicates', 5);
  cm{i} = confusionmat(trueclass, pred{i});
  ce{i} = condentropy(trueclass, pred{i});
  %[coeff, score] = princomp(features{i}, 'econ');
  %pred1 = kmeans(score(:,1:2), 2);
  %cm1{i} = confusionmat(pred1, trueclass);
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
[coeff, score] = princomp(X, 'econ');
figure;
hold all;
scatter(score(trueclass==2, 1), score(trueclass==2, 2), 'bd');
scatter(score(trueclass==3, 1), score(trueclass==3, 2), 'r*');
axis equal;
xlabel('PC1');
ylabel('PC2');
pred_pca = kmeans(score(:, 1:2), 2, 'distance', 'correlation', 'replicates', 10);
cm_pca = confusionmat(trueclass, pred_pca);
ce_pca = condentropy(trueclass, pred_pca);


%% FFT
xft = fft(X');
xft = xft';
[coeff_ft, score_ft] = princomp(abs(xft));
figure;
hold all;
scatter(score_ft(trueclass==2, 1), score_ft(trueclass==2, 2), 'd', 'DisplayName', 'walking');
scatter(score_ft(trueclass==3, 1), score_ft(trueclass==3, 2), 'r*', 'DisplayName', 'running');
legend('show', 'Location', 'best');

% aditya's version
y = X';
L = length(y);
xx = abs(fft(y, [], 1));
freq = (1:L/2)/L; 

% Plot single-sided amplitude spectrum.
X1 = freq;
Y1 = xx(1:floor(L/2), :);

% direct SVD/PCA, then Kmeans using cityblock distance.
[coeff_ft, score_ft] = princomp(Y1', 'econ');
ggg = kmeans(score_ft(:, 1:2), 2, 'Distance', 'correlation', 'Display','final', 'replicates', 10);
%ggg = kmeans(score(:, 1:2), 2, 'Display','final', 'replicates', 10);
cm2 = confusionmat(trueclass, ggg);
cmh2 = condentropy(cm2);

figure;
hold all;
scatter(score_ft(trueclass==2, 1), score_ft(trueclass==2, 2), 'd', 'DisplayName', 'walking');
scatter(score_ft(trueclass==3, 1), score_ft(trueclass==3, 2), 'r*', 'DisplayName', 'running');
%legend('show', 'Location', 'best');
export_fig 'scatter-mocap-fft-100.pdf' '-pdf'


%% LDS
[model_lds, LL_lds ] = learn_lds(X, 'Hidden', 8, 'MaxIter', 100);
[coeff_lds, score_lds] = princomp(model_lds.C, 'econ');
pred_lds = kmeans(score_lds(:,1:2), 2, 'distance', 'correlation');
cm_lds= confusionmat(trueclass, pred_lds);
ce_lds = condentropy(trueclass, pred_lds);

figure;
hold all;
scatter(score_lds(trueclass==2, 1), score_lds(trueclass==2, 2), 'bd');
scatter(score_lds(trueclass==3, 1), score_lds(trueclass==3, 2), 'r*');

export_fig 'scatter-mocap-lds.pdf' '-pdf'


%% DTW
[f_dtw, cm_dtw, cmh_dtw] = dtw_cluster(X', trueclass);



%% test on the amc data
clear;
load('../motion35-amc-labeled.mat');
classind = [find(classlabel==1); find(classlabel==2)];
% 1 = walking, 2 = running
classind = classind(randperm(length(classind)));

% build data in X
% rfoot
X52 = mocap35_dim{52}(:, classind);


X = X52';
X = X52(1:60, :)';

trueclass = classlabel(classind);
%[group, features, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 8, 'MaxIter', 100, 'Class', trueclass, 'IsotropicQ', 'IsotropicR', 'IsotropicQ0');
[model_train, LL] = learn_clds(X, 'Hidden', 4, 'MaxIter', 100);
features = abs(model_train.C);
pred = kmeans(features, 2);
cm = confusionmat(trueclass, pred);
ce = condentropy(trueclass, pred);

% check PLiF
H = 4;
[Feature, P, D, mu0, zhat, model] = fingerprint(X, 'Hidden', H, 'MaxIter', 100, 'Class', trueclass, 'IsotropicQ', 'IsotropicR', 'IsotropicQ0');
 
% use PLiF to initialize
model0.mu0 = mu0;
H = length(mu0);
model0.Q0 = eye(H) * 0.01;
model0.A = diag(D);
model0.Q = model0.Q0;
model0.C = P;
model0.R = eye(size(X, 1)) * 0.01;
[model_train1, LL1] = learn_clds(X, 'Model', model0, 'Hidden', H, 'MaxIter', 100);
features1 = abs(model_train1.C(:, 3:end));
pred1 = kmeans(features1, 2, 'replicates', 10);
%pred1 = kmeans(features1, 2, 'distance', 'correlation', 'replicates', 10, 'Options', statset('Display','final'));
cm1 = confusionmat(trueclass, pred1);
ce1 = condentropy(trueclass, pred1);

[coeff1, score1] = princomp(features1(:, 1:end), 'econ');
pred1_pca = kmeans(score1(:,1:2), 2, 'replicates', 10);
cm1_pca = confusionmat(trueclass, pred1_pca);
ce1_pca = condentropy(trueclass, pred1_pca);


figure;
hold all;
scatter(score1(trueclass==1, 1), score1(trueclass==1, 2), 'd', 'DisplayName', 'walking');
scatter(score1(trueclass==2, 1), score1(trueclass==2, 2), 'r*', 'DisplayName', 'running');
legend('show', 'Location', 'best');
export_fig 'scatter-mocap35-52-clds.pdf' '-pdf'

%% FFT for amc data
y = X';
L = length(y);
xx = abs(fft(y, [], 1));
freq = (1:L/2)/L; 

% Plot single-sided amplitude spectrum.
X1 = freq;
Y1 = xx(1:floor(L/2), :);

% direct SVD/PCA, then Kmeans using cityblock distance.
[coeff_ft, score_ft] = princomp(Y1', 'econ');
%ggg = kmeans(score_ft(:, 1:2), 2, 'Distance', 'correlation', 'Display','final', 'replicates', 1);
ggg = kmeans(score_ft(:, 1:2), 2, 'Display','final', 'replicates', 10);
cm2 = confusionmat(trueclass, ggg);
cmh2 = condentropy(cm2);

figure;
hold all;
scatter(score_ft(trueclass==1, 1), score_ft(trueclass==1, 2), 'd', 'DisplayName', 'walking');
scatter(score_ft(trueclass==2, 1), score_ft(trueclass==2, 2), 'r*', 'DisplayName', 'running');
%legend('show', 'Location', 'best');
export_fig 'scatter-mocap-fft-100.pdf' '-pdf'

%% PCA + kmeans on amc data
[coeff, score] = princomp(X, 'econ');
figure;
hold all;
scatter(score(trueclass==1, 1), score(trueclass==1, 2), 'bd');
scatter(score(trueclass==2, 1), score(trueclass==2, 2), 'r*');
axis equal;
xlabel('PC1');
ylabel('PC2');
pred_pca = kmeans(score(:, 1:2), 2,  'replicates', 10);
cm_pca = confusionmat(trueclass, pred_pca);
ce_pca = condentropy(trueclass, pred_pca);


%% LDS
[model_lds, LL_lds ] = learn_lds(X, 'Hidden', 4, 'MaxIter', 100);
[coeff_lds, score_lds] = princomp(model_lds.C, 'econ');
pred_lds = kmeans(score_lds(:,1:2), 2);
cm_lds= confusionmat(trueclass, pred_lds);
ce_lds = condentropy(trueclass, pred_lds);

figure;
hold all;
scatter(score_lds(trueclass==2, 1), score_lds(trueclass==2, 2), 'bd');
scatter(score_lds(trueclass==3, 1), score_lds(trueclass==3, 2), 'r*');

export_fig 'scatter-mocap-lds.pdf' '-pdf'


%% DTW
[f_dtw, cm_dtw, cmh_dtw] = dtw_cluster(X', trueclass);

