clear;
load('../motion35-amc-labeled.mat');
classind = [find(classlabel==1); find(classlabel==2)];
% 1 = walking, 2 = running
classind = classind(randperm(length(classind)));

% build data in X
% rfoot
X52 = mocap35_dim{52}(:, classind);


%X = X52';
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
cm_clds = confusionmat(trueclass, pred1);
ce_clds = condentropy(trueclass, pred1);

[coeff1, score1] = princomp(features1(:, 1:end), 'econ');
pred_clds_pca = kmeans(score1(:,1:2), 2, 'replicates', 10);
cm_clds_pca = confusionmat(trueclass, pred_clds_pca);
ce_clds_pca = condentropy(trueclass, pred_clds_pca);


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
pred_fft = kmeans(score_ft(:, 1:2), 2, 'Display','final', 'replicates', 10);
cm_fft = confusionmat(trueclass, pred_fft);
ce_fft = condentropy(cm_fft);

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
[f_dtw, cm_dtw, ce_dtw] = dtw_cluster(X', trueclass);

