% sample code to generate features on mocap data and compare with other
% methods
% use with caution
% modify only when you understand each step

clear;
DATAFILE='../motion16-labeled.1.mat';

%% test on the mocap data
load(DATAFILE);
classind = [find(class==2); find(class==3)];
% 2 = walking, 3 = running
% random permutation
classind = classind(randperm(length(classind)));
trueclass = class(classind);

% build data in X
X9 = motion_dim{9}(:, classind);  %lhipjoint.z
X15 = motion_dim{15}(:, classind);  %ltibia.z
X18 = motion_dim{18}(:, classind);  %lfoot.z
X24 = motion_dim{24}(:, classind);  %rhipjoint.z
X33 = motion_dim{33}(:, classind);  % rfoot.z


X = X33';
%X = X33(1:50, :)';
%X = X33(1:100, :)';

%[model_clds, LL] = learn_clds(X, 'Hidden', 4, 'MaxIter', 10000);
[model_clds, LL] = learn_clds(X, 'Hidden', 8, 'MaxIter', 100);
features_clds = abs(model_clds.C);
%pred = kmeans(features_clds, 2, 'distance', 'correlation', 'replicates', 10);
pred = kmeans(features_clds, 2, 'replicates', 10);
cm_clds = confusionmat(trueclass, pred);
ce_clds = condentropy(trueclass, pred);

disp('CLDS conditional entropy');
disp(ce_clds);

[coeff1, score1] = princomp(features_clds(:, 1:end), 'econ');

figure;
hold all;
scatter(score1(trueclass==2, 1), score1(trueclass==2, 2), 'd', 'DisplayName', 'walking');
scatter(score1(trueclass==3, 1), score1(trueclass==3, 2), 'r*', 'DisplayName', 'running');
legend('show', 'Location', 'best');
export_fig 'scatter-mocap16-rfootz-clds.pdf' '-pdf'


%% for baseline: PCA + kmean
[coeff, score] = princomp(X, 'econ');
figure;
hold all;
scatter(score(trueclass==2, 1), score(trueclass==2, 2), 'bd');
scatter(score(trueclass==3, 1), score(trueclass==3, 2), 'r*');
axis equal;
xlabel('PC1');
ylabel('PC2');
export_fig 'scatter-mocap16-rfootz-pca.pdf' '-pdf'

%pred_pca = kmeans(score(:, 1:2), 2, 'distance', 'correlation', 'replicates', 10);
pred_pca = kmeans(score(:, 1:2), 2, 'replicates', 10);
cm_pca = confusionmat(trueclass, pred_pca);
ce_pca = condentropy(trueclass, pred_pca);

disp('PCA+Kmeans conditional entropy');
disp(ce_pca);

%% FFT
xft = fft(X');
xft = xft';
[coeff_ft, score_ft] = princomp(abs(xft));
%ggg = kmeans(score_ft(:, 1:2), 2, 'Distance', 'correlation', 'Display','final', 'replicates', 10);
ggg = kmeans(score_ft(:, 1:2), 2, 'Display','final', 'replicates', 10);
cm_fft = confusionmat(trueclass, ggg);
ce_fft = condentropy(cm_fft);

figure;
hold all;
scatter(score_ft(trueclass==2, 1), score_ft(trueclass==2, 2), 'd', 'DisplayName', 'walking');
scatter(score_ft(trueclass==3, 1), score_ft(trueclass==3, 2), 'r*', 'DisplayName', 'running');
legend('show', 'Location', 'best');
export_fig 'scatter-mocap16-rfootz-fft.pdf' '-pdf'

disp('FFT+Kmeans conditional entropy');
disp(ce_fft);


%% LDS
addpath('../dynammo');
[model_lds, LL_lds ] = learn_lds(X, 'Hidden', 8, 'MaxIter', 100);
[coeff_lds, score_lds] = princomp(model_lds.C, 'econ');
%pred_lds = kmeans(score_lds(:,1:2), 2, 'distance', 'correlation', 'Display','final', 'replicates', 10);
pred_lds = kmeans(score_lds(:,1:2), 2, 'Display','final', 'replicates', 10);
cm_lds= confusionmat(trueclass, pred_lds);
ce_lds = condentropy(trueclass, pred_lds);

figure;
hold all;
scatter(score_lds(trueclass==2, 1), score_lds(trueclass==2, 2), 'bd');
scatter(score_lds(trueclass==3, 1), score_lds(trueclass==3, 2), 'r*');

export_fig 'scatter-mocap16-rfootz-lds.pdf' '-pdf'
disp('LDS conditional entropy');
disp(ce_fft);

%% DTW
[f_dtw, cm_dtw, cmh_dtw] = dtw_cluster(X', trueclass);
