% testing the fingerprint clustering on mocap data

clear;
load('motion16-labeled.1.mat');
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
[group, entrop, P, D, mu0, features, components] = fingerprint_classify(X, 'Hidden', 7, 'MaxIter', 100, 'Class', trueclass, 'IsotropicQ', 'IsotropicR', 'IsotropicQ0');
figure;
hold all;
scatter(features(trueclass==2, 1), features(trueclass==2, 2));
scatter(features(trueclass==3, 1), features(trueclass==3, 2));

figure; plot(D, '*');
hold on;
plot(cos([1:100]*2*pi/100), sin([1:100]*2*pi/100), 'black')
axis equal;

figure;
[nouse, tmpidx] = sort(group, 'descend');
imagesc(features(tmpidx, 1:2));

figure; 
imagesc(abs(P(tmpidx, 2:2:6)));