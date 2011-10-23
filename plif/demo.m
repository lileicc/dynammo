% $Author: leili $@cs.cmu.edu
% $Date: 2011-10-08 22:30:16 -0400 (Sat, 08 Oct 2011) $
% $Rev: 335 $

% demo
N = 500;
t = (1:N)';
f = 1 / 100;
b = [sin(2 * pi * (1/100) * t) cos(2 * pi * (1/100) * t) sin(2 * pi * (1/98) * t + pi/6)];

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

class = [1 1 1 2 2];

[group, feature, entrop, P, D, mu0, model] = fingerprint_classify(X, 'Hidden', 6, 'MaxIter', 100, 'Class', class);

% ploting the LDS
figure;
imagesc(model.C);
set(gca,'YTick', 1:5, 'YTickLabel',1:5);
%set(gca,'XTick', 1:3, 'XTickLabel', 1:3);
colormap(autumn);
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([2.5, 2.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([3.5, 3.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([4.5, 4.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([5.5, 5.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([6.5, 6.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([7.5, 7.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');

figure;
imagesc(feature(:, 1:2));
set(gca,'YTick', 1:5, 'YTickLabel',1:5);
set(gca,'XTick', 1:2, 'XTickLabel',['FP1'; 'FP2']);
colormap(autumn);
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');

[coeff, score] = princomp(X);
figure;
imagesc(score(:, 1:2));
set(gca,'YTick', 1:5, 'YTickLabel',1:5);
set(gca,'XTick', 1:2, 'XTickLabel',['PC1'; 'PC2']);
colormap(autumn);
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');

save('demo.mat');
