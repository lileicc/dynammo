% generate a sample with two frequencies
%
N = 500;
t = (1:N)';
a = zeros(N, 2);
f = 1 / 100;
t1 = 2 * pi * f * t;
b = [cos(t1) sin(t1) sin(t1)+cos(t1)];
f1 = 1 / 110;
%t2 = 2 * pi * f1 * t + pi / 4 + t1;
t2 = 2 * pi * f1 * t;
c = [0.5 * cos(t2) 0.5 * sin(t2)];
X = [b c]';
M = size(X, 1);
figure;
for i = 1:M
  subplot(M, 1, i);
  hold on;
  plot(X(i, :));
  if (i > 1)
    plot(X(1, :), 'black--');
  end
end

class = [1 1 1 2 2];

[group, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 4, 'MaxIter', 100, 'Class', class);

HIDDEN = length(D);

ind = find(abs(imag(D)) > 1e-10); % assume 0==1e-10
num_real = HIDDEN;
% further eliminate the conjugate ones
if (~ isempty(ind))
  num_real = ind(1) - 1;
  %D = D(subind);
  
  Q = abs(P(:, (num_real+1):2:end)); %only the magnitude
end

figure;
imagesc(Q);
colormap(flipud(gray));
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');


%% harmonics grouping example
N = 500;
t = (1:N)';
a = zeros(N, 2);
f1 = 1 / 100;
f2 = 1 / 70;
t1 = 2 * pi * f1 * t;
t2 = 2 * pi * f2 * t;
b = [sin(t1) , cos(t1), sin(t1) + cos(t1)];

f1 = 1 / 110;
f2 = 1 / 30;
t1 = 2 * pi * f1 * t;
t2 = 2 * pi * f2 * t;
c = [sin(t1) + 0.2 * sin(t2), cos(t1) + 0.2 * sin(t2 + pi/4)];

X = [b c]';
M = size(X, 1);
h = figure('position', [568   287   368   419]);
for i = 1:M
  subplot('position', [0.05, 1 - 0.15*i - 0.05 * i + 0.05, 0.92, 0.14]);
  hold on;
  ylim([-2, 2]);
  plot(X(i, :));
  if (i > 1)
    plot(X(1, :), 'black--');
  end
end

class = [1 1 1 2 2];

[group, entrop, P, D, mu0, feature] = fingerprint_classify(X, 'Hidden', 6, 'MaxIter', 100, 'Class', class);

periods = 2 * pi ./ angle(D);

HIDDEN = length(D);

ind = find(abs(imag(D)) > 1e-10); % assume 0==1e-10
num_real = HIDDEN;
% further eliminate the conjugate ones
if (~ isempty(ind))
  num_real = ind(1) - 1;
  %D = D(subind);
  
  Q = abs(P(:, (num_real+1):2:end)); %only the magnitude
end

figure;
imagesc(Q);
set(gca,'YTick', 1:5, 'YTickLabel',1:5);
set(gca,'XTick', 1:3, 'XTickLabel', 1:3);
colormap(flipud(gray));
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([2.5, 2.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');


figure;
imagesc(feature(:, 1:2));
set(gca,'YTick', 1:5, 'YTickLabel',1:5);
set(gca,'XTick', 1:2, 'XTickLabel',['FP1'; 'FP2']);
colormap(flipud(gray));
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
colormap(flipud(gray));
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
