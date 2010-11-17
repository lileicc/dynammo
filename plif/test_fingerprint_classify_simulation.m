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
[group, fp, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 4, 'MaxIter', 100, 'Class', class);

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

[group, feature, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 6, 'MaxIter', 100, 'Class', class);

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



%%
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


class = [1 1 1 2 2];

%[group, entrop, P, D, mu0, feature] = fingerprint_classify(X, 'Hidden', 6, 'MaxIter', 50, 'Class', class);
[group, feature, entrop, P, D, mu0, model] = fingerprint_classify(X, 'Hidden', 6, 'MaxIter', 100, 'Class', class);
[u_k, V_k, P_k] = forward(X, model);
[Ex] = cell2mat(backward(u_k, V_k, P_k, model));
z0 = model.mu0;
zz = [];
for i = 1:500
  zz = [zz, z0];
  z0 = model.A * z0;
end
figure;
for i = 1:length(z0);
  subplot(length(z0), 1, i);
  plot(zz(i, :));
end

z0 = mu0;
zzz = [];
for i = 1:500
  zzz = [zzz, z0];
  z0 = D .* z0;
end
figure;
for i = 1:length(z0);
  subplot(length(z0), 1, i);
  plot(angle(zzz(i, :)));
end

figure; 
plot3(1:500, real(zzz(3, :)), imag(zzz(3, :)));
plot3(1:500, real(zzz(7, :)), imag(zzz(7, :)));

figure;
%colormap(autumn);
quiverc(1:size(P, 2), 1:size(P,1), real(P), imag(P));
axis ij
xlim([0.5, size(P,2)+0.5]);
ylim([0.5, size(P,1)+0.5]);
set(gca,'YTick', 1:size(P,1), 'YTickLabel',1:size(P,1));
xa = xlim;
ya = ylim;
for i = 1.5 : 1: (size(P, 1)-0.5)
  line(xa, [i, i],'LineStyle', '--', 'Color', 'black', 'LineWidth', 1);
end
for i = 1.5 : 1: (size(P, 2)-0.5)
  line([i, i], ya,'LineStyle', '--', 'Color', 'black', 'LineWidth', 1);
end



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

% plot the magnitude
figure;
imagesc(Q);
set(gca,'YTick', 1:5, 'YTickLabel',1:5);
set(gca,'XTick', 1:3, 'XTickLabel', 1:3);
colormap(autumn);
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
line([2.5, 2.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');

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


% svd version of the fingerprints
[uu, ss, vv] = svd(Q, 'econ');
fp_svd = uu * ss;
figure; 
imagesc(fp_svd(:, 1:2));
set(gca,'YTick', 1:5, 'YTickLabel',1:5);
set(gca,'XTick', 1:2, 'XTickLabel',['FP1'; 'FP2']);
colormap(autumn);
xa = xlim;
line(xa, [1.5, 1.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [2.5, 2.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [3.5, 3.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line(xa, [4.5, 4.5],'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
line([1.5, 1.5], [0.5, 5.5], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');

% FFT
figure; 
for i = 1 : size(X, 1)
  a = fft(X(i, :));
  subplot(5, 1, i);
  plot(500./(0:499), abs(a));
end
a = fft(X(4, :));
plot(500./(0:499), abs(a));
xlabel