%% compare against FFT

clear;
N = 32;
t = (0:(N-1))';
f1 = 1 / 32;
f2 = 1 / 16;
f3 = 1 / 128;
t1 = 2 * pi * f1 * t;
t2 = 2 * pi * f2 * t;
X = sin(2 * pi * f3 * t);

X = X';
H = N;

A = diag([exp(2i*pi*(0:(N-1)) / H)]); % set the frequency to be same as FFT
mu0 = ones(H, 1);
Q0 = 0.001 * eye(H);
Q = 0.001 * eye(H);
[model_train, LL] = learn_clds(X, 'model.A', A, 'model.mu0', mu0, 'model.Q0', Q0, 'model.Q', Q, 'Hidden', N, 'MaxIter', 500);

xft = fft(X);
xft_sp = abs(xft) / N;
x_clds_sp = abs(model_train.C); % spectrum for clds

figure('Position', [126, 184, 493, 218]);
set(gca, 'FontSize', 20, 'XTick', [0:0.2:1], 'YTick', [0:0.2:0.4]);
box on
colormap colorGray;
hold all;
ttt = x_clds_sp > 1e-4;
aaa = (0:(N-1))/N;
bar(aaa(ttt), [xft_sp(ttt)', x_clds_sp(ttt)'], 0.98);
ylabel('spectrum');
xlabel('frequency');
xlim([-0.02, 1]);
ylim([0, 0.55]);
legend('DFT', 'CLDS', 'Location', 'Best');
%export_fig 'clds-vs-fft.pdf' '-pdf'
saveas(gcf, 'clds-vs-fft.fig')



figure('Position', [126, 184, 493, 218]);
set(gca, 'FontSize', 18, 'XTick', [0:0.2:1], 'YTick', [0:0.1:0.4]);
box on
colormap colorGray;
hold all;
bar((0:(N-1))/N, [xft_sp', x_clds_sp'], 0.98);
ylabel('spectrum');
xlabel('frequency');
xlim([-0.02, 1]);
ylim([0, 0.42]);
legend('DFT', 'CLDS', 'Location', 'Best');
%export_fig 'clds-vs-fft3.pdf' '-pdf' %if you want to export to pdf, uncomment this 
saveas(gcf, 'clds-vs-fft3.fig')


figure;
colormap colorGray;
hold all;
bar((0:(N-1))/N, xft_sp);
bar([1/127.1904, -1/127.1904], x_clds_sp', 'g');
legend('DFT', 'CLDS', 'Location', 'Best');
ylabel('spectrum');
xlabel('frequency');

figure;
colormap colorGray;
semilogx(LL);
ylabel('log-likelihood');
xlabel('iteration');

%% spectrum by fft
clear;
N = 32;
t = (0:(N-1))';
f3 = 1 / 128;
X = sin(2 * pi * f3 * t);
X = X';
xx = sin(2*pi*f3*(N:128));

xft = fft(X);
xft_sp = abs(xft) / N;

figure('Position', [126, 184, 493, 218]);
set(gca, 'FontSize', 18, 'Box', 'on', 'XTick', [0:20:140], 'YTick', [-1:0.5:1]);
hold all;
plot(X, '.');
plot(32:128, xx, '--r', 'LineWidth', 2);
legend('observation', 'true signal', 'Location', 'Best');
%export_fig 'signal3.pdf' '-pdf'
saveas(gcf, 'signal3.fig');

figure('Position', [126, 184, 493, 218]);
colormap colorGray;
set(gca, 'FontSize', 18, 'Box', 'on', 'XTick', [0:0.2:1], 'YTick', [0:0.2:0.6]);
bar((0:(N-1))/N, xft_sp);
ylabel('spectrum');
xlabel('frequency');
xlim([-0.02, 1]);
ylim([0, 0.65]);
%export_fig 'fft-signal3.pdf' '-pdf'
saveas(gcf, 'fft-signal3.fig');



