function [Y1, cm2, cmh2] = linearfft(y, classind, class, s)
L = length(y);
xx = abs(fft(y, [], 1));
freq = (1:L/2)/L; 

% Plot single-sided amplitude spectrum.
X1 = freq;
Y1 = xx(1:floor(L/2), :);

% direct SVD/PCA, then Kmeans using cityblock distance.
[coeff, score] = princomp(Y1', 'econ');
ggg = kmeans(score(:, 1:2), 2, 'Distance', 'correlation', 'Display','final', 'replicates', 10);
cm2 = confusionmat(class(classind), ggg);
cmh2 = condentropy(cm2);

% Create figure
%figure1 = figure('XVisual',...
%    '0x22 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
%axes('Parent',figure1,'FontWeight','bold','FontSize',16);
%box('on');
%hold('all');

% Create title
%title(strcat('Amplitude Specturm of ', s),'FontWeight','bold',...
%    'FontSize',16);

% Create xlabel
%xlabel('Frequency (cycles/bin)','FontWeight','bold','FontSize',16);

% Create ylabel
%ylabel('|Y(f)|','FontWeight','bold','FontSize',16);

% Create stem
%stem(X1,Y1,'MarkerSize',3,'LineStyle','none','Color',[0 0 1]);

%print 
