function [Y1, cm2, cmh2] = dtw_cluster(y, class, s)

Y1 = zeros(size(y, 2), size(y, 2));
for i=1:size(y, 2) 
		i
	for j=1:size(y, 2) 	
			% for each pair
			[Dist,D,k,w]=dtw(y(:, i)',y(:, j)');
			Y1(i, j) = Dist/k;
	end
end


% direct SVD/PCA, then Kmeans using euclidean.
[coeff, score] = princomp(Y1', 'econ');
ggg = kmeans(score(:, 1:2), 2, 'Display','final', 'replicates', 10);
cm2 = confusionmat(ggg, class);
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
