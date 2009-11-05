filenames = dir('*.csv');

data = [];
for i = 1:length(filenames)
  [U1 F1] = textread(filenames(i).name, '%s %d', 'headerlines', 1);
  data = [data; F1'];
end

X = data;
clear data;
N = size(X, 2);
M = size(X, 1);
candM = 10;
H = 10;
GapTick = ceil(N / 20);
Trials = N : (-GapTick) : GapTick;
time = zeros(length(candM), length(Trials));
time_slow = zeros(length(candM), length(Trials));
i = 1;
j = 1;
for M = candM
  i = 1;
  for LEN = Trials
    fprintf('@%d, %d \n', i, j);
    
    tic;
    fingerprint(X(1:M, 1:LEN), 'Hidden', H, 'MaxIter', 20);
    time_slow(j, i) = toc;
    
    clear ans;
    pack;
    
    tic;
    fingerprint(X(1:M, 1:LEN), 'Hidden', H, 'MaxIter', 20, 'Fast');
    time(j, i) = toc;   
    
    clear ans;
    pack;
    
    i = i + 1;
  end
  j = j + 1;
end

figure;
plot(Trials, time(1, :), 'LineWidth', 2, 'DisplayName', '10 sequences');
legend('show', 'Location', 'Best');
xlabel('sequence length (ticks)');
ylabel('wall clock time (s)');
xlim([0,N+1]);


% figure;
% hold all;
% plot(candM, time(:, 10), 'DisplayName', 'length=4310');
% plot(candM, time(:, 7), 'DisplayName', 'length=3017');
% plot(candM, time(:, 4), 'DisplayName', 'length=1724');
% xlabel('# of sequences');
% ylabel('wall clock time(s)');
