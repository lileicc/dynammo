% test scalability of Kalman Finger Printing (PLF).
clear;
data = load('chlorine_level_data_cl2fullLarge.dat');

X = data';
N = size(X, 2);
M = size(X, 1);
candM = ceil(M ./ 10 .* [1:10]);
H = 15;
GapTick = ceil(N / 10);
Trials = GapTick : GapTick : N;
time = zeros(length(candM), length(Trials));
i = 1;
j = 1;
for M = candM
  i = 1;
  for LEN = Trials
    tic;
    fingerprint(X(1:M, 1:LEN), 'Hidden', H, 'MaxIter', 20);
    time(j, i) = toc;
    i = i + 1;
  end
  j = j + 1;
end

figure;
plot(Trials, time(1, :), 'DisplayName', '17 sequences');
hold all;
plot(Trials, time(3, :), 'DisplayName', '50 sequences');
plot(Trials, time(5, :), 'DisplayName', '83 sequences');
plot(Trials, time(7, :), 'DisplayName', '117 sequences');
plot(Trials, time(9, :), 'DisplayName', '150 sequences');
legend('show', 'Location', 'Best');
xlabel('sequence length (ticks)');
ylabel('wall clock time(s)');

figure;
hold all;
plot(candM, time(:, 10), 'DisplayName', 'length=4310');
plot(candM, time(:, 7), 'DisplayName', 'length=3017');
plot(candM, time(:, 4), 'DisplayName', 'length=1724');
xlabel('# of sequences');
ylabel('wall clock time(s)');
