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
   
    tic;
    fingerprint(X(1:M, 1:LEN), 'Hidden', H, 'MaxIter', 20, 'Fast');
    time(j, i) = toc;   
    
    i = i + 1;
  end
  j = j + 1;
end

figure;
plot(Trials, time(1, :), 'LineWidth', 2, 'DisplayName', '17 sequences');
hold all;
plot(Trials, time(3, :), '-+', 'LineWidth', 2, 'DisplayName', '50 sequences');
plot(Trials, time(5, :), '-o', 'LineWidth', 2, 'DisplayName', '83 sequences');
plot(Trials, time(7, :), '--', 'LineWidth', 2, 'DisplayName', '117 sequences');
plot(Trials, time(9, :), '-.', 'LineWidth', 2, 'DisplayName', '150 sequences');
legend('show', 'Location', 'Best');
xlabel('sequence length (ticks)');
ylabel('wall clock time (s)');
xlim([0,4500]);

figure;
hold all;
plot(candM, time_slow(:, 10), '--', 'LineWidth', 2, 'DisplayName', 'PLiF-basic');
plot(candM, time(:, 10), 'LineWidth', 2, 'DisplayName', 'PLiF');
xlabel('# of sequences');
ylabel('wall clock time (s)');
legend('show', 'Location', 'Best');
xlim([0,170]);

% figure;
% hold all;
% plot(candM, time(:, 10), 'DisplayName', 'length=4310');
% plot(candM, time(:, 7), 'DisplayName', 'length=3017');
% plot(candM, time(:, 4), 'DisplayName', 'length=1724');
% xlabel('# of sequences');
% ylabel('wall clock time(s)');
