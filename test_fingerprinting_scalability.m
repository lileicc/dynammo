% test scalability of Kalman Finger Printing (PLF).
clear;
data = load('chlorine_level_data_cl2fullLarge.dat');

X = data';
N = size(X, 2);
M = size(X, 1);
candM = ceil(M ./ 5 .* [1:5]);
H = 15;
GapTick = ceil(N / 10);
Trials = 10 : GapTick : N;
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
