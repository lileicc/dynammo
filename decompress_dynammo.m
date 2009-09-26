function [X] = decompress_dynammo(data)
% decompress from compressed data using dynammo learned features,
% there are three compression methods: fixed hop, adaptive, and optimal by
% dynamic programming.
% Args:
%   data: a vectore of numbers
%     the first number identifies the compression type: 
%     1 --> fixed
%     2 --> adaptive
%     3 --> optimal
%     the second, third, firth numbers are N, M, H
%     followed by mu0, A, C
%     the remaining will be 
%     for fixed, a number for hop, else, each H numbers for H
%     hidden variables for the time tick.
%     for adaptive and optimal, each H+1 numbers will be time tick, and H
%     hidden variables for the time tick.
%
% Returns:
%   X: decompressed data, M * N matrix
%
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
%

tp = data(1);
N = data(2);
M = data(3);
H = data(4);
tt = 5;
model = struct;
model.mu0 = reshape(data(tt : (tt + H - 1)), [], 1);
tt = tt + H;
model.A = reshape(data(tt : (tt + H * H -1)), [], H);
tt = tt + H * H;
model.C = reshape(data(tt : (tt + M * H - 1)), [], H);
tt = tt + M * H;
X = zeros(M, N);
z = model.mu0;
X(:, 1) = model.C * z;
if (tp == 1)
  if (tt <= length(data))
    hop = data(tt);
    tt = tt + 1;
  end
end
currentTick = 1;

for i = 2:N
  if (i > currentTick)
    if (tp == 1)      
        currentTick = currentTick + hop;
    else
      if (tt <= length(data))
        currentTick = data(tt);
        tt = tt + 1;
      end
    end
    if (tt <= length(data))
      z = data(tt : (tt + H - 1));
      tt = tt + H;
    end
  else
    z = model.A * z;
  end
  X(:, i) = model.C * z;
end