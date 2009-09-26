function [errors, ratios, comp_data] = compress_dynammo(X, varargin)
% compression using dynammo learned features,
% there are three compression methods: fixed hop, adaptive, and optimal by
% dynamic programming.
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
% Optional Args:
%   'method': followed by a string in {'optimal', 'fixed', 'adaptive'}
%   'MaxSpaceRatio': followed by a number in [0, 1] denoting the
%   maximum space preserved. (only used in 'optimal', and 'adaptive')
%   'Hop': followed  by a number used in 'fixed'.
%   'model', followed by a struct with the following attributes:
%     A: transition matrix (also named as F in papers), H * H
%     C: transmission  matrix(also called projection matrix, named as G
% in papers), M * H
%     Q: transition covariance, H * H 
%     R: transmission covariance, M * M
%     mu0: initial states (also named as z_0 sometimes), H * 1
%     Q0: initial covariance, H * H
%
% Returns:
%   errors: squared error
%   ratios: compression ratio
%   comp_data: compressed data, a cell array, each with a vector of numbers
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
% $Author$@cs.cmu.edu
% $Date$
% $Rev$
%

N = size(X, 2);
M = size(X, 1);

a = find(strcmp('model', varargin));
if (isempty(a))
  [model, X] = learn_lds_dynammo(X, varargin{:});
  [Ex] = backward(forward(X, model), model);
else 
  model = varargin{a+1};
  if (iscell(varargin{a+2})) 
    Ex = varargin{a+2};
  else
    [Ex] = backward(forward(X, model), model);
  end
end

baseerr = zeros(N + 1, 1);
for i = 1 : N
  baseerr(i) = norm(X(:, i) - real(C * Ex{i}), 'fro') ^ 2;
end

H = length(model.mu0);
ORIGINAL = numel(X) + 2;

a = find(strcmp('method', varargin));
if (isempty(a) || strcmp('optimal', varargin{a+1}))
  % dynammo_optimal compression using dynamic programming
  a = find(strcmp('MaxSpaceRatio', varargin));
  if (isempty(a))
    maxRatio = 1;
  else
    maxRatio = varargin{a+1};
  end
  thresh = ceil(N*maxRatio);
  errs = zeros(N + 1, thresh);
  last = zeros(N + 1, thresh);
  errs(1, 1) = baseerr(1);
  
  for ss = 1 : N
    mu = Ex{ss};
    tt = ss + 1;
    curerr = 0;
    while (tt <= N + 1)
      mu = model.A * mu;
      if (ss > 1)
        i = 2;
      else
        i = 1;
      end      
      while ((i <= ss) && (i <= thresh))
        tmperr = curerr + baseerr(tt) + errs(ss, i);
        if ((last(tt, i+1) == 0) || (tmperr < errs(tt, i + 1)))
          errs(tt, i+1) = tmperr;
          last(tt, i+1) = ss;
        end
        i = i + 1;
      end
      if (tt <= N)
        curerr = curerr + norm(X(:, tt) - real(model.C * mu), 'fro')^2;
      end
      tt = tt+1;
    end
  end
  
  errors = errs(N+1, :);
  OVERHEAD = numel(model.mu0) + numel(model.A) + numel(model.C) + 4;
  K = length(errors);
  store_space = ((1:K) - 1) .* (H + 1) + OVERHEAD;
  ratios = ORIGINAL ./ store_space;  
  comp_data = cell(K, 1);
  TAG = 3;
  for i = 1 : K
    j = i;
    tt = N + 1;
    temp = [];
    while (last(tt, j) > 1)
      tt = last(tt, j);
      temp = [temp, tt];
    end
    comp_data{i} = zeros(store_space, 1);
    comp_data{1}(1:OVERHEAD) = [TAG; N; M; H; model.mu0; reshape(model.A, [], 1); reshape(model.C, [], 1)];
    ss = OVERHEAD + 1;
    for tt = temp(end:(-1):1)
      comp_data{i}(ss) = tt;
      comp_data{i}((ss+1) : (ss+H)) = Ex{tt};
      ss = ss + H + 1;
    end
  end
elseif (strcmp('adaptive', varargin{a+1}))
  % dynammo_adaptive compression
  a = find(strcmp('MaxSpaceRatio', varargin));
  if (isempty(a))
    maxRatio = 1;
  else
    maxRatio = varargin{a+1};
  end
  total_base_err = sum(baseerr);
  thresh = 4 * total_base_err / N / maxRatio;
  mark = zeros(N, 1);
  mark(1) = 1;
  mu = Ex{1};
  curerr = baseerr(1);
  ss = 1;
  totalerr = curerr;
  comp_data = cell(1, 1);
  TAG = 2;
  
  for i = 2 : N
    mark(i) = 0;
    mu = model.A * mu;
    temperr = norm(X(:, i) - real(model.C * mu), 'fro') ^ 2;
    curerr = curerr + temperr;
    if (curerr > thresh)      
      temperr = baseerr(i);
      curerr = baseerr(i);
      ss = ss + 1;
      mark(ss) = i;      
      mu = Ex{i};
    end
    totalerr = totalerr + temperr;
  end
    
  OVERHEAD = numel(model.mu0) + numel(model.A) + numel(model.C) + 4;
  errors = totalerr;
  store_space = (ss - 1) .* (H + 1) + OVERHEAD;
  ratios = ORIGINAL / store_space;
  comp_data{1} = zeros(store_space, 1);
  comp_data{1}(1:OVERHEAD) = [TAG; N; M; H; model.mu0; reshape(model.A, [], 1); reshape(model.C, [], 1)];
  tt = OVERHEAD + 1;
  for i = 2 : ss
    comp_data{1}(tt : (tt + H)) = [mark(i); Ex{mark(i)}];
    tt = tt + H + 1;
  end
elseif (strcmp('fixed', varargin{a+1}))
  % dynammo_fixed compression using fixed hop
  a = find(strcmp('Hop', varargin));
  if (isempty(a))
    hop = 4;
  else
    hop = varargin{a+1};
  end
  len = ceil(N/hop);
  TAG = 1;
  OVERHEAD = numel(model.mu0) + numel(model.A) + numel(model.C) + 5;
  store_space = (len - 1) .* H + OVERHEAD;
  ratios = ORIGINAL / store_space;
  comp_data{1} = zeros(store_space, 1);
  comp_data{1}(1:OVERHEAD) = [TAG; N; M; H; model.mu0; reshape(model.A, [], 1); reshape(model.C, [], 1); hop];
  tt = OVERHEAD + 1;
  for i = 2 : len
    comp_data{1}(tt : (tt + H - 1)) = Ex{(i-1)*hop + 1};
    tt = tt + H;
  end
end
  

