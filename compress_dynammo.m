function [errors, ratios, comp_data] = compress_dynammo(X, varargin)
% compression using dynammo learned features,
% there are three compression methods: fixed hop, adaptive, and optimal by
% dynamic programming.
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
% Optional Args:
%   'method': followed by a string in {'optimal', 'fixed', 'adaptive'}
%   'MaxCompressionRatio': followed by a number in [0, 1] denoting the
%   maximum space preserved. (only used in 'optimal')
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
%   comp_data: compressed data

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

a = find(strcmp('method', varargin));
if (isempty(a) || strcmp('optimal', varargin{a+1}))
  % dynammo_optimal compression using dynamic programming
  a = find(strcmp('MaxCompressionRatio', varargin));
  if (isempty(a))
    maxRatio = 1;
  else
    maxRatio = varargin{a+1};
  end
  thresh = ceil(N*maxRatio);
  errs = zeros(N + 1, thresh);
  last = zeros(N + 1, thresh);
  mark = zeros(1, N + 1);
  mark(1) = 1;
  mu = Ex{1};
  baseerr = zeros(N + 1, 1);
  curerr = 0;
  for i = 1 : N
    baseerr(i) = norm(X(:, i) - real(C * Ex{i}), 'fro') ^ 2;
    curerr = curerr + baseerr(i);
  end
  errs(1, 1) = baseerr(1);
  
  for ss = 1 : N
    mu = Ex{ss};
    tt = ss + 1;
    curerr = 0;
    while (tt <= N + 1)
      mu = A * mu;
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
        curerr = curerr + norm(X(:, tt) - real(C * mu), 'fro')^2;
      end
      tt = tt+1;
    end
  end
  H = length(model.mu0);
  ORIGINAL = numel(X) + 2;
  errors = errs(N+1, :);
  OVERHEAD = numel(model.mu0) + numel(model.A) + numel(model.C) + 4;
  K = length(errors);
  store_space = ((1:K) - 1) .* (H + 1) + OVERHEAD;
  ratios = ORIGINAL ./ store_space;  
  comp_data = cell(K, 1);
  for i = 1 : K
    TAG = 3;
    j = i;
    tt = N + 1;
    temp = [];
    while (last(tt, j) > 1)
      tt = last(tt, j);
      temp = [temp, tt];
    end
    comp_data{i} = zeros(store_space, 1);
    comp_data{i}(1:(4+H)) = [TAG, N, M, H, model.mu0'];
    ss = 5 + H;
    for tt = temp(end:(-1):1)
      comp_data{i}(ss) = tt;
      comp_data{i}((ss+1) : (ss+H)) = Ex{tt};
      ss = ss + H + 1;
    end
  end
elseif (strcmp('adaptive', varargin{a+1}))
  % dynammo_adaptive compression
  
  
elseif (strcmp('fixed', varargin{a+1}))
  % dynammo_fixed compression using fixed hop
  
end
  

