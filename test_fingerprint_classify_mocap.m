% generate a sample with two frequencies
%
N = 300;
t = (1:N)';
a = zeros(N, 2);
f = 1 / 100;
t1 = 2 * pi * f * t;
b = [cos(t1) sin(t1) sin(t1)+cos(t1)];
f1 = 1 / 110;
%t2 = 2 * pi * f1 * t + pi / 4 + t1;
t2 = 2 * pi * f1 * t;
c = [0.5 * cos(t2) 0.5 * sin(t2)];
X = [b c];
M = size(X, 2);
figure;
for i = 1:M
    subplot(M, 1, i);
    plot(X(:, i));
end

class = [1 1 1 2 2];

[group, entrop, P, D, mu0] = fingerprint_classify(X, 'Hidden', 4, 'Iteration', 100, 'Class', class);

HIDDEN = length(D);

ind = find(abs(imag(D)) > 1e-10); % assume 0==1e-10
num_real = HIDDEN;
% further eliminate the conjugate ones
if (~ isempty(ind)) 
    num_real = ind(1) - 1;
    %D = D(subind);

    Q = abs(P(:, (num_real+1):2:end)); %only the magnitude
end

figure;
imagesc(Q);