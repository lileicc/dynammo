t = 1 : 300;
f1 = 1 / 100;
u = [t / 150 - 1; t / 150 - 1];
x = u + [sin(2 * pi * f1 * t); cos(2 * pi * f1 * t)];
f2 = 1 / 80;
y = x + [sin(2 * pi * f2 * t); cos(2 * pi * f2 * t)];

scatter(u(1,:), u(2,:));
hold all;
scatter(x(1, :), x(2,:));
scatter(y(1,:), y(2,:));

observed = true(3, 300);
observed(2, 200:300) = false;

data = [u; x; y];
bones = [1, 2, 1; 2, 1, 1; 2, 3, 1; 3, 2, 1];
[model, Xhat, LL] = learn_lds_dynammop_bone_newton(data, 'Bone', bones, 'MaxIter', 100, 'Hidden', 3, 'Observed', observed);


