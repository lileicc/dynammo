t = 1 : 300;
f1 = 1 / 100;
u = [t / 150 - 1; t / 150 - 1];
%u = zeros(2, 300);
x = u + [sin(2 * pi * f1 * t); cos(2 * pi * f1 * t)];
f2 = 1 / 80;
y = x + [sin(2 * pi * f2 * t); cos(2 * pi * f2 * t)];

scatter(u(1,:), u(2,:));
hold all;
scatter(x(1, :), x(2,:));
scatter(y(1,:), y(2,:));

observed = true(3, 300);
%observed(2, 200:300) = false;
observed(3, 200:300) = false;

data = [u; x; y];
bones = [1, 2, 1; 2, 1, 1; 2, 3, 1; 3, 2, 1];
[model, Xhat, LL] = learn_lds_dynammop_bone_newton_direct(data, 'Bone', bones, 'MaxIter', 100, 'Hidden', 4, 'Observed', observed);
save('test_simulated_solar_multi_bone_direct.mat');

[model, Xhat, LL] = learn_lds_dynammop_bone_newton(data, 'Bone', bones, 'MaxIter', 100, 'Hidden', 4, 'Observed', observed);
save('test_simulated_solar_multi_bone.mat');

[model, Xhat, LL] = learn_lds_dynammop(data, 'MaxIter', 100, 'Hidden', 4, 'Observed', observed);
save('test_simulated_solar_fly.mat');
%% test basic learn_lds
[model, LL] = learn_lds(data, 'Bone', bones, 'MaxIter', 1000, 'Hidden', 4, 'Observed', observed);
[mu, V, P] = forward(data, model);
[Ez, Ezz, Ez1z] = backward(mu, V, P, model);
plot((model.C * cell2mat(Ez))' - data');

