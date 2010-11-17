clear;
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
observed(2, 200:300) = false;
%observed(3, 200:300) = false;
%observed(3, 200:290) = false;

data = [u; x; y];
bones = [1, 2, 1; 2, 1, 1; 2, 3, 1; 3, 2, 1];
%bones = [1, 2, 0.8; 2, 1, 0.8; 2, 3, 0.5; 3, 2, 0.5];
[model_direct, Xhat_newton_direct, LL] = learn_lds_dynammop_bone_newton_direct(data, 'Bone', bones, 'MaxIter', 100, 'Hidden', 4, 'Observed', observed, 'PlotFun', @(x)plot(x'));
%save('test_simulated_solar_multi_bone_direct.mat');

[model_multi, Xhat_multi_bone, LL] = learn_lds_dynammop_bone_newton(data, 'Bone', bones, 'MaxIter', 100, 'Hidden', 4, 'Observed', observed, 'PlotFun', @(x)plot(x'));
%save('test_simulated_solar_multi_bone.mat');

[model_dynammop, Xhat_dynammop, LL] = learn_lds_dynammop(data, 'MaxIter', 100, 'Hidden', 4, 'Observed', observed, 'PlotFun', @plot);
%save('test_simulated_solar_fly.mat');

[model_dynammo, Xhat_dynammo, LL] = learn_lds_dynammo(data, 'MaxIter', 100, 'Hidden', 4, 'Observed', observed, 'PlotFun', @plot);
save('test_simulated_solar_5-6_missing200-300.mat');

%% test basic learn_lds
[model, LL] = learn_lds(data, 'Bone', bones, 'MaxIter', 1000, 'Hidden', 4, 'Observed', observed);
[mu, V, P] = forward(data, model);
[Ez, Ezz, Ez1z] = backward(mu, V, P, model);
plot((model.C * cell2mat(Ez))' - data');


%% plot
%%
Xhat = Xhat_dynammop;
Xhat = data;
figure;
plot(Xhat(5:6, :)');
x1 = 199.5;
x2 = 300;
y1 = min(ylim);
y2 = -1;
lv = [x1, y1; x1, y2];
lv2 = [x2, y1; x2, y2];
xlabel('frame');
%ylabel('bone length (m)'); 
legend('show');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);


Xhat = Xhat_newton_direct;
[bbb, bbv, bbs] = get_bones(Xhat, 'Dim', 2, 'Threshold', 1);
figure;
plot(bbs(:, 2, 3));
ylim([0, 1.2]);
x1 = 199.5;
x2 = 300;
y1 = min(ylim);
y2 = 0.5;
lv = [x1, y1; x1, y2];
lv2 = [x2, y1; x2, y2];
xlabel('frame');
ylabel('bone length (m)'); 
legend('show');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);