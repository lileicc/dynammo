

Y_multi_bone = Y;
%% plot bone length
% fly
Y_fly = Y;
[b, bv, bls] = get_bones(Y_fly);
figure;
hold on;
plot(bls(:, 25, 37), 'b', 'DisplayName', 'RELB-RUPA');
plot(bls(:, 25, 28), 'r--', 'DisplayName', 'RELB-RFRM');
plot(bls(:, 28, 38), 'm-.', 'DisplayName', 'RFRM-RWRB');
ylim([0 0.6]);
lv = [100, min(ylim);100, max(ylim)];
lv2 = [500, min(ylim);500, max(ylim)];
xlabel('frame');
ylabel('bone length (m)'); 
legend('show');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);


% multi_bone
%Y_fly = Y;
[b, bv, bls] = get_bones(Y_multi_bone);
figure;
hold on;
plot(bls(:, 25, 37), 'b', 'DisplayName', 'RELB-RUPA');
plot(bls(:, 25, 28), 'r--', 'DisplayName', 'RELB-RFRM');
plot(bls(:, 28, 38), 'm-.', 'DisplayName', 'RFRM-RWRB');
ylim([0 0.6]);
lv = [100, min(ylim);100, max(ylim)];
lv2 = [500, min(ylim);500, max(ylim)];
xlabel('frame');
ylabel('bone length (m)'); 
legend('show');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);


%original
[b, bv, bls] = get_bones(X);
figure;
hold on;
plot(bls(:, 25, 37), 'b', 'DisplayName', 'RELB-RUPA');
plot(bls(:, 25, 28), 'r--', 'DisplayName', 'RELB-RFRM');
plot(bls(:, 28, 38), 'm-.', 'DisplayName', 'RFRM-RWRB');
ylim([0 0.6]);
lv = [100, min(ylim);100, max(ylim)];
lv2 = [500, min(ylim);500, max(ylim)];
xlabel('frame');
ylabel('bone length (m)'); 
legend('show');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);


%% plot frame for 132_43 100 - 500, 
% for frame 282
%% plot frame
markers = X';
figure;
t = 282;
draw_skel(markers(t, :), lab.colheaders);
hold on;
xyz = [markers(t,109), markers(t,110), markers(t,111); markers(t,73), markers(t,74), markers(t,75); markers(t,82), markers(t,83), markers(t,84); markers(t,112), markers(t,113), markers(t,114); markers(t,115), markers(t,116), markers(t,117); markers(t,79), markers(t,80), markers(t,81)];
line(xyz(1:2, 1), xyz(1:2, 2), xyz(1:2, 3), 'Color', 'b', 'LineWidth', 6);
line(xyz(2:3, 1), xyz(2:3, 2), xyz(2:3, 3), 'Color', 'r', 'LineWidth', 6);
line(xyz(3:4, 1), xyz(3:4, 2), xyz(3:4, 3), 'Color', 'g', 'LineWidth', 6);
line(xyz(4:5, 1), xyz(4:5, 2), xyz(4:5, 3), 'Color', 'c', 'LineWidth', 6);
line(xyz([3,5], 1), xyz([3,5], 2), xyz([3,5], 3), 'Color', 'y', 'LineWidth', 6);
line(xyz(5:6, 1), xyz(5:6, 2), xyz(5:6, 3), 'Color', 'm', 'LineWidth', 6);
plot3(markers(t,73), markers(t,74), markers(t,75), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,82), markers(t,83), markers(t,84), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,112), markers(t,113), markers(t,114), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,115), markers(t,116), markers(t,117), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,79), markers(t,80), markers(t,81), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');


markers = Y_fly';
figure;
t = 282;
draw_skel(markers(t, :), lab.colheaders);
hold on;
xyz = [markers(t,109), markers(t,110), markers(t,111); markers(t,73), markers(t,74), markers(t,75); markers(t,82), markers(t,83), markers(t,84); markers(t,112), markers(t,113), markers(t,114); markers(t,115), markers(t,116), markers(t,117); markers(t,79), markers(t,80), markers(t,81)];
line(xyz(1:2, 1), xyz(1:2, 2), xyz(1:2, 3), 'Color', 'b', 'LineWidth', 6);
line(xyz(2:3, 1), xyz(2:3, 2), xyz(2:3, 3), 'Color', 'r', 'LineWidth', 6);
line(xyz(3:4, 1), xyz(3:4, 2), xyz(3:4, 3), 'Color', 'g', 'LineWidth', 6);
line(xyz(4:5, 1), xyz(4:5, 2), xyz(4:5, 3), 'Color', 'c', 'LineWidth', 6);
line(xyz([3,5], 1), xyz([3,5], 2), xyz([3,5], 3), 'Color', 'y', 'LineWidth', 6);
line(xyz(5:6, 1), xyz(5:6, 2), xyz(5:6, 3), 'Color', 'm', 'LineWidth', 6);
plot3(markers(t,73), markers(t,74), markers(t,75), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,82), markers(t,83), markers(t,84), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,112), markers(t,113), markers(t,114), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,115), markers(t,116), markers(t,117), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,79), markers(t,80), markers(t,81), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');


markers = Y_multi_bone';
figure;
t = 282;
draw_skel(markers(t, :), lab.colheaders);
hold on;
xyz = [markers(t,109), markers(t,110), markers(t,111); markers(t,73), markers(t,74), markers(t,75); markers(t,82), markers(t,83), markers(t,84); markers(t,112), markers(t,113), markers(t,114); markers(t,115), markers(t,116), markers(t,117); markers(t,79), markers(t,80), markers(t,81)];
line(xyz(1:2, 1), xyz(1:2, 2), xyz(1:2, 3), 'Color', 'b', 'LineWidth', 6);
line(xyz(2:3, 1), xyz(2:3, 2), xyz(2:3, 3), 'Color', 'r', 'LineWidth', 6);
line(xyz(3:4, 1), xyz(3:4, 2), xyz(3:4, 3), 'Color', 'g', 'LineWidth', 6);
line(xyz(4:5, 1), xyz(4:5, 2), xyz(4:5, 3), 'Color', 'c', 'LineWidth', 6);
line(xyz([3,5], 1), xyz([3,5], 2), xyz([3,5], 3), 'Color', 'y', 'LineWidth', 6);
line(xyz(5:6, 1), xyz(5:6, 2), xyz(5:6, 3), 'Color', 'm', 'LineWidth', 6);
plot3(markers(t,73), markers(t,74), markers(t,75), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,82), markers(t,83), markers(t,84), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,112), markers(t,113), markers(t,114), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,115), markers(t,116), markers(t,117), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,79), markers(t,80), markers(t,81), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');


%%
% bone
[b, bv, bls] = get_bones(Y_bone);
figure;
hold on;
plot(bls(:, 35, 31), 'b', 'DisplayName', 'RTHI-RKNE');
plot(bls(:, 33, 31), 'r', 'DisplayName', 'RKNE-RSHN');
plot(bls(:, 33, 21), 'g', 'DisplayName', 'RSHN-RANK');
ylim([0 0.5]);
xlabel('frame');
ylabel('bone length (m)'); 
legend('show');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);


% original
[b, bv, bls] = get_bones(X);
figure;
hold on;
plot(bls(:, 35, 31), 'b', 'DisplayName', 'RTHI-RKNE');
plot(bls(:, 33, 31), 'r', 'DisplayName', 'RKNE-RSHN');
plot(bls(:, 33, 21), 'g', 'DisplayName', 'RSHN-RANK');
ylim([0 0.5]);
xlabel('frame');
ylabel('bone length (m)'); 
legend('show');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);

%% plot frame
markers = X;
figure;
t = 70;
draw_skel(markers(t, :), lab.colheaders);
hold on;
xyz = [markers(t,103), markers(t,104), markers(t,105); markers(t,91), markers(t,92), markers(t,93); markers(t,97), markers(t,98), markers(t,99); markers(t,61), markers(t,62), markers(t,63)];
line(xyz(1:2, 1), xyz(1:2, 2), xyz(1:2, 3), 'Color', 'b', 'LineWidth', 6);
line(xyz(2:3, 1), xyz(2:3, 2), xyz(2:3, 3), 'Color', 'r', 'LineWidth', 6);
line(xyz(3:4, 1), xyz(3:4, 2), xyz(3:4, 3), 'Color', 'g', 'LineWidth', 6);
plot3(markers(t,91), markers(t,92), markers(t,93), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,97), markers(t,98), markers(t,99), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

markers = Y_fly;
figure;
t = 70;
draw_skel(markers(t, :), lab.colheaders);
hold on;
xyz = [markers(t,103), markers(t,104), markers(t,105); markers(t,91), markers(t,92), markers(t,93); markers(t,97), markers(t,98), markers(t,99); markers(t,61), markers(t,62), markers(t,63)];
line(xyz(1:2, 1), xyz(1:2, 2), xyz(1:2, 3), 'Color', 'b', 'LineWidth', 6);
line(xyz(2:3, 1), xyz(2:3, 2), xyz(2:3, 3), 'Color', 'r', 'LineWidth', 6);
line(xyz(3:4, 1), xyz(3:4, 2), xyz(3:4, 3), 'Color', 'g', 'LineWidth', 6);
plot3(markers(t,91), markers(t,92), markers(t,93), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,97), markers(t,98), markers(t,99), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');


markers = Y_bone;
figure;
t = 70;
draw_skel(markers(t, :), lab.colheaders);
hold on;
xyz = [markers(t,103), markers(t,104), markers(t,105); markers(t,91), markers(t,92), markers(t,93); markers(t,97), markers(t,98), markers(t,99); markers(t,61), markers(t,62), markers(t,63)];
line(xyz(1:2, 1), xyz(1:2, 2), xyz(1:2, 3), 'Color', 'b', 'LineWidth', 6);
line(xyz(2:3, 1), xyz(2:3, 2), xyz(2:3, 3), 'Color', 'r', 'LineWidth', 6);
line(xyz(3:4, 1), xyz(3:4, 2), xyz(3:4, 3), 'Color', 'g', 'LineWidth', 6);
plot3(markers(t,91), markers(t,92), markers(t,93), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot3(markers(t,97), markers(t,98), markers(t,99), 'ob', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

% line(xyz(1:2, 1), xyz(1:2, 2), xyz(1:2, 3), 'Color', 'b', 'LineWidth', 3);
% line(xyz(2:3, 1), xyz(2:3, 2), xyz(2:3, 3), 'Color', 'r', 'LineWidth', 3);
% line(xyz(3:4, 1), xyz(3:4, 2), xyz(3:4, 3), 'Color', 'g', 'LineWidth', 3);
% plot3(markers(t,91), markers(t,92), markers(t,93), 'ob', 'MarkerFaceColor', 'k');
% plot3(markers(t,97), markers(t,98), markers(t,99), 'ob', 'MarkerFaceColor', 'k');


%% plot signal
Y_spline = spline_interp(X, observed);
Y_svd = EMSVD(X, observed, 15, 1000);

figure;
subplot(5, 1, 1);
Y = X;
idx = 91:93;
plot(Y(:, idx), 'DisplayName', 'Original');
lv = [25, min(ylim);25, max(ylim)];
lv2 = [90, min(ylim);90, max(ylim)];
% line(lv(:,1), lv(:,2), 'Color', 'k');
% line(lv2(:,1), lv2(:,2), 'Color', 'k');
xlim([10, 100]);
%text(50, -1, 'abc');
ylabel('coordinate');
%xlabel('time');
title('Original');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);



subplot(5, 1, 2);
Y = Y_spline;
idx = 91:93;
plot(Y(:, idx));
lv = [25, min(ylim);25, max(ylim)];
lv2 = [90, min(ylim);90, max(ylim)];
% line(lv(:,1), lv(:,2), 'Color', 'k');
% line(lv2(:,1), lv2(:,2), 'Color', 'k');
%text(50, -1, 'abc');
xlim([10, 100]);
ylabel('coordinate');
%xlabel('time');
title('Spline');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);


subplot(5, 1, 3);
Y = Y_svd;
idx = 91:93;
plot(Y(:, idx));
lv = [25, min(ylim);25, max(ylim)];
lv2 = [90, min(ylim);90, max(ylim)];
% line(lv(:,1), lv(:,2), 'Color', 'k');
% line(lv2(:,1), lv2(:,2), 'Color', 'k');
%text(50, -1, 'abc');
xlim([10, 100]);
ylabel('coordinate');
%xlabel('time');
title('MSVD');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);


subplot(5, 1, 4);
Y = Y_fly;
idx = 91:93;
plot(Y(:, idx));
lv = [25, min(ylim);25, max(ylim)];
lv2 = [90, min(ylim);90, max(ylim)];
% line(lv(:,1), lv(:,2), 'Color', 'k');
% line(lv2(:,1), lv2(:,2), 'Color', 'k');
%text(50, -1, 'abc');
xlim([10, 100]);
ylabel('coordinate');
%xlabel('time');
title('LDS (baseline)');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);

subplot(5, 1, 5);
Y = Y_bone;
idx = 91:93;
plot(Y(:, idx));
lv = [25, min(ylim);25, max(ylim)];
lv2 = [90, min(ylim);90, max(ylim)];
% line(lv(:,1), lv(:,2), 'Color', 'k');
% line(lv2(:,1), lv2(:,2), 'Color', 'k');
%text(50, -1, 'abc');
xlim([10, 100]);
ylabel('coordinate');
xlabel('time');
title('BoLeRO');
[lx, ly] = dsxy2figxy(gca, lv(:, 1), lv(:, 2));
annotation('line', lx, ly);
[lx, ly] = dsxy2figxy(gca, lv2(:, 1), lv2(:, 2));
annotation('line', lx, ly);

%% plot data matrix

figure;
oo = true(size(observed));
oo(1,1) = false;
imagesc(oo');
colormap(gray);

figure;
imagesc(observed');
colormap(gray);

figure;
oo = ~observed;
imagesc(oo');
colormap(gray);


