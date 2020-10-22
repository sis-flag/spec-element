clear;
rng(0);

% peremeters
K = 8000;
V = rand(20);
V = K * V;

h = 0;
beta = 100;

% plot potential
figure();
% bar3(V);
M = size(V, 1);
pV = [V, zeros(M,1); zeros(1,M+1)];
px = linspace(0, 1, M+1);
[px2, px1] = meshgrid(px, px);
cf = pcolor(px1, px2, pV);
cf.LineStyle = 'none';
caxis([0, K])
axis square
colorbar;

% Robin boundary
[U, lam] = eigR2d(V, h, 6);
W = solveN2d(V, h/beta);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

% plot landscape
figure();
mesh(x1, x2, w);

% plot valley line
figure();
hold on
cf = pcolor(x1, x2, w);
cf.LineStyle = 'none';
axis square
colorbar;
[vx1, vx2] = valley_line(w);
cf = plot(vx1, vx2, 'k.');
cf.MarkerSize = 3;

% plot eigen mode
figure();
for k = 1:4
    subplot(2, 2, k);
    uk = my_nmlz(getval2d(U(:,k)));
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
end

% Dirichlet boundary
[U, lam] = eigD2d(V, 6);
W = solveD2d(V);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

% plot landscape
figure();
mesh(x1, x2, w);

% plot valley line
figure();
hold on
cf = pcolor(x1, x2, w);
cf.LineStyle = 'none';
axis square
colorbar;
[vx1, vx2] = valley_line(w);
cf = plot(vx1, vx2, 'k.');
cf.MarkerSize = 3;
title('landscape')

% plot eigen mode
figure();
for k = 1:4
    subplot(2, 2, k);
    uk = my_nmlz(getval2d(U(:,k)));
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
    title(sprintf('eigenmode %d', k));
end