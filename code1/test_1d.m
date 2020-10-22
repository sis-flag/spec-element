clear;
rng(0);

% peremeters
% K = 8000;
% V = rand(1, 20);
% V = K * V;

V = 10*[9, 4, 8, 2, 7, 6, 1, 3, 4, 5];
h = 0;
beta = 100;

% plot potential
figure();
cf = bar(((1:length(V))-0.5)/length(V), V);
cf.BarWidth = 1;

% Robin boundary
[U, lam] = eigR1d(V, h, 6);
W = solveN1d(V, h/beta);

w = getval1d(W);
x = linspace(0, 1, length(w));

figure();
hold on
plot(x, w, 'k')
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / (lam(k) + 1/max(w) + beta) )
end

% Dirichlet boundary
[U, lam] = eigD1d(V, 6);
W = solveD1d(V);

w = getval1d(W);
x = linspace(0, 1, length(w));

figure();
hold on
% plot(x, w, 'k')
for k = 1:1
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / (lam(k) + 1/max(w) + beta) )
end