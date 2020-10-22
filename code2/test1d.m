clear

% prob = test_prob1d('D', 0.5);
prob = test_prob1d('R', 1);
mesh = mesh1d([0, 0.1, 0.2, 0.5, 0.8, 1.0]);

% test source problem
U = solve1d(prob, mesh);
[x, u, du] = getval1d(U, mesh);

figure
plot(x, u, x, prob.u(x))
disp(max(abs(u - prob.u(x))))

figure
plot(x, du, x, prob.du(x))
disp(max(abs(du - prob.du(x))))

% test eigen problem
[U, lam] = eig1d(prob, mesh);

figure
hold on
for k = 1:6
    [x, u, du] = getval1d(U(:,k), mesh);
    
%     h = diff(x);
%     dpu = (u(3:end) - u(1:end-2))  ./ (h(1:end-1) + h(2:end));
%     adpu = prob.a(x(2:end-1)) .* dpu;
%     dadpu = (adpu(3:end) - adpu(1:end-2)) ./ (h(2:end-2) + h(3:end-1));
%     bdpu = prob.b(x(3:end-2)) .* dpu(2:end-1);
%     cppu = prob.c(x(3:end-2)) .* u(3:end-2);
    
    plot(x, u)
%     plot(x(3:end-2), (-dadpu + bdpu + cppu)/lam(k));
end