clear

prob = test_prob2d_1(); % problem
% mesh = mesh2d([0, 0.1, 0.2, 0.5, 0.8, 1.0], ...
%             [0, 0.1, 0.5, 0.7, 1.0]); % ununiform mesh
mesh = mesh2d(10);

% test source problem
U = solve2d(prob, mesh);

[px, py, pu, pux, puy] = getval2d(U, mesh);

figure
surf(px, py, pu);
figure
surf(px, py, prob.u(px, py));
figure
surf(px, py, pu-prob.u(px, py));

disp(max(max(abs(pu-prob.u(px, py)))))
disp(max(max(abs(pux-prob.dux(px, py)))))
disp(max(max(abs(puy-prob.duy(px, py)))))

% test eigen problem
[U, lam] = eig2d(prob, mesh);

disp(lam)

for k = 1:6
    [px, py, pu, pux, puy] = getval2d(U(:,k), mesh, 20);
    
    figure
    surf(px, py, abs(pu));
    
    ph = px(2,1) - px(1,1);
    lap_pu = conv2(pu, [0,-1,0;-1,4,-1;0,-1,0], 'valid') / (ph)^2;
    pb1 = arrayfun(prob.b1, px, py);
    pb2 = arrayfun(prob.b2, px, py);
    pc = arrayfun(prob.c, px, py);
    err = lap_pu...
        + pb1(2:end-1,2:end-1).*pux(2:end-1,2:end-1)...
        + pb2(2:end-1,2:end-1).*puy(2:end-1,2:end-1)...
        + pc(2:end-1,2:end-1) .*pu(2:end-1,2:end-1)...
        - lam(k) * pu(2:end-1,2:end-1);
    
    figure
    surf(px(2:end-1,2:end-1), py(2:end-1,2:end-1), abs(err));
    
    disp(max(max(abs(err))))
end