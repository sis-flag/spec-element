function U = solve1d(prob, mesh, N)

% defult input
if nargin < 3
    N = 6;
end
M = mesh.M;
h = mesh.h;

% Gauss nodes
[xn, wn] = lgpoints(N);

% basis on Gauss nodes
[phi, dphi] = basis(N, xn);

% local to global
if prob.bd == 'P'
    l2g = @(m, n) mod((m-1) * N + (n-1), M*N) + 1;
else
    l2g = @(m, n) (m-1) * N + n;
end

A = zeros(M*N+1, M*N+1);
F = zeros(M*N+1, 1);

for m =1:M
    lxn = (xn + 1)/2 * h(m) + mesh.x(m);
    lan = arrayfun(prob.a, lxn);
    lbn = arrayfun(prob.b, lxn);
    lcn = arrayfun(prob.c, lxn);
    lfn = arrayfun(prob.f, lxn);
    
    lA = (2/h(m)) * dphi' * diag(wn.*lan) * dphi ...
       + phi' * diag(wn.*lbn) * dphi ...
       + (h(m)/2) * phi' * diag(wn.*lcn) * phi;
       
    lF = (h(m)/2) * phi' * (lfn.*wn);
        
%     [liA, ljA, lvA] = find(lA);
%     iA(nnzA+1: nnzA+(N+1)^2) = l2g(m, liA);
%     jA(nnzA+1: nnzA+(N+1)^2) = l2g(m, ljA);
%     vA(nnzA+1: nnzA+(N+1)^2) = lvA;
%     nnzA = nnzA + (N+1)^2;
    
    gi = l2g(m, 1:N+1);
    A(gi, gi) = A(gi, gi) + lA;
    F(gi) = F(gi) + lF;
end


if prob.bd == 'D'
    A(1,:) = 0; A(end,:) = 0;
    A(1,1) = 1; A(end,end) = 1;
    F(1) = prob.u0(0);
    F(end) = prob.u0(1);
    
    U = A \ F;
    
elseif prob.bd == 'R'
    A(1,1) = A(1,1) + prob.h(0);
    A(end,end) = A(end,end) + prob.h(1);
    F(1,1) = F(1,1) + prob.g(0);
    F(end,end) = F(end,end) + prob.g(1);
    
    U = A \ F;
    
elseif prob.bd == 'P'
    
    U = A \ F;
    U = [U; U(1)];
end
