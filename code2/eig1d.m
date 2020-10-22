function [U, lam] = eig1d(prob, mesh, num, N)

% defult input
if nargin < 4
    N = 6;
end
if nargin < 3
    num = 6;
end
M = mesh.M;
h = mesh.h;

% Gauss nodes
[xn, wn] = lgpoints(N);

% basis on Gauss nodes
[phi, dphi] = basis(N, xn);
lP0 = phi' * diag(wn) * phi;

% local to global
if prob.bd == 'P'
    l2g = @(m, n) mod((m-1) * N + (n-1), M*N) + 1;
else
    l2g = @(m, n) (m-1) * N + n;
end

A = zeros(M*N+1, M*N+1);
B = zeros(M*N+1, M*N+1);

for m =1:M
    lxn = (xn + 1)/2 * h(m) + mesh.x(m);
    lan = arrayfun(prob.a, lxn);
    lbn = arrayfun(prob.b, lxn);
    lcn = arrayfun(prob.c, lxn);
    
    lA = (2/h(m)) * dphi' * diag(wn.*lan) * dphi ...
       + phi' * diag(wn.*lbn) * dphi ...
       + (h(m)/2) * phi' * diag(wn.*lcn) * phi;
        
    gi = l2g(m, 1:N+1);
    A(gi, gi) = A(gi, gi) + lA;
    B(gi, gi) = B(gi, gi) + (h(m)/2) * lP0;
end


if prob.bd == 'D'
    edge = [1, M*N+1];
    A(edge,:) = []; A(:,edge) = [];
    B(edge,:) = []; B(:,edge) = [];
    
    [U, lam] = eigs(A, B, num, 0);
    [~, ind] = sort(abs(diag(lam)));
    lam = diag(lam(ind, ind));
    U = [zeros(1,num); U(:,ind); zeros(1,num)];
    
elseif prob.bd == 'R'
    A(1,1) = A(1,1) + prob.h(0);
    A(end,end) = A(end,end) + prob.h(1);
    
    [U, lam] = eigs(A, B, num, 0);
    [~, ind] = sort(abs(diag(lam)));
    lam = diag(lam(ind, ind));
    U = U(:,ind);
    
elseif prob.bd == 'P'
    
    [U, lam] = eigs(A, B, num, 0);
    [~, ind] = sort(abs(diag(lam)));
    lam = diag(lam(ind, ind));
    U = [U(:,ind); U(1,ind)];
end
