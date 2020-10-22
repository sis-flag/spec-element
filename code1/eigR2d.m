function [U, lam] = eigR2d(V, h, num, N)
% solve 2-d anderson eigen problem
% - u''(x) + V(x) u(x) = lam u(x) for x in [0,1] x [0,1]
% u'(x) + h u(x) = 0 for x on boundary
% V(x) is piecewise constant
% input:
%     V(2-d array):   piecewise constant of V(x)
%     h(number):      paremeter in equation
%     num(integer):   number of eigenvalues required
%     N(integer):     degree of polynomials (default N = 10)
% output:
%     U(2-d array):   array with size (num, (M*N+1)^2)
%                     each column represents projection on polynomial basis
%                     use function getval2d to get value of sulution
%     lam(1-d array): array with size (num, 1), eigenvalues

% V must be square !!!

% default input
if nargin < 4
    N = 6;
end

M = size(V, 1);
hm = 1 / M;

[Ahat, Bhat, ~, KHB0, KHB1, KBH0, KBH1, ~, ~, ~, ~] = lgmat2d(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);

    function u = l2g(m1, m2, n)
        n1 = mod(n-1, N+1) +1;
        n2 = ceil(n / (N+1));
        u1 = (m1-1)*N + n1;
        u2 = (m2-1)*N + n2;
        u = (u2-1)*(M*N+1) + u1;
    end

nnzA = length(iAhat); nnzB = length(iBhat);
iA = zeros(1, M*M*(nnzA+nnzB));
jA = zeros(1, M*M*(nnzA+nnzB));
vA = zeros(1, M*M*(nnzA+nnzB));
iB = zeros(1, M*M*(nnzB));
jB = zeros(1, M*M*(nnzB));
vB = zeros(1, M*M*(nnzB));

kA = 0; kB = 0;
for m1 =1:M
    for m2 = 1:M
        iA(kA+1:kA+nnzA) = l2g(m1, m2, iAhat);
        jA(kA+1:kA+nnzA) = l2g(m1, m2, jAhat);
        vA(kA+1:kA+nnzA) = vAhat;
        kA = kA+nnzA;
        
        iA(kA+1:kA+nnzB) = l2g(m1, m2, iBhat);
        jA(kA+1:kA+nnzB) = l2g(m1, m2, jBhat);
        vA(kA+1:kA+nnzB) = hm*hm/4 * V(m1,m2) * vBhat;
        kA = kA+nnzB;
        
        iB(kB+1:kB+nnzB) = l2g(m1, m2, iBhat);
        jB(kB+1:kB+nnzB) = l2g(m1, m2, jBhat);
        vB(kB+1:kB+nnzB) = hm*hm/4 * vBhat;
        kB = kB+nnzB;
    end
end

A = sparse(iA, jA, vA, (M*N+1)^2, (M*N+1)^2);
B = sparse(iB, jB, vB, (M*N+1)^2, (M*N+1)^2);

% left boundary
m2 = 1;
for m1 = 1:M
    ind = l2g(m1, m2, 1:(N+1)^2);
    A(ind, ind) = A(ind, ind) + hm/2 * h * KHB0;
end

% right boundary
m2 = M;
for m1 = 1:M
    ind = l2g(m1, m2, 1:(N+1)^2);
    A(ind, ind) = A(ind, ind) + hm/2 * h * KHB1;
end

% bottom boundary
m1 = 1;
for m2 = 1:M
    ind = l2g(m1, m2, 1:(N+1)^2);
    A(ind, ind) = A(ind, ind) + hm/2 * h * KBH0;
end

% top boundary
m1 = M;
for m2 = 1:M
    ind = l2g(m1, m2, 1:(N+1)^2);
    A(ind, ind) = A(ind, ind) + hm/2 * h * KBH1;
end

[U, lam] = eigs(A, B, num, 0);
lam = diag(lam);
end