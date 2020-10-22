function [U] = solveD2d(V, N)
% solve 2-d anderson source problem
% - u''(x) + V(x) u(x) = 1 for x in [0,1] x [0,1]
% u(x) = 0 for x on boundary
% V(x) is piecewise constant
% input:
%     V(2-d array):   piecewise constant of V(x)
%     N(integer):     degree of polynomials (default N = 10)
% output:
%     U(1-d array):   array with size (1, (M*N+1)^2)
%                     each column represents projection on polynomial basis
%                     use function getval2d to get value of sulution

% V must be square !!!

% default input
if nargin < 2
    N = 6;
end

M = size(V, 1);
hm = 1 / M;

[Ahat, Bhat, Fhat] = lgmat2d(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);
[iFhat, ~, vFhat] = find(Fhat);

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
F = zeros((M*N+1)^2, 1);

kA = 0;
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

        ind = l2g(m1, m2, iFhat);
        F(ind) = F(ind) + hm*hm/4 * vFhat;
    end
end

A = sparse(iA, jA, vA, (M*N+1)^2, (M*N+1)^2);

tmp = reshape(1:(M*N+1)^2, (M*N+1), (M*N+1));
edge = [tmp(1,:), tmp(end,:), tmp(:,1)', tmp(:,end)'];
A(edge,:) = []; A(:,edge) = [];
F(edge) = [];

UU = A \ F;

U = zeros((M*N+1)^2, 1);
ind = setdiff(tmp(:), edge);
U(ind) = UU;
end