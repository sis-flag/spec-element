function [U, lam] = eigP1d(V, num, N)
% solve 1-d anderson eigen problem
% - u''(x) + V(x) u(x) = lam u(x) for x in [0, 1]
% parabolic boundary condation
% V(x) is piecewise constant
% input:
%     V(1-d array):   piecewise constant of V(x)
%     num(integer):   number of eigenvalues required
%     N(integer):     degree of polynomials (default N = 10)
% output:
%     U(2-d array):   array with size (num, M*N+1)
%                     each column represents projection on polynomial basis
%                     use function getval1d to get value of sulution
%     lam(1-d array): array with size (num, 1), eigenvalues

% default input
if nargin < 3
    N = 6;
end

M = length(V);
hm = 1 / M;

[Ahat, Bhat] = lgmat(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);

    function ind = l2g(m, n)
        ind = (m-1)*N + n;
        ind(ind==M*N+1) = 1;
    end

nnzA = length(iAhat); nnzB = length(iBhat);
iA = zeros(1, M*(nnzA+nnzB));
jA = zeros(1, M*(nnzA+nnzB));
vA = zeros(1, M*(nnzA+nnzB));
iB = zeros(1, M*(nnzB));
jB = zeros(1, M*(nnzB));
vB = zeros(1, M*(nnzB));

kA = 0; kB = 0;
for m =1:M
    iA(kA+1:kA+nnzA) = l2g(m, iAhat);
    jA(kA+1:kA+nnzA) = l2g(m, jAhat);
    vA(kA+1:kA+nnzA) = 2/hm * vAhat;
    kA = kA+nnzA;
    
    iA(kA+1:kA+nnzB) = l2g(m, iBhat);
    jA(kA+1:kA+nnzB) = l2g(m, jBhat);
    vA(kA+1:kA+nnzB) = hm/2 * V(m) * vBhat;
    kA = kA+nnzB;
    
    iB(kB+1:kB+nnzB) = l2g(m, iBhat);
    jB(kB+1:kB+nnzB) = l2g(m, jBhat);
    vB(kB+1:kB+nnzB) = hm/2 * vBhat;
    kB = kB+nnzB;
end

A = sparse(iA, jA, vA, M*N, M*N);
B = sparse(iB, jB, vB, M*N, M*N);

[UU, lam] = eigs(A, B, num, 0);
lam = diag(lam);

U = [UU; UU(1,:)];
end