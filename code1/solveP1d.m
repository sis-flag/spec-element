function [U] = solveP1d(V, N)
% solve 1-d anderson source problem
% - u''(x) + V(x) u(x) = 1 for x in [0, 1]
% parabolic boundary condation
% V(x) is piecewise constant
% input:
%     V(1-d array):   piecewise constant of V(x)
%     N(integer):     degree of polynomials (default N = 10)
% output:
%     U(1-d array):   array with size (1, M*N+1)
%                     each column represents projection on polynomial basis
%                     use function getval1d to get value of sulution

% default input
if nargin < 2
    N = 6;
end

M = length(V);
hm = 1 / M;

[Ahat, Bhat, Fhat] = lgmat(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);
[iFhat, ~, vFhat] = find(Fhat);

    function ind = l2g(m, n)
        ind = (m-1)*N + n;
        ind(ind==M*N+1) = 1;
    end

nnzA = length(iAhat); nnzB = length(iBhat);
iA = zeros(1, M*(nnzA+nnzB));
jA = zeros(1, M*(nnzA+nnzB));
vA = zeros(1, M*(nnzA+nnzB));
F = zeros(M*N, 1);

kA = 0;
for m =1:M
    iA(kA+1:kA+nnzA) = l2g(m, iAhat);
    jA(kA+1:kA+nnzA) = l2g(m, jAhat);
    vA(kA+1:kA+nnzA) = 2/hm * vAhat;
    kA = kA+nnzA;
    
    iA(kA+1:kA+nnzB) = l2g(m, iBhat);
    jA(kA+1:kA+nnzB) = l2g(m, jBhat);
    vA(kA+1:kA+nnzB) = hm/2 * V(m) * vBhat;
    kA = kA+nnzB;

    ind = l2g(m, iFhat);
    F(ind) = F(ind) + hm/2 * vFhat;
end

A = sparse(iA, jA, vA, M*N, M*N);

UU = A \ F;

U = [UU; UU(1)];
end