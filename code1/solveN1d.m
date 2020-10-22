function [U] = solveN1d(V, g, x0, u0, N)
% solve 1-d anderson source problem
% - u''(x) + V(x) u(x) = 1 for x in [0, 1]
% u'(x) = g for x = 0 or x = 1
% enforce u(x0) = u0
% V(x) is piecewise constant
% input:
%     V(1-d array):   piecewise constant of V(x)
%     g(number):      paremeter of boundary condition
%     x0, u0(number): enforce condition (0 < x0 < 1) (default no-enforce)
%     N(integer):     degree of polynomials (default N = 10)
% output:
%     U(1-d array):   array with size (1, M*N+1)
%                     each column represents projection on polynomial basis
%                     use function getval1d to get value of sulution

% default input
if nargin < 5
    N = 6;
end
if nargin < 4
    x0 = []; u0 = [];
end

M = length(V);
hm = 1 / M;

[Ahat, Bhat, Fhat] = lgmat(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);
[iFhat, ~, vFhat] = find(Fhat);

l2g = @(m, n) (m-1) * N + n;

nnzA = length(iAhat); nnzB = length(iBhat);
iA = zeros(1, M*(nnzA+nnzB));
jA = zeros(1, M*(nnzA+nnzB));
vA = zeros(1, M*(nnzA+nnzB));
F = zeros(M*N+1, 1);

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

A = sparse(iA, jA, vA, M*N+1, M*N+1);

F(1) = F(1) + g;
F(end) = F(end) + g;

if ~isempty(x0)
    m = round(x0/hm);
    mm = l2g(m, 1);
    A(mm, :) = 0;
    A(mm, mm) = 1;
    F(mm) = u0;
end

U = A \ F;
end