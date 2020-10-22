function [Ahat, Bhat, Fhat, Hhat0, Hhat1, Ghat0, Ghat1] = lgmat(N)
% base matrix of legrend polynomials

iA = [0, N, 0, N, 1:N-1] +1;
jA = [0, N, N, 0, 1:N-1] +1;
vA = [1/2, 1/2, -1/2, -1/2, ones(1,N-1)];
Ahat = sparse(iA, jA, vA, N+1, N+1);

iB = [0, N, 0, 1, 0, 2, 0, N, 1, N, 2, N] +1;
jB = [0, N, 1, 0, 2, 0, N, 0, N, 1, N, 2] +1;
vB = [2/3, 2/3, -1/sqrt(6), -1/sqrt(6), 1/sqrt(90), 1/sqrt(90), ...
    1/3, 1/3, -1/sqrt(6), -1/sqrt(6), -1/sqrt(90), -1/sqrt(90) ];
k1 = 1:N-1;  k2 = 1:N-3;
iB = [iB, k1+1, k2+3, k2+1];
jB = [jB, k1+1, k2+1, k2+3];
vB = [vB, 2./(2*k1-1)./(2*k1+3),...
    -1./(2*k2+3)./sqrt(2*k2+1)./sqrt(2*k2+5),...
    -1./(2*k2+3)./sqrt(2*k2+1)./sqrt(2*k2+5)];
Bhat = sparse(iB, jB, vB, N+1, N+1);

Fhat = sparse([1,N+1,2],[1,1,1],[1,1,-2/sqrt(6)], N+1, 1);

Hhat0 = sparse(1,1,1, N+1, N+1);
Hhat1 = sparse(N+1,N+1, 1, N+1, N+1);

Ghat0 = sparse(1,1,1, N+1, 1);
Ghat1 = sparse(N+1,1,1, N+1, 1);
