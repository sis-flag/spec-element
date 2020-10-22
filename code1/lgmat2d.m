function [KA, KB, KF, KHB0, KHB1, KBH0, KBH1, KGF0, KGF1, KFG0, KFG1] = lgmat2d(N)
% base matrix of legrend polynomials (2-d)

[Ahat, Bhat, Fhat, Hhat0, Hhat1, Ghat0, Ghat1] = lgmat(N);

KA = kron(Ahat, Bhat) + kron(Bhat, Ahat);
KB = kron(Bhat, Bhat);
KF = kron(Fhat, Fhat);

KHB0 = kron(Hhat0, Bhat);
KHB1 = kron(Hhat1, Bhat);
KBH0 = kron(Bhat, Hhat0);
KBH1 = kron(Bhat, Hhat1);

KGF0 = kron(Ghat0, Fhat);
KGF1 = kron(Ghat1, Fhat);
KFG0 = kron(Fhat, Ghat0);
KFG1 = kron(Fhat, Ghat1);
