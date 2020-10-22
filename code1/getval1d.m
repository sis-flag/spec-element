function u = getval1d(U, Ns, N)
% get value from projection on polynomial basis (1-d)
% input:
%     U(1-d array): projection on polynomial basis
%     Ns(integer):  number of sample points in each interval (default Ns = 50)
%     N(integer):   degree of polynomials (default N = 10)
% output:
%     u(1-d array): function values on sample points

% default input
if nargin < 3
    N = 6;
end
if nargin < 2
    Ns = 50;
end

M = round(length(U)-1) / N;

U = reshape(U, M*N+1, 1);

    function y = phi(n, x)
        if n == 1
            y = (1-x)/2;
        elseif n == N+1
            y = (1+x)/2;
        else
            p1 = legendre(n,x); p1 = p1(1,:);
            p2 = legendre(n-2,x); p2 = p2(1,:);
            y = (p1 - p2) / sqrt(4*n-2);
        end
    end

xhat = linspace(-1, 1, Ns +1);
yhat = zeros(N+1, Ns+1);
for n = 1:N+1
    yhat(n,:) = phi(n, xhat);
end

u = zeros(1, Ns*M+1);
for m = 1:M
    Uloc = U((m-1)*N+1: m*N+1);
    u((m-1)*Ns+1: m*Ns+1) = (Uloc') * yhat;
end
end
