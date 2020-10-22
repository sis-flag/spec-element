function u = getval2d(U, Ns, N)
% get value from projection on polynomial basis (2-d)
% input:
%     U(1-d array): projection on polynomial basis
%     Ns(integer):  number of sample points in each interval (default Ns = 20)
%     N(integer):   degree of polynomials (default N = 10)
% output:
%     u(2-d array): function values on sample points

% default input
if nargin < 3
    N = 6;
end
if nargin < 2
    Ns = 20;
end

M = round(sqrt(length(U))-1) / N;

U = reshape(U, M*N+1, M*N+1);

    function y = phi(n, x)
        if n == 1
            y = (1-x)/2;
        elseif n == N+1
            y = (1+x)/2;
        else
            p1 = legendre(n,x); p1 = p1(1,:,:);
            p2 = legendre(n-2,x); p2 = p2(1,:,:);
            y = (p1 - p2) / sqrt(4*n-2);
        end
    end

xhat = linspace(-1, 1, Ns +1);
yhat = zeros(N+1, Ns+1);
for n = 1:N+1
    yhat(n,:) = phi(n, xhat);
end

u = zeros(1, Ns*M+1);
for m1 = 1:M
    for m2 = 1:M
        Uloc = U((m1-1)*N+1: m1*N+1, (m2-1)*N+1: m2*N+1);
        u((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1) =...
            (yhat') * Uloc * yhat;
    end
end

end
