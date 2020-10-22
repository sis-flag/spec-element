function [x, w] = lgpoints(N)
% return Legendre-Gauss points on [-1, 1]

% Jacobi matrix
j = 1:N;
A = diag(j./sqrt((2*j-1).*(2*j+1)),1)...
  + diag(j./sqrt((2*j-1).*(2*j+1)),-1);
x = sort(eig(sparse(A)));

if nargout < 2
    return
end

% legendre polynomial
Lx = zeros(length(x), N);
Lx(:,1) = x;
Lx(:,2) = 1.5 * x.*x - 0.5;
for k = 3:N
    % k L_k(x) = (2k-1) x L_{k-1}(x) - (k-1) L_{k-2}(x)
    Lx(:,k) = ((2*k-1)*x.*Lx(:,k-1) - (k-1)*Lx(:,k-2)) / k;
end

% derivative of legendre polynomial
dLx = zeros(length(x), N+1);
dLx(:,1) = 1;
dLx(:,2) = 3 * x;
for k = 3:N+1
    % L_k'(x) = L_{k-2}'(x) + (2k-1) L_{k-1}(x)
    dLx(:,k) = dLx(:,k-2) + (2*k-1) * Lx(:,k-1);
end

w = 2./ ((1-x.^2).*dLx(:,N+1).^2);