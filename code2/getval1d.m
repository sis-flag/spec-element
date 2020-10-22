function [x, u, du] = getval1d(U, mesh, Ns, N)
% get value from projection on polynomial basis (1-d)
% input:
%     U(1-d array): projection on polynomial basis
%     Ns(integer):  number of sample points in each interval
%     N(integer):   degree of polynomials
% output:
%     x(1-d array):  sample points
%     u(1-d array):  function values on sample points
%     du(1-d array): derivative of function on sample points

% default input
if nargin < 4
    N = 6;
end
if nargin < 3
    Ns = 20;
end
M = mesh.M;
h = mesh.h;

U = reshape(U, M*N+1, 1);

% basis on refenence interval
xhat = linspace(-1, 1, Ns +1);
[uhat, duhat] = basis(N, xhat);

% scale to current interval
u = zeros(1, Ns*M+1);
du = zeros(1, Ns*M+1);
x = zeros(1, Ns*M+1);
for m = 1:M
    x((m-1)*Ns+1: m*Ns+1) = (xhat + 1)/2 * h(m) + mesh.x(m);
    
    Uloc = U((m-1)*N+1: m*N+1);
    u((m-1)*Ns+1: m*Ns+1) = uhat * Uloc;
    du((m-1)*Ns+1: m*Ns+1) = duhat * Uloc / (h(m)/2);
end

end
