function [x, y, u, dxu, dyu] = getval2d(U, mesh, Ns, N)
% get value from projection on polynomial basis (2-d)
% input:
%     U(1-d array): projection on polynomial basis
%     Ns(integer):  number of sample points in each interval (default Ns = 20)
%     N(integer):   degree of polynomials (default N = 10)
% output:
%     x,y(2-d array):      sample points
%     u(2-d array):        function values on sample points
%     dxu,dyu(2-d array):  derivative of function on sample points

% default input
if nargin < 4
    N = 6;
end
if nargin < 3
    Ns = 10;
end
Mx = mesh.Mx; My = mesh.My;
hx = mesh.hx; hy = mesh.hy;

U = reshape(U, Mx*N+1, My*N+1);

% basis on refenence interval
xhat = linspace(-1, 1, Ns +1);
[uhat, duhat] = basis(N, xhat);

x = zeros(1, Ns*Mx+1);
for m = 1:Mx
    x((m-1)*Ns+1: m*Ns+1) = (xhat + 1)/2 * hx(m) + mesh.x(m);
end
y = zeros(1, Ns*My+1);
for m = 1:My
    y((m-1)*Ns+1: m*Ns+1) = (xhat + 1)/2 * hy(m) + mesh.y(m);
end
[y, x] = meshgrid(y, x); % caution!

u = zeros(Ns*Mx+1, Ns*My+1);
dxu = zeros(Ns*Mx+1, Ns*My+1);
dyu = zeros(Ns*Mx+1, Ns*My+1);
for mx = 1:Mx
    for my = 1:My
        Uloc = U((mx-1)*N+1: mx*N+1, (my-1)*N+1: my*N+1);
        u((mx-1)*Ns+1: mx*Ns+1, (my-1)*Ns+1: my*Ns+1) =...
            uhat * Uloc * uhat';
        dxu((mx-1)*Ns+1: mx*Ns+1, (my-1)*Ns+1: my*Ns+1) =...
            duhat * Uloc * uhat' / (hx(mx)/2);
        dyu((mx-1)*Ns+1: mx*Ns+1, (my-1)*Ns+1: my*Ns+1) =...
            uhat * Uloc * duhat' / (hy(my)/2);
    end
end

end
