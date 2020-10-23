function [U, lam] = eig2d(prob, mesh, num, N)

% defult input
if nargin < 4
    N = 6;
end
if nargin < 3
    num = 6;
end
Mx = mesh.Mx; My = mesh.My;
hx = mesh.hx; hy = mesh.hy;

% Gauss nodes
[xn, wn] = lgpoints(N);
[yn2d, xn2d] = meshgrid(xn, xn); % caution!
xn2d = reshape(xn2d, [], 1);
yn2d = reshape(yn2d, [], 1);
wn2d = kron(wn, wn);

% basis on Gauss nodes
[phi, dphi] = basis(N, xn);
phi0 = kron(phi, phi);
phix = kron(phi, dphi);
phiy = kron(dphi, phi);

phi2d = phi0' * diag(wn2d) * phi0;

% local to global

    function gn = l2g(mx, my, nx, ny)
        
        if nargin < 4
            n = nx;
            nx = mod(n-1, N+1) + 1;
            ny = floor((n-1)/(N+1)) + 1;
        end
        
        gnx = (mx-1) * N + nx;
        gny = (my-1) * N + ny;
        if prob.bd == 'P'
            if gnx == Mx*N + 1
                gnx = 1;
            end
            if gny == My*N + 1
                gny = 1;
            end
            gn = gnx + (gny-1)* (Mx*N);
        else
            gn = gnx + (gny-1)* (Mx*N+1);
        end
    end

nnzA = 0;
iA = zeros(Mx*My*(N+1)^4, 1);
jA = zeros(Mx*My*(N+1)^4, 1);
vA = zeros(Mx*My*(N+1)^4, 1);
nnzB = 0;
iB = zeros(Mx*My*(N+1)^4, 1);
jB = zeros(Mx*My*(N+1)^4, 1);
vB = zeros(Mx*My*(N+1)^4, 1);
for mx =1:Mx
    for my = 1:My
        lxn = (xn2d + 1)/2 * hx(mx) + mesh.x(mx);
        lyn = (yn2d + 1)/2 * hy(my) + mesh.y(my);
        la11n = arrayfun(prob.a11, lxn, lyn);
        la12n = arrayfun(prob.a12, lxn, lyn);
        la21n = arrayfun(prob.a21, lxn, lyn);
        la22n = arrayfun(prob.a22, lxn, lyn);
        lb1n = arrayfun(prob.b1, lxn, lyn);
        lb2n = arrayfun(prob.b2, lxn, lyn);
        lcn = arrayfun(prob.c, lxn, lyn);
        
        lA = (hy(my)/hx(mx)) * phix' * diag(wn2d.*la11n) * phix ...
            + phix' * diag(wn2d.*la12n) * phiy ...
            + phiy' * diag(wn2d.*la21n) * phix ...
            + (hx(mx)/hy(my)) * phiy' * diag(wn2d.*la22n) * phiy ...
            + (hy(my)/2) * phi0' * diag(wn2d.*lb1n) * phix ...
            + (hx(mx)/2) * phi0' * diag(wn2d.*lb2n) * phiy ...
            + (hx(mx)*hy(my)/4) * phi0' * diag(wn2d.*lcn) * phi0;
        
        lB = (hx(mx)*hy(my)/4) * phi2d;
        
%         gind = l2g(mx, my, 1:(N+1)^2);
%         A(gind, gind) = A(gind, gind) + lA;
%         B(gind, gind) = B(gind, gind) + lB;

        gind = l2g(mx, my, 1:(N+1)^2);
        [gj, gi] = meshgrid(gind, gind);
        
        iA(nnzA+1: nnzA+(N+1)^4) = gi(:);
        jA(nnzA+1: nnzA+(N+1)^4) = gj(:);
        vA(nnzA+1: nnzA+(N+1)^4) = lA(:);
        nnzA = nnzA + (N+1)^4;
        
        iB(nnzB+1: nnzB+(N+1)^4) = gi(:);
        jB(nnzB+1: nnzB+(N+1)^4) = gj(:);
        vB(nnzB+1: nnzB+(N+1)^4) = lB(:);
        nnzB = nnzB + (N+1)^4;
    end
end

A = sparse(iA, jA, vA, (Mx*N+1)*(My*N+1), (Mx*N+1)*(My*N+1));
B = sparse(iB, jB, vB, (Mx*N+1)*(My*N+1), (Mx*N+1)*(My*N+1));

if prob.bd == 'D'
    temp = reshape(1:(Mx*N+1)*(My*N+1), Mx*N+1, My*N+1);
    edge = [temp(:,1); temp(:,end); temp(1,:)'; temp(end,:)'];
    A(edge,:) = []; A(:,edge) = [];
    B(edge,:) = []; B(:,edge) = [];
    
    [UU, lam] = eigs(A, B, num, 0);
    [~, ind] = sort(abs(diag(lam)));
    lam = diag(lam(ind, ind));
    U = zeros((Mx*N+1)*(My*N+1), num);
    o_o = setdiff(1:(Mx*N+1)*(My*N+1), edge);
    U(o_o,:) = UU(:,ind);
    
elseif prob.bd == 'R'
    
    % bottom boundary
    my = 1;
    for mx = 1:Mx
        lxn = (xn + 1)/2 * hx(mx) + mesh.x(mx);
        lyn = ones(N+1, 1) * mesh.y(my);
        lhn = arrayfun(prob.h, lxn, lyn);
        
        lA = (hx(mx)/2) * phi' * diag(wn.*lhn) * phi;
        
        gind = l2g(mx, my, 1:N+1, 1);
        A(gind, gind) = A(gind, gind) + lA;
    end
    
    % top boundary
    my = My;
    for mx = 1:Mx
        lxn = (xn + 1)/2 * hx(mx) + mesh.x(mx);
        lyn = ones(N+1, 1) * mesh.y(my+1);
        lhn = arrayfun(prob.h, lxn, lyn);
        
        lA = (hx(mx)/2) * phi' * diag(wn.*lhn) * phi;
        
        gind = l2g(mx, my, 1:N+1, N+1);
        A(gind, gind) = A(gind, gind) + lA;
    end
    
    % left boundary
    mx = 1;
    for my = 1:My
        lxn = ones(N+1, 1) * mesh.x(mx);
        lyn = (xn + 1)/2 * hy(my) + mesh.y(my);
        lhn = arrayfun(prob.h, lxn, lyn);
        
        lA = (hy(my)/2) * phi' * diag(wn.*lhn) * phi;
        
        gind = l2g(mx, my, 1, 1:N+1);
        A(gind, gind) = A(gind, gind) + lA;
    end
    
    % right boundary
    mx = Mx;
    for my = 1:My
        lxn = ones(N+1, 1) * mesh.x(mx+1);
        lyn = (xn + 1)/2 * hy(my) + mesh.y(my);
        lhn = arrayfun(prob.h, lxn, lyn);
        
        lA = (hy(my)/2) * phi' * diag(wn.*lhn) * phi;
        
        gind = l2g(mx, my, N+1, 1:N+1);
        A(gind, gind) = A(gind, gind) + lA;
    end
    
    [U, lam] = eigs(A, B, num, 0);
    [~, ind] = sort(abs(diag(lam)));
    lam = diag(lam(ind, ind));
    U = U(:,ind);
    
elseif prob.bd == 'P'
    % TODO
end

end
