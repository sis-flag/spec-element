function U = solve2d(prob, mesh, N)

% defult input
if nargin < 3
    N = 6;
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
F = zeros((Mx*N+1)*(My*N+1), 1);

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
        lfn = arrayfun(prob.f, lxn, lyn);
        
        lA = (hy(my)/hx(mx)) * phix' * diag(wn2d.*la11n) * phix ...
            + phix' * diag(wn2d.*la12n) * phiy ...
            + phiy' * diag(wn2d.*la21n) * phix ...
            + (hx(mx)/hy(my)) * phiy' * diag(wn2d.*la22n) * phiy ...
            + (hy(my)/2) * phi0' * diag(wn2d.*lb1n) * phix ...
            + (hx(mx)/2) * phi0' * diag(wn2d.*lb2n) * phiy ...
            + (hx(mx)*hy(my)/4) * phi0' * diag(wn2d.*lcn) * phi0;
        
        lF = (hx(mx)*hy(my)/4) * phi0' * (lfn.*wn2d);
        
%         gind = l2g(mx, my, 1:(N+1)^2);
%         A(gind, gind) = A(gind, gind) + lA;
%         F(gind) = F(gind) + lF;
        
        gind = l2g(mx, my, 1:(N+1)^2);
        [gj, gi] = meshgrid(gind, gind);
        
        iA(nnzA+1: nnzA+(N+1)^4) = gi(:);
        jA(nnzA+1: nnzA+(N+1)^4) = gj(:);
        vA(nnzA+1: nnzA+(N+1)^4) = lA(:);
        nnzA = nnzA + (N+1)^4;
        
        F(gind) = F(gind) + lF;
    end
end

A = sparse(iA, jA, vA, (Mx*N+1)*(My*N+1), (Mx*N+1)*(My*N+1));

if prob.bd == 'D'
    
    lP0 = phi' * diag(wn) * phi;
    
    % bottom boundary
    my = 1;
    for mx = 1:Mx
        lxn = (xn + 1)/2 * hx(mx) + mesh.x(mx);
        lyn = ones(N+1, 1) * mesh.y(my);
        lu0n = arrayfun(prob.u0, lxn, lyn);
        
        lU = lP0 \ (phi' * (lu0n.*wn));
        
        gind = l2g(mx, my, 1:N+1, 1);
        A(gind, :) = 0;
        for k = 1:N+1
            A(gind(k), gind(k)) = 1;
        end
        F(gind) = lU;
    end
    
    % top boundary
    my = My;
    for mx = 1:Mx
        lxn = (xn + 1)/2 * hx(mx) + mesh.x(mx);
        lyn = ones(N+1, 1) * mesh.y(my+1);
        lu0n = arrayfun(prob.u0, lxn, lyn);
        
        lU = lP0 \ (phi' * (lu0n.*wn));
        
        gind = l2g(mx, my, 1:N+1, N+1);
        A(gind, :) = 0;
        for k = 1:N+1
            A(gind(k), gind(k)) = 1;
        end
        F(gind) = lU;
    end
    
    % left boundary
    mx = 1;
    for my = 1:My
        lxn = ones(N+1, 1) * mesh.x(mx);
        lyn = (xn + 1)/2 * hy(my) + mesh.y(my);
        lu0n = arrayfun(prob.u0, lxn, lyn);
        
        lU = lP0 \ (phi' * (lu0n.*wn));
        
        gind = l2g(mx, my, 1, 1:N+1);
        A(gind, :) = 0;
        for k = 1:N+1
            A(gind(k), gind(k)) = 1;
        end
        F(gind) = lU;
    end
    
    % right boundary
    mx = Mx;
    for my = 1:My
        lxn = ones(N+1, 1) * mesh.x(mx+1);
        lyn = (xn + 1)/2 * hy(my) + mesh.y(my);
        lu0n = arrayfun(prob.u0, lxn, lyn);
        
        lU = lP0 \ (phi' * (lu0n.*wn));
        
        gind = l2g(mx, my, N+1, 1:N+1);
        A(gind, :) = 0;
        for k = 1:N+1
            A(gind(k), gind(k)) = 1;
        end
        F(gind) = lU;
    end
    
    U = A \ F;
    
elseif prob.bd == 'R'
    
    % bottom boundary
    my = 1;
    for mx = 1:Mx
        lxn = (xn + 1)/2 * hx(mx) + mesh.x(mx);
        lyn = ones(N+1, 1) * mesh.y(my);
        lhn = arrayfun(prob.h, lxn, lyn);
        lgn = arrayfun(prob.g, lxn, lyn);
        
        lA = (hx(mx)/2) * phi' * diag(wn.*lhn) * phi;
        lF = (hx(mx)/2) * phi' * (lgn.*wn);
        
        gind = l2g(mx, my, 1:N+1, 1);
        A(gind, gind) = A(gind, gind) + lA;
        F(gind) = F(gind) + lF;
    end
    
    % top boundary
    my = My;
    for mx = 1:Mx
        lxn = (xn + 1)/2 * hx(mx) + mesh.x(mx);
        lyn = ones(N+1, 1) * mesh.y(my+1);
        lhn = arrayfun(prob.h, lxn, lyn);
        lgn = arrayfun(prob.g, lxn, lyn);
        
        lA = (hx(mx)/2) * phi' * diag(wn.*lhn) * phi;
        lF = (hx(mx)/2) * phi' * (lgn.*wn);
        
        gind = l2g(mx, my, 1:N+1, N+1);
        A(gind, gind) = A(gind, gind) + lA;
        F(gind) = F(gind) + lF;
    end
    
    % left boundary
    mx = 1;
    for my = 1:My
        lxn = ones(N+1, 1) * mesh.x(mx);
        lyn = (xn + 1)/2 * hy(my) + mesh.y(my);
        lhn = arrayfun(prob.h, lxn, lyn);
        lgn = arrayfun(prob.g, lxn, lyn);
        
        lA = (hy(my)/2) * phi' * diag(wn.*lhn) * phi;
        lF = (hy(my)/2) * phi' * (lgn.*wn);
        
        gind = l2g(mx, my, 1, 1:N+1);
        A(gind, gind) = A(gind, gind) + lA;
        F(gind) = F(gind) + lF;
    end
    
    % right boundary
    mx = Mx;
    for my = 1:My
        lxn = ones(N+1, 1) * mesh.x(mx+1);
        lyn = (xn + 1)/2 * hy(my) + mesh.y(my);
        lhn = arrayfun(prob.h, lxn, lyn);
        lgn = arrayfun(prob.g, lxn, lyn);
        
        lA = (hy(my)/2) * phi' * diag(wn.*lhn) * phi;
        lF = (hy(my)/2) * phi' * (lgn.*wn);
        
        gind = l2g(mx, my, N+1, 1:N+1);
        A(gind, gind) = A(gind, gind) + lA;
        F(gind) = F(gind) + lF;
    end
    
    U = A \ F;
    
elseif prob.bd == 'P'
    U = reshape(A \ F, Mx*N, My*N);
    U = [U, U(:,1)];
    U = [U; U(1,:)];
    U = reshape(U, (Mx*N+1)*(My*N+1), 1);
end

end
