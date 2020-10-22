function prob = test_prob2d_2()

% boundary codition type
% 'D' for Dirichlet, 'R' for Robin, 'P' for periodic
bd = 'D';

    function a = a11(x, ~)
        a = 1 + cos(x);
    end

    function a = a12(x, y)
        a = x * y;
    end

    function a = a22(~, y)
        a = 1 + exp(y);
    end

    function b = b1(x, ~)
        b = exp(x);
    end

    function b = b2(~, y)
        b = cos(y);
    end

    function c = c(x, y)
        c = 1 + x + y;
    end

    function u = u(x, y)
        u = sin(pi*x) .* sin(pi*y);
    end

    function f = f(x, y)
        f = sin(pi*x)*sin(pi*y)*(x+y+1)...
          - x*pi*cos(pi*x)*sin(pi*y)...
          - y*pi*cos(pi*y)*sin(pi*x)...
          + pi^2*sin(pi*x)*sin(pi*y)*(cos(x) + 1)...
          + pi^2*sin(pi*x)*sin(pi*y)*(exp(y) + 1)...
          + pi*cos(pi*y)*sin(pi*x)*cos(y)...
          + pi*exp(x)*cos(pi*x)*sin(pi*y)...
          - pi*exp(y)*cos(pi*y)*sin(pi*x)...
          + pi*cos(pi*x)*sin(pi*y)*sin(x)...
          - 2*x*y*pi^2*cos(pi*x)*cos(pi*y);
    end

prob = struct(...
    'a11', @a11, 'a12', @a12, 'a21', @a12, 'a22', @a22,...
    'b1', @b1, 'b2', @b2, 'c', @c, 'f', @f, ...
    'u', @u, 'u0', @u, 'bd', bd);
end
