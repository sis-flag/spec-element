function prob = test_prob1d(bd, s)

% boundary codition type
% 'D' for Dirichlet, 'R' for Robin, 'P' for periodic

% bd = 'D';
% s = 0;
% bd = 'R';
% s = 1;

    function a = a(x)
        a = 1 + sin(x);
    end

    function b = b(x)
        b = cos(x);
    end

    function c = c(x)
        c = 1 + x;
    end

    function f = f(x)
        f = (pi*pi + 1 + pi*pi*sin(x) + x).*sin(pi*x + s);
    end

    function u = u(x)
        u = sin(pi*x + s);
    end

    function du = du(x)
        du = pi*cos(pi*x + s);
    end

    function h = h(x)
        if x == 0
            h = pi;
        elseif x == 1
            h = -pi;
        end
    end

    function g = g(x)
        if x == 0
            g = -a(0)*pi*cos(s) + h(0)*sin(s);
        elseif x == 1
            g = -a(1)*pi*cos(s) - h(1)*sin(s);
        end
    end

prob = struct(...
    'a', @a, 'b', @b, 'c', @c, 'f', @f, ...
    'u0', @u, 'h', @h, 'g', @g, 'bd', bd, ...
    'u', @u, 'du', @du);
end
