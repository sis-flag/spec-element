function prob = test_prob2d_1()


% boundary codition type
% 'D' for Dirichlet, 'R' for Robin, 'P' for periodic
bd = 'R';
sx = 0.3; sy = 0.6;

    function a = a11(~, ~)
        a = 1;
    end

    function a = a12(~, ~)
        a = 0;
    end

    function a = a22(~, ~)
        a = 1;
    end

    function b = b1(x, y)
        b = exp(x) + cos(y);
    end

    function b = b2(x, y)
        b = sin(x) + log(1+y);
    end

    function c = c(x, y)
        c = 1 + x + y;
    end

    function u = u(x, y)
        u = sin(pi*x + sx) .* sin(pi*y + sy);
    end

    function ux = dux(x, y)
        ux = pi * cos(pi*x + sx) .* sin(pi*y + sy);
    end

    function uy = duy(x, y)
        uy = pi * sin(pi*x + sx) .* cos(pi*y + sy);
    end

    function f = f(x, y)
        f = (2*pi*pi + c(x,y)) * u(x,y) ...
          + b1(x,y) * dux(x,y) + b2(x,y) * duy(x,y);
    end

    function h = h(x, y)
        h = sin(1+x) * exp(y);
    end

    function g = g(x, y)
        if x == 0
            g = -dux(x,y) + h(x,y) * u(x,y);
        elseif x == 1
            g = dux(x,y) + h(x,y) * u(x,y);
        elseif y == 0
            g = -duy(x,y) + h(x,y) * u(x,y);
        elseif y == 1
            g = duy(x,y) + h(x,y) * u(x,y);
        end
    end

prob = struct(...
    'a11', @a11, 'a12', @a12, 'a21', @a12, 'a22', @a22,...
    'b1', @b1, 'b2', @b2, 'c', @c, 'f', @f, ...
    'u0', @u, 'h', @h, 'g', @g, 'bd', bd, ...
    'u', @u, 'dux', @dux, 'duy', @duy);
end
