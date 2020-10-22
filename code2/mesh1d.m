function mesh = mesh1d(x)

if length(x)==1
    M = x;
    x = linspace(0, 1, M+1);
end

mesh = struct('x', x, 'M', length(x)-1, 'h', diff(x));
