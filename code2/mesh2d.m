function mesh = mesh2d(x, y)

if nargin < 2
    Mx = x; My = x;
    x = linspace(0, 1, Mx +1);
    y = linspace(0, 1, My +1);
end
    
if length(x)==1 && length(y)==1
    Mx = x; My = y;
    x = linspace(0, 1, Mx +1);
    y = linspace(0, 1, My +1);
end

mesh = struct('x', x, 'y', y, ...
    'Mx', length(x)-1, 'My', length(y)-1, ...
    'hx', diff(x), 'hy', diff(y));