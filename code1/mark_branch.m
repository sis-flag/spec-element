function mark = mark_branch(V)
% find connected branches of V == 0
% input:
%     V(2-d array): elements are all 0 or 1
% output:
%     mark(2-d array): different branches are marked by different integers

M = size(V, 1);
mark = zeros(size(V));
num_mark = 1;

    % DFS mark all connected elements
    function set_mark(k1, k2)
        if mark(k1, k2) > 0
            return
        end
        
        mark(k1, k2) = num_mark;
        
        if k1 > 1 && V(k1-1, k2) == 0
            set_mark(k1-1, k2);
        end
        if k1 < M && V(k1+1, k2) == 0
            set_mark(k1+1, k2);
        end
        if k2 > 1 && V(k1, k2-1) == 0
            set_mark(k1, k2-1);
        end
        if k2 < M && V(k1, k2+1) == 0
            set_mark(k1, k2+1);
        end
    end

% loop all V=0 elements
[ind1, ind2] = find(V == 0);
for k = 1:length(ind1)
    k1 = ind1(k); k2 = ind2(k);
    
    % skip marked element
    if mark(k1, k2) > 0
        continue
    else
        set_mark(k1, k2);
        num_mark = num_mark + 1;
    end
end

end