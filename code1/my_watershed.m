function mark = my_watershed(w)

N = size(w,1);
mark = zeros(size(w));
num_mark = 1;

[~, ind] = sort(w(:));
for k = 1:length(ind)
    indi = mod(ind(k)-1, N) +1;
    indj = ceil(ind(k) / N);
    
    if indi == 1 || indi == N || indj == 1 || indj == N
        continue
    end
    
    mloc = mark(indi-1:indi+1, indj-1:indj+1);
    
    mloc = mloc(:);
    mloc(mloc == 0) = [];
    mloc(isnan(mloc)) = [];
    mrange = range(mloc);
    
    if isempty(mloc)
        mark(indi, indj) = num_mark + 1;
        num_mark = num_mark + 1;
    elseif mrange == 0
        mark(indi, indj) = mloc(1);
    else
        mark(indi, indj) = nan;
    end
    
%     if mod(k, floor(N*N/10)) == 0
%         figure();
%         hold on
%         s = pcolor(x1, x2, mark);
%         s.LineStyle = 'none';
%         colorbar;
%     end
end