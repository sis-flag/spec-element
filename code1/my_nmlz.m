function u = my_nmlz(u)
% normalize u such that |u| <= 1 and max |u| = 1
maxu = max(max(u));
minu = min(min(u));
if maxu < -minu
    u = u / minu;
else
    u = u / maxu;
end
