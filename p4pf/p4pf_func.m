function [min]=p4pf_func(M,m,fgt)

[f R t] = P4Pf_m(m, M);

min=9999;
% solutions test
for i=1:length(f)
    if min>abs(fgt-f(i))/fgt
        min=abs(fgt-f(i))/fgt;
    end
end