function [A] = csr2dense(v,c,r)

A = zeros(length(r),length(r));

itv = 1;

for i = 1:1:length(r)
    
    for j = 1:1:length(r)
        
        if c(itv) == j
            A(i,j) = v(itv);
            itv = itv+1;
        else
            A(i,j) = 0;
        end
        
    end

end
