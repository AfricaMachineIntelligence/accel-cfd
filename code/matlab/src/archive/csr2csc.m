function [A] = csr2csc(v,c,r)

idi = zeros(1,length(c));

itv = 1;

for i = 1:1:length(r)
    
    for j = 1:1:length(r)
        
        if c(itv) == j
            idi(itv) = i;
            itv = itv+1;
        end
    
    end

end


A = sparse(idi,c,v,length(r),length(r));
