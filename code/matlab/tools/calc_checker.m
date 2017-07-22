close all
clear variables
clc

m = 10;
n = 10; 
l = 10;
indx = zeros(m,n,l);

for i = 1:1:m
    for j = 1:1:n
        for k = 1:1:l
            indx(i,j,k) = calc_index_3d(i,j,k,m,n,l);
        end
    end
end
