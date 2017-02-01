% CSR Matrix Creator for 2D Heat Transfer Problems
% Matt Blomquist - Jan 30, 2017

close all
clear all
clc

% Set grid geo
m = 100;
n = 100; 

% Determine Number of Variables
i_nodes = (m-2)*(n-2);
s_nodes = 4*((n-2)+(m-2));
c_nodes = 4;

% Calculate Total A matrix Variables
i_vars = i_nodes*5;
s_vars = s_nodes*4;
c_vars = c_nodes*3;

t_vars = i_vars+s_vars+c_vars;

n_perc = t_vars/((m*n)^2);

line1_text = 'Total number of A matrix (2D) vars for %3.0f is : %3.0f.\n';
line2_text = 'Sparsity is: %0.5f\n';

fprintf(line1_text,(m*n)^2,t_vars)
fprintf(line2_text,n_perc)

% Set BC
bcw = 1;
bce = 0;
bcn = 0;
bcs = 0;


% Create T matrix
T = zeros(m,n);

% Create A matrix
A = zeros(m*n,m*n);

% Create coefficients (interior)
aw_i = .1;
ae_i = .1;
as_i = .1;
an_i = .1;
Su_i = 0;
ap_i = aw_i*4+Su_i;

% i = 1 sweep (top row)
i = 1;

for j = 1:1:n
    
    idp = j;
    idw = idp-1;
    ide = idp+1;
    ids = m+j;
    
    if j == 1
        A(idp,idp) = ap_i;
        A(idp,ide) = -ae_i;
        A(idp,ids) = -as_i;
    elseif j == n
        
        A(idp,idw) = -aw_i;
        A(idp,idp) = ap_i;
        A(idp,ids) = -as_i;      
        
    else
        A(idp,idw) = -aw_i;
        A(idp,idp) = ap_i;
        A(idp,ide) = -ae_i;
        A(idp,ids) = -as_i;    
    end
    
    
end 
% j = 1 sweep
j = 1;

for i = 1:1:m
    
    idp = (i-1)*m+j;
    ide = idp+j;
    idn = (i-2)*m+j;
    ids = (i)*m+j;
    
    if i == 1
       A(idp,idp) = ap_i;
       A(idp,ide) = -ae_i;
       A(idp,ids) = -as_i;
    elseif i == m
       A(idp,idp) = ap_i;
       A(idp,ide) = -ae_i;
       A(idp,idn) = -an_i;
       
    else
       A(idp,idp) = ap_i;
       A(idp,ide) = -ae_i;
       A(idp,idn) = -an_i;
       A(idp,ids) = -as_i;
    end
end 

% Interior Nodes
for j = 2:1:n-1
    for i = 2:1:m-1
        idp = (j-1)*m+i;
        idw = idp-1;
        ide = idp+1;
        idn = (j-2)*m+i;
        ids = (j)*m+i;
        
        A(idp,idn) = -an_i;
        A(idp,idw) = -aw_i;
        A(idp,idp) = ap_i;
        A(idp,ide) = -ae_i;
        A(idp,ids) = -as_i;
    end
end

% j = n sweep
j = n;

for i = 1:1:m
    
    idp = (i-1)*m+j;
    idw = idp-1;
    idn = (i-2)*m+j;
    ids = (i)*m+j;
    
    if i == 1
       A(idp,idp) = ap_i;
       A(idp,idw) = -aw_i;
       A(idp,ids) = -as_i;
    elseif i == m
       A(idp,idp) = ap_i;
       A(idp,idw) = -aw_i;
       A(idp,idn) = -an_i;
       
    else
       A(idp,idp) = ap_i;
       A(idp,idw) = -aw_i;
       A(idp,idn) = -an_i;
       A(idp,ids) = -as_i;
    end
end 
    
% i = m sweep (top row)
i = m;

for j = 1:1:n
    
    idp = (n-1)*m+j;
    idw = idp-1;
    ide = idp+1;
    idn = (n-2)*m+j;
    
    if j == 1
        A(idp,idp) = ap_i;
        A(idp,ide) = -ae_i;
        A(idp,idn) = -an_i;
    elseif j == n
        
        A(idp,idw) = -aw_i;
        A(idp,idp) = ap_i;
        A(idp,idn) = -an_i;      
        
    else
        A(idp,idw) = -aw_i;
        A(idp,idp) = ap_i;
        A(idp,ide) = -ae_i;
        A(idp,idn) = -an_i;    
    end
    
    
end 

% Create CSR
count = 1;
count2 = 1;
first_check = 1;

for j = 1:1:m*n
    
    for i = 1:1:m*n
        
        if A(i,j) ~= 0
            
            val(count) = A(i,j);
            col(count) = i;
            
            if first_check == 1
                row(count2) = i;
                
                first_check = 2;
                count2 = count2+1;
            end
            
            count = count+1;
            
        end
        
    end
    
    first_check = 1;
    
end

d = zeros(1,n*m);

% Create random d matrix
%for z = 1:1:n*m
%    d(z) = 3;
%end

for z = 1:1:m
    d(z) = 1;
end

for z = m*n:-1:m*n-m
    d(z) = 1;
end

% Create Solution File
x = d/A;
itr = 0;

% Plot T field
for i = 1:1:m
    for j = 1:1:n
        T(i,j) = x((i-1)*m+j);
    end
end

contourf(T,10)
grid minor

% Write CSV to File
display('Writing to file')

formatSpec = '%4.8f\n';
fileID = fopen('csr_200_val.txt','w');
fprintf(fileID,formatSpec,val);
fclose(fileID);

formatSpec = '%8.0f\n';
fileID = fopen('csr_200_col.txt','w');
fprintf(fileID,formatSpec,col);
fclose(fileID);

formatSpec = '%8.0f\n';
fileID = fopen('csr_200_row.txt','w');
fprintf(fileID,formatSpec,row);
fclose(fileID);

formatSpec = '%4.8f\n';
fileID = fopen('csr_200_x.txt','w');
fprintf(fileID,formatSpec,x);
fclose(fileID);

formatSpec = '%4.8f\n';
fileID = fopen('csr_200_d.txt','w');
fprintf(fileID,formatSpec,d);
fclose(fileID);

display('All done')
