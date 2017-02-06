% CSR Matrix Creator for 2D Heat Transfer Problems
% Matt Blomquist - Feb 5, 2017

close all
clear all
clc

% Set grid geo
m = 100;
n = 100; 

% Determine Number of Variables
i_nodes = (m-2)*(n-2);
s_nodes = 2*((n-2)+(m-2));
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

% Create coefficients (interior)
aw_i = .1;
ae_i = .1;
as_i = .1;
an_i = .1;
Su_i = 0;
ap_i = aw_i*4+Su_i;

% Initialize CSR Arrays
val = zeros(1,t_vars);
col = zeros(1,t_vars);
row = zeros(1,m*n);

% Initialize Counters
idv = 1;
idr = 1;

% Start Row Loop
for i = 1:1:m
    
    % North Row Sweep
    if i == 1
        
        % West to East Sweep
        for j = 1:1:n
            
            idp = (i-1)*m+j;
            idw = idp-1;
            ide = idp+1;
            ids = (i)*m+j;
            
            if j == 1
                
                row(idr) = idp;
                idr = idr+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = -ae_i;
                col(idv) = ide;
                idv = idv+1;
                
                val(idv) = -as_i;
                col(idv) = ids;
                idv = idv+1;
                
            elseif j == n
                
                row(idr) = idw;
                idr = idr+1;
                
                val(idv) = -aw_i;
                col(idv) = idw;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = -as_i;
                col(idv) = ids;
                idv = idv+1;
                
            else
                
                row(idr) = idw;
                idr = idr+1;
                
                val(idv) = -aw_i;
                col(idv) = idw;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = ae_i;
                col(idv) = ide;
                idv = idv+1;
                
                val(idv) = -as_i;
                col(idv) = ids;
                idv = idv+1;
                
            end
        end
        
    % South Row Sweep
    elseif i == m
        
        % West to East Sweep
        for j = 1:1:n
            
            idp = (i-1)*m+j;
            idw = idp-1;
            ide = idp+1;
            idn = (i-2)*m+j;
            
            if j == 1
                
                row(idr) = idn;
                idr = idr+1;
                
                val(idv) = -an_i;
                col(idv) = idn;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = -ae_i;
                col(idv) = ide;
                idv = idv+1;
                              
            elseif j == n
                
                row(idr) = idn;
                idr = idr+1;
                
                val(idv) = -an_i;
                col(idv) = idn;
                idv = idv+1;
                
                val(idv) = -aw_i;
                col(idv) = idw;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                               
            else
                
                row(idr) = idn;
                idr = idr+1;
                
                val(idv) = -an_i;
                col(idv) = idn;
                idv = idv+1;
                
                val(idv) = -aw_i;
                col(idv) = idw;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = -ae_i;
                col(idv) = ide;
                idv = idv+1;
                
                
            end
        end
        
    % Top to Bottom Sweep
    else
        
        % West to East Sweep
        for j = 1:1:n
            
            idp = (i-1)*m+j;
            idw = idp-1;
            ide = idp+1;
            idn = (i-2)*m+j;
            ids = (i)*m+j;
            
            if j == 1
                
                row(idr) = idn;
                idr = idr+1;
                
                val(idv) = -an_i;
                col(idv) = idn;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = -ae_i;
                col(idv) = ide;
                idv = idv+1;
                
                val(idv) = -as_i;
                col(idv) = ids;
                idv = idv+1;
                              
            elseif j == n
                
                row(idr) = idn;
                idr = idr+1;
                
                val(idv) = -an_i;
                col(idv) = idn;
                idv = idv+1;
                
                val(idv) = -aw_i;
                col(idv) = idw;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = -as_i;
                col(idv) = ids;
                idv = idv+1;
                               
            else
                
                row(idr) = idn;
                idr = idr+1;
                
                val(idv) = -an_i;
                col(idv) = idn;
                idv = idv+1;
                
                val(idv) = -aw_i;
                col(idv) = idw;
                idv = idv+1;
                
                val(idv) = ap_i;
                col(idv) = idp;
                idv = idv+1;
                
                val(idv) = -ae_i;
                col(idv) = ide;
                idv = idv+1;
                
                val(idv) = -as_i;
                col(idv) = ids;
                idv = idv+1;
                
            end
        end
        
    end
    
end

% Generate Solution Array

d = zeros(1,m*n);

for z = 1:1:m*n
    d(z) = rand;
end

% % Write CSV to File
% display('Writing to file')
% 
% format_title = 'csr_%dx%d_val.txt';
% title = sprintf(format_title,m,n);
% formatSpec = '%4.8f\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,val);
% fclose(fileID);
% 
% format_title = 'csr_%dx%d_col.txt';
% title = sprintf(format_title,m,n);
% formatSpec = '%d\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,col);
% fclose(fileID);
% 
% format_title = 'csr_%dx%d_row.txt';
% title = sprintf(format_title,m,n);
% formatSpec = '%d\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,row);
% fclose(fileID);
% 
% format_title = 'csr_%dx%d_d.txt';
% title = sprintf(format_title,m,n);
% formatSpec = '%4.8f\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,d);
% fclose(fileID);
% 
% display('All done')
