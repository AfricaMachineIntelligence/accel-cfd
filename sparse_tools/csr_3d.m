% CSR Matrix Creator for 3D Heat Transfer Problems
% Matt Blomquist - Feb 5, 2017

close all
clear all
clc

% Set 3D Grid Geometry
m = 100;
n = 100;
l = 100;

% Determine Number of Nodes
i_nodes = (m-2)*(n-2)*(l-2);
l_nodes = 4*((m-2)+(n-2)+(l-2));
s_nodes = 2*((m-2)*(n-2)+(n-2)*(l-2)+(l-2)*(m-2));
c_nodes = 8;

% Calculate Total A matrix Variables
i_vars = i_nodes*7;
s_vars = s_nodes*6;
l_vars = l_nodes*5;
c_vars = c_nodes*4;

t_vars = i_vars+l_vars+s_vars+c_vars;

n_perc = t_vars/((m*n*l)^2);

line1_text = 'Total number of A matrix (3D) vars for %3.0f is : %3.0f.\n';
line2_text = 'Sparsity is: %0.5f\n';

fprintf(line1_text,(m*n*l)^2,t_vars)
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
at_i = .1;
ab_i = .1;
Su_i = 0;
ap_i = aw_i*6+Su_i;

% Initialize CSR Arrays
val = zeros(1,t_vars);
col = zeros(1,t_vars);
row = zeros(1,m*n*l);

% Initialize Counters
idv = 1;
idr = 1;

% Start Depth Loop
for k = 1:1:l

    % Start k = 1 loop
    if k == 1
        
        % Start Row Sweep
        for i = 1:1:m

            % North Row Sweep
            if i == 1

                % West to East Sweep
                for j = 1:1:n

                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    ids = (k-1)*m*n+(i)*m+j;
                    idb = (k)*m*n+(i-1)*m+j;

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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    end
                end

            % South Row Sweep
            elseif i == m

                % West to East Sweep
                for j = 1:1:n

                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    idn = (k-1)*m*n+(i-2)*m+j;
                    idb = (k)*m*n+(i-1)*m+j;


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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    end
                end

            % Top to Bottom Sweep
            else

                % West to East Sweep
                for j = 1:1:n

                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    idn = (k-1)*m*n+(i-2)*m+j;
                    ids = (k-1)*m*n+(i)*m+j;
                    idb = (k)*m*n+(i-1)*m+j;

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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    end
                end

            end

        end
    
    % Start k = l loop
    elseif k == l
        
        % Start Row Sweep
        for i = 1:1:m

            % North Row Sweep
            if i == 1

                % West to East Sweep
                for j = 1:1:n

                    idt = (k-2)*m*n+(i-1)*m+j;
                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    ids = (k-1)*m*n+(i)*m+j;

                    if j == 1

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
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

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
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

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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

                    idt = (k-2)*m*n+(i-1)*m+j;
                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    idn = (k-1)*m*n+(i-2)*m+j;

                    if j == 1

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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

                    idt = (k-2)*m*n+(i-1)*m+j;
                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    idn = (k-1)*m*n+(i-2)*m+j;
                    ids = (k-1)*m*n+(i)*m+j;

                    if j == 1

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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
        
    % Start k 2:1:l-1 loop
    else
        
        % Start Row Sweep
        for i = 1:1:m

            % North Row Sweep
            if i == 1

                % West to East Sweep
                for j = 1:1:n

                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    ids = (k-1)*m*n+(i)*m+j;
                    idt = (k-2)*m*n+(i-1)*m+j;
                    idb = (k)*m*n+(i-1)*m+j;

                    if j == 1

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    end
                end

            % South Row Sweep
            elseif i == m

                % West to East Sweep
                for j = 1:1:n

                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    idn = (k-1)*m*n+(i-2)*m+j;
                    idt = (k-2)*m*n+(i-1)*m+j;
                    idb = (k)*m*n+(i-1)*m+j;

                    if j == 1

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -an_i;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = ap_i;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -ae_i;
                        col(idv) = ide;
                        idv = idv+1;
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -an_i;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -aw_i;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = ap_i;
                        col(idv) = idp;
                        idv = idv+1;
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;


                    end
                end

            % Top to Bottom Sweep
            else

                % West to East Sweep
                for j = 1:1:n

                    idp = (k-1)*m*n+(i-1)*m+j;
                    idw = idp-1;
                    ide = idp+1;
                    idn = (k-1)*m*n+(i-2)*m+j;
                    ids = (k-1)*m*n+(i)*m+j;
                    idt = (k-2)*m*n+(i-1)*m+j;
                    idb = (k)*m*n+(i-1)*m+j;

                    if j == 1

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        row(idr) = idt;
                        idr = idr+1;

                        val(idv) = -at_i;
                        col(idv) = idt;
                        idv = idv+1;

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
                        
                        val(idv) = -ab_i;
                        col(idv) = idb;
                        idv = idv+1;

                    end
                end

            end

        end
    end
end

% Generate Solution Array

d = zeros(1,m*n*l);

for z = 1:1:m*n*l
    d(z) = rand;
end

% % Write CSV to File
% display('Writing to file')
% 
% format_title = 'csr_%dx%dx%d_val.txt';
% title = sprintf(format_title,m,n,l);
% formatSpec = '%4.8f\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,val);
% fclose(fileID);
% 
% format_title = 'csr_%dx%dx%d_col.txt';
% title = sprintf(format_title,m,n,l);
% formatSpec = '%d\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,col);
% fclose(fileID);
% 
% format_title = 'csr_%dx%dx%d_row.txt';
% title = sprintf(format_title,m,n,l);
% formatSpec = '%d\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,row);
% fclose(fileID);
% 
% format_title = 'csr_%dx%dx%d_d.txt';
% title = sprintf(format_title,m,n,l);
% formatSpec = '%4.8f\n';
% fileID = fopen(title,'w');
% fprintf(fileID,formatSpec,d);
% fclose(fileID);
% 
% display('All done')
