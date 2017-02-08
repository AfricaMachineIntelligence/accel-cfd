function [A,d] = make_sparse_matrix(Ap,Ae,Aw,An,As,At,Ab,bp)

[m,n,l] = size(Ap);

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

val = zeros(1,t_vars);
row = zeros(1,t_vars);
col = zeros(1,t_vars);

d = zeros(1,m*n*l);

% Initialize Counters
idv = 1;

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

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n
                        
                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    else
                        
                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;
                        
                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;

                    elseif j == n

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;

                    else

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                    elseif j == n

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                    else

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;

                    elseif j == n

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;

                    else

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
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

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    elseif j == n
                        
                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    else

                        val(idv) = -At(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idt;
                        idv = idv+1;

                        val(idv) = -An(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idn;
                        idv = idv+1;

                        val(idv) = -Aw(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idw;
                        idv = idv+1;

                        val(idv) = Ap(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idp;
                        idv = idv+1;

                        val(idv) = -Ae(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ide;
                        idv = idv+1;

                        val(idv) = -As(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = ids;
                        idv = idv+1;
                        
                        val(idv) = -Ab(j,i,k);
                        row(idv) = j+(i-1)*n+(k-1)*n*m;
                        col(idv) = idb;
                        idv = idv+1;

                    end
                end

            end

        end
    end
end

for k = 1:1:l
    for i = 1:1:m
        for j = 1:1:n
            d(j+(i-1)*n+(k-1)*n*m) = bp(j,i,k);
        end
    end
end

A = sparse(row,col,val);
return
