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

indx = zeros(m,n,l);

% Calculate indx Values
for i = 1:1:m
    for j = 1:1:n
        for k = 1:1:l
            indx(i,j,k) = calc_index_3d(i,j,k,m,n,l);
        end
    end
end

% -------------------------------------------------------------------------
% Fill Interior Nodes
for k = 2:1:l-1
    for j = 2:1:n-1
        for i = 2:1:m-1
            
            
            % Bottom Node
            val(indx(i,j,k)-3) = -Ab(i,j,k);
            row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
            col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);
            
            % South Node
            val(indx(i,j,k)-2) = -As(i,j,k); 
            row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
            col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);
            
            % West Node
            val(indx(i,j,k)-1) = -Aw(i,j,k);
            row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
            col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 
            
            % Center Node (p)
            val(indx(i,j,k)) = Ap(i,j,k);
            row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
            col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
            
            % East Node
            val(indx(i,j,k)+1) = -Ae(i,j,k); 
            row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
            col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);
            
            % North Node
            val(indx(i,j,k)+2) = -An(i,j,k);
            row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
            col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);
            
            % Top Node
            val(indx(i,j,k)+3) = -At(i,j,k);
            row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
            col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
            
            % Solution Value
            d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
            
        end
    end
end
% -------------------------------------------------------------------------
% Fill West Wall (i = 1)
for k = 2:1:l-1
    for j = 2:1:n-1
        
        i = 1;
        
        % Bottom Node
        val(indx(i,j,k)-2) = -Ab(i,j,k);
        row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

        % South Node
        val(indx(i,j,k)-1) = -As(i,j,k); 
        row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1);

        % Center Node (p)
        val(indx(i,j,k)) = Ap(i,j,k);
        row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

        % East Node
        val(indx(i,j,k)+1) = -Ae(i,j,k); 
        row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

        % North Node
        val(indx(i,j,k)+2) = -An(i,j,k);
        row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

        % Top Node
        val(indx(i,j,k)+3) = -At(i,j,k);
        row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
        
        % Solution Value
        d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
        
    end
end
% -------------------------------------------------------------------------
% Fill East Wall (i = m)
for k = 2:1:l-1
    for j = 2:1:n-1
            
        i = m;
        
        % Bottom Node
        val(indx(i,j,k)-3) = -Ab(i,j,k);
        row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);

        % South Node
        val(indx(i,j,k)-2) = -As(i,j,k); 
        row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

        % West Node
        val(indx(i,j,k)-1) = -Aw(i,j,k);
        row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
        col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

        % Center Node (p)
        val(indx(i,j,k)) = Ap(i,j,k);
        row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

        % North Node
        val(indx(i,j,k)+1) = -An(i,j,k);
        row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1);

        % Top Node
        val(indx(i,j,k)+2) = -At(i,j,k);
        row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
        % Solution Value
        d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
            
    end
end
% -------------------------------------------------------------------------
% Fill South Wall (j = 1)
for k = 2:1:l-1
    for i = 2:1:m-1
        
        j = 1;
        
        % Bottom Node
        val(indx(i,j,k)-2) = -Ab(i,j,k);
        row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

        % West Node
        val(indx(i,j,k)-1) = -Aw(i,j,k);
        row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
        col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

        % Center Node (p)
        val(indx(i,j,k)) = Ap(i,j,k);
        row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

        % East Node
        val(indx(i,j,k)+1) = -Ae(i,j,k); 
        row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

        % North Node
        val(indx(i,j,k)+2) = -An(i,j,k);
        row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

        % Top Node
        val(indx(i,j,k)+3) = -At(i,j,k);
        row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
            
        % Solution Value
        d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);

    end
end
% -------------------------------------------------------------------------
% Fill North Wall (j = n)
for k = 2:1:l-1
    for i = 2:1:m-1
        
        j = n;

        % Bottom Node
        val(indx(i,j,k)-3) = -Ab(i,j,k);
        row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);

        % South Node
        val(indx(i,j,k)-2) = -As(i,j,k); 
        row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

        % West Node
        val(indx(i,j,k)-1) = -Aw(i,j,k);
        row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
        col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

        % Center Node (p)
        val(indx(i,j,k)) = Ap(i,j,k);
        row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

        % East Node
        val(indx(i,j,k)+1) = -Ae(i,j,k); 
        row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

        % Top Node
        val(indx(i,j,k)+2) = -At(i,j,k);
        row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
        % Solution Value
        d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);

    end
end
% -------------------------------------------------------------------------
% Fill Bottom Wall (k = 1)
for j = 2:1:n-1
    for i = 2:1:m-1

        k = 1;
            
        % South Node
        val(indx(i,j,k)-2) = -As(i,j,k); 
        row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

        % West Node
        val(indx(i,j,k)-1) = -Aw(i,j,k);
        row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
        col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

        % Center Node (p)
        val(indx(i,j,k)) = Ap(i,j,k);
        row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

        % East Node
        val(indx(i,j,k)+1) = -Ae(i,j,k); 
        row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

        % North Node
        val(indx(i,j,k)+2) = -An(i,j,k);
        row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

        % Top Node
        val(indx(i,j,k)+3) = -At(i,j,k);
        row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
            
        % Solution Value
        d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);

    end
end
% -------------------------------------------------------------------------
% Fill Top Wall (k = l)
for j = 2:1:n-1
    for i = 2:1:m-1

        k = l;

        % Bottom Node
        val(indx(i,j,k)-3) = -Ab(i,j,k);
        row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);

        % South Node
        val(indx(i,j,k)-2) = -As(i,j,k); 
        row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

        % West Node
        val(indx(i,j,k)-1) = -Aw(i,j,k);
        row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
        col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

        % Center Node (p)
        val(indx(i,j,k)) = Ap(i,j,k);
        row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

        % East Node
        val(indx(i,j,k)+1) = -Ae(i,j,k); 
        row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

        % North Node
        val(indx(i,j,k)+2) = -An(i,j,k);
        row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
        col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);
            
        % Solution Value
        d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);

    end
end
% -------------------------------------------------------------------------
% Fill in Corner (W,S,B)
i = 1;
j = 1;
k = 1;

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

% East Node
val(indx(i,j,k)+1) = -Ae(i,j,k); 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

% North Node
val(indx(i,j,k)+2) = -An(i,j,k);
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

% Top Node
val(indx(i,j,k)+3) = -At(i,j,k);
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
            
% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Corner (E,S,B)
i = m;
j = 1;
k = 1;

% West Node
val(indx(i,j,k)-1) = -Aw(i,j,k);
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

% North Node
val(indx(i,j,k)+1) = -An(i,j,k);
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1);

% Top Node
val(indx(i,j,k)+2) = -At(i,j,k);
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Corner (W,N,B)
i = 1;
j = n;
k = 1;

% South Node
val(indx(i,j,k)-1) = -As(i,j,k); 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1);

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

% East Node
val(indx(i,j,k)+1) = -Ae(i,j,k); 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

% Top Node
val(indx(i,j,k)+2) = -At(i,j,k);
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Corner (E,N,B)
i = m;
j = n;
k = 1;

% South Node
val(indx(i,j,k)-2) = -As(i,j,k); 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

% West Node
val(indx(i,j,k)-1) = -Aw(i,j,k);
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

% Top Node
val(indx(i,j,k)+1) = -At(i,j,k);
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k);
            
% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Corner (W,S,T)
i = 1;
j = 1;
k = l;

% Bottom Node
val(indx(i,j,k)-1) = -Ab(i,j,k);
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-2);

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

% East Node
val(indx(i,j,k)+1) = -Ae(i,j,k); 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

% North Node
val(indx(i,j,k)+2) = -An(i,j,k);
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Corner (E,S,T)
i = m;
j = 1;
k = l;

% Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k);
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

% West Node
val(indx(i,j,k)-1) = -Aw(i,j,k);
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

% North Node
val(indx(i,j,k)+1) = -An(i,j,k);
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1);
            
% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Corner (W,N,T)
i = 1;
j = n;
k = l;

% Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k);
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

% South Node
val(indx(i,j,k)-1) = -As(i,j,k); 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1);

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

% East Node
val(indx(i,j,k)+1) = -Ae(i,j,k); 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);
            
% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Corner (E,N,T)
i = m;
j = n;
k = l;

% Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k);
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);

% South Node
val(indx(i,j,k)-2) = -As(i,j,k); 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

% West Node
val(indx(i,j,k)-1) = -Aw(i,j,k);
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

% Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k);
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
            
% Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
% -------------------------------------------------------------------------
% Fill in Line (W2E - S,B)
for i = 2:1:m-1
    j = 1;
    k = 1;

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+2) = -An(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

    % Top Node
    val(indx(i,j,k)+3) = -At(i,j,k);
    row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
end
% -------------------------------------------------------------------------
% Fill in Line (W2E - N,B)
for i = 2:1:m-1
    j = n;
    k = 1;

    % South Node
    val(indx(i,j,k)-2) = -As(i,j,k); 
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);
    
    % Top Node
    val(indx(i,j,k)+2) = -At(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
end
% -------------------------------------------------------------------------
% Fill in Line (W2E - S,T)
for i = 2:1:m-1
    j = 1;
    k = l;

    % Bottom Node
    val(indx(i,j,k)-2) = -Ab(i,j,k);
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+2) = -An(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
end
% -------------------------------------------------------------------------
% Fill in Line (W2E - N,T)
for i = 2:1:m-1
    j = n;
    k = l;
    
    % Bottom Node
    val(indx(i,j,k)-3) = -Ab(i,j,k);
    row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);

    % South Node
    val(indx(i,j,k)-2) = -As(i,j,k); 
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
end
% -------------------------------------------------------------------------
% Fill in Line (S2N - W,B)
for j = 2:1:n-1
    i = 1;
    k = 1;

    % South Node
    val(indx(i,j,k)-1) = -As(i,j,k); 
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1);

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+2) = -An(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

    % Top Node
    val(indx(i,j,k)+3) = -At(i,j,k);
    row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------
% Fill in Line (S2N - E,B)
for j = 2:1:n-1
    i = m;
    k = 1;

    % South Node
    val(indx(i,j,k)-2) = -As(i,j,k); 
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+1) = -An(i,j,k);
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1);

    % Top Node
    val(indx(i,j,k)+2) = -At(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------
% Fill in Line (S2N - W,T)
for j = 2:1:n-1
    i = 1;
    k = l;

    % Bottom Node
    val(indx(i,j,k)-2) = -Ab(i,j,k);
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

    % South Node
    val(indx(i,j,k)-1) = -As(i,j,k); 
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1);

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+2) = -An(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------
% Fill in Line (S2N - E,T)
for j = 2:1:n-1
    i = m;
    k = l;
    
    % Bottom Node
    val(indx(i,j,k)-3) = -Ab(i,j,k);
    row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);

    % South Node
    val(indx(i,j,k)-2) = -As(i,j,k); 
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+1) = -An(i,j,k);
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------
% Fill in Line (B2T - W,S)
for k = 2:1:l-1
    i = 1;
    j = 1;

    % Bottom Node
    val(indx(i,j,k)-1) = -Ab(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-2);

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+2) = -An(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1);

    % Top Node
    val(indx(i,j,k)+3) = -At(i,j,k);
    row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------
% Fill in Line (B2T - E,S)
for k = 2:1:l-1
    i = m;
    j = 1;
    
    % Bottom Node
    val(indx(i,j,k)-2) = -Ab(i,j,k);
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % North Node
    val(indx(i,j,k)+1) = -An(i,j,k);
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1);

    % Top Node
    val(indx(i,j,k)+2) = -At(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------
% Fill in Line (B2T - W,N)
for k = 2:1:l-1
    i = 1;
    j = n;
    
    % Bottom Node
    val(indx(i,j,k)-2) = -Ab(i,j,k);
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2);

    % South Node
    val(indx(i,j,k)-1) = -As(i,j,k); 
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1);

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % East Node
    val(indx(i,j,k)+1) = -Ae(i,j,k); 
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1);

    % Top Node
    val(indx(i,j,k)+2) = -At(i,j,k);
    row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------
% Fill in Line (B2T - E,N)
for k = 2:1:l-1
    i = m;
    j = n;
    
    % Bottom Node
    val(indx(i,j,k)-3) = -Ab(i,j,k);
    row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2);

    % South Node
    val(indx(i,j,k)-2) = -As(i,j,k); 
    row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1);

    % West Node
    val(indx(i,j,k)-1) = -Aw(i,j,k);
    row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1); 
    col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1); 

    % Center Node (p)
    val(indx(i,j,k)) = Ap(i,j,k);
    row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1);

    % Top Node
    val(indx(i,j,k)+1) = -At(i,j,k);
    row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1);
    col(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k);
            
    % Solution Value
    d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k);
    
end
% -------------------------------------------------------------------------

A = sparse(row,col,val);
return
