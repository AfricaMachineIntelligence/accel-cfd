function [Ap,Ae,Aw,An,As,At,Ab,bp] = make_dense_coeffs(m,n,l,...
    dx,dy,dz,k_c,bc_mat,Tn)

% Create Coefficient Matrices
An = zeros(m,n,l);
As = zeros(m,n,l);
Ae = zeros(m,n,l);
Aw = zeros(m,n,l);
At = zeros(m,n,l);
Ab = zeros(m,n,l);

Sc = zeros(m,n,l);
Sp = zeros(m,n,l);

Apn = zeros(m,n,l);
bp = zeros(m,n,l);
Ap = zeros(m,n,l);

% Assign Linearized Source Terms to Sides
for i = 2:1:m-1
    for j = 2:1:n-1
        % Bottom
        Sc(i,j,1) = bc_mat(6,1)*(dx*dy/dz)/dz-...
            bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
        if bc_mat(6,3) ~= 0
            Sp(i,j,1) = -(2*k_c*(dx*dy/dz))/dz;
        end
        
        % Top
        Sc(i,j,l) = bc_mat(5,1)*(dx*dy/dz)/dz-...
            bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
        if bc_mat(5,3) ~= 0
            Sp(i,j,l) = -(2*k_c*(dx*dy/dz))/dz;
        end
    end
end

for i = 2:1:m-1
    for k = 2:1:n-1
        % South
        Sc(i,1,k) = bc_mat(4,1)*(dx*dz/dy)/dy-...
            bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy;
        if bc_mat(4,3) ~= 0
            Sp(i,1,k) = -(2*k_c*(dx*dz/dy))/dy;
        end
        
        % North
        Sc(i,n,k) = bc_mat(3,1)*(dx*dz/dy)/dy-...
            bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy;
        if bc_mat(3,3) ~= 0
            Sp(i,n,k) = -(2*k_c*(dx*dz/dy))/dy;
        end
    end
end

for k = 2:1:l-1
    for j = 2:1:n-1
        % West
        Sc(1,j,k) = bc_mat(2,1)*(dz*dy/dx)/dx-...
            bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx;
        if bc_mat(2,3) ~= 0
            Sp(1,j,k) = -(2*k_c*(dz*dy/dx))/dx;
        end
        
        % East
        Sc(m,j,k) = bc_mat(1,1)*(dz*dy/dx)/dx-...
            bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx;
        if bc_mat(1,3) ~= 0
            Sp(m,j,k) = -(2*k_c*(dz*dy/dx))/dx;
        end
    end
end

% Add Linearized Source Terms to Edges
% West South
Sc(1,1,:) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy;
if bc_mat(2,3) ~= 0 && bc_mat(4,3) ~= 0
    Sp(1,1,:) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(2,3) ~= 0
    Sp(1,1,:) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(4,3) ~= 0
    Sp(1,1,:) = -(2*k_c*(dx*dz/dy))/dy;
end

% West North
Sc(1,n,:) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy;
if bc_mat(2,3) ~= 0 && bc_mat(3,3) ~= 0
    Sp(1,n,:) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(2,3) ~= 0
    Sp(1,n,:) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(3,3) ~= 0
    Sp(1,n,:) = -(2*k_c*(dx*dz/dy))/dy;
end

% East South
Sc(m,1,:) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy;
if bc_mat(1,3) ~= 0 && bc_mat(4,3) ~= 0
    Sp(m,1,:) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(1,3) ~= 0
    Sp(m,1,:) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(4,3) ~= 0
    Sp(m,1,:) = -(2*k_c*(dx*dz/dy))/dy;
end

% East North
Sc(m,n,:) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy;
if bc_mat(1,3) ~= 0 && bc_mat(3,3) ~= 0
    Sp(m,n,:) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(1,3) ~= 0
    Sp(m,n,:) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(3,3) ~= 0
    Sp(m,n,:) = -(2*k_c*(dx*dz/dy))/dy;
end

% West Bottom
Sc(1,:,1) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(2,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(1,:,1) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0
    Sp(1,:,1) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(6,3) ~= 0
    Sp(1,:,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% West Top
Sc(1,:,l) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(2,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(1,:,l) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0
    Sp(1,:,l) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(5,3) ~= 0
    Sp(1,:,l) = -(2*k_c*(dx*dy/dz))/dz;
end

% East Bottom
Sc(m,:,1) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(1,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,:,1) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0
    Sp(m,:,1) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(6,3) ~= 0
    Sp(m,:,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% East Top
Sc(m,:,l) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(1,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(m,:,l) = -(2*k_c*(dz*dy/dx))/dx-...
        (2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0
    Sp(m,:,l) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(5,3) ~= 0
    Sp(m,:,l) = -(2*k_c*(dx*dy/dz))/dz;
end

% South Bottom
Sc(:,1,1) = bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(4,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(:,1,1) = -(2*k_c*(dx*dz/dy))/dy-...
        (2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(4,3) ~= 0
    Sp(:,1,1) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(6,3) ~= 0
    Sp(:,1,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% South Top
Sc(:,1,l) = bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(4,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(:,1,l) = -(2*k_c*(dx*dz/dy))/dy-...
        (2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(4,3) ~= 0
    Sp(:,1,l) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(5,3) ~= 0
    Sp(:,1,l) = -(2*k_c*(dx*dy/dz))/dz;
end

% North Bottom
Sc(:,n,1) = bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(3,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(:,n,1) = -(2*k_c*(dx*dz/dy))/dy-...
        (2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(3,3) ~= 0
    Sp(:,n,1) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(6,3) ~= 0
    Sp(:,n,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% North Top
Sc(:,n,l) = bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(3,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(:,n,l) = -(2*k_c*(dx*dz/dy))/dy-...
        (2*k_c*(dx*dy/dz))/dz;
end

% Add Linearized Source Terms to Corners
% West South Bottom
Sc(1,1,1) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(2,3) ~= 0 && bc_mat(4,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(1,1,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(4,3) ~= 0
    Sp(1,1,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(4,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(1,1,1) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(1,1,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0
    Sp(1,1,1) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(4,3) ~= 0
    Sp(1,1,1) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(6,3) ~= 0
    Sp(1,1,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% West South Top
Sc(1,1,l) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(2,3) ~= 0 && bc_mat(4,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(1,1,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(4,3) ~= 0
    Sp(1,1,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(4,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(1,1,l) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(1,1,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0
    Sp(1,1,l) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(4,3) ~= 0
    Sp(1,1,l) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(5,3) ~= 0
    Sp(1,1,l) = -(2*k_c*(dx*dy/dz))/dz;
end

% West North Bottom
Sc(1,n,1) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(2,3) ~= 0 && bc_mat(3,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(1,n,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(3,3) ~= 0
    Sp(1,n,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(3,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(1,n,1) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(1,n,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0
    Sp(1,n,1) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(3,3) ~= 0
    Sp(1,n,1) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(6,3) ~= 0
    Sp(1,n,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% West North Top
Sc(1,n,l) = bc_mat(2,1)*(dz*dy/dx)/dx-bc_mat(2,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(2,3) ~= 0 && bc_mat(3,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(1,n,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(3,3) ~= 0
    Sp(1,n,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(3,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(1,n,l) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(1,n,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(2,3) ~= 0
    Sp(1,n,l) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(3,3) ~= 0
    Sp(1,n,l) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(5,3) ~= 0
    Sp(1,n,l) = -(2*k_c*(dx*dy/dz))/dz;
end

% East North Bottom
Sc(m,n,1) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(1,3) ~= 0 && bc_mat(3,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,n,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(3,3) ~= 0
    Sp(m,n,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(3,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,n,1) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,n,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0
    Sp(m,n,1) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(3,3) ~= 0
    Sp(m,n,1) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(6,3) ~= 0
    Sp(m,n,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% East North Top
Sc(m,n,l) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(3,1)*(dx*dz/dy)/dy-bc_mat(3,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(1,3) ~= 0 && bc_mat(3,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(m,n,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(3,3) ~= 0
    Sp(m,n,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(3,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(m,n,l) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,n,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0
    Sp(m,n,l) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(3,3) ~= 0
    Sp(m,n,l) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(5,3) ~= 0
    Sp(m,n,l) = -(2*k_c*(dx*dy/dz))/dz;
end

% East South Bottom
Sc(m,1,1) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(6,1)*(dx*dy/dz)/dz-bc_mat(6,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(1,3) ~= 0 && bc_mat(4,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,1,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(4,3) ~= 0
    Sp(m,1,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(4,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,1,1) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(6,3) ~= 0
    Sp(m,1,1) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0
    Sp(m,1,1) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(4,3) ~= 0
    Sp(m,1,1) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(6,3) ~= 0
    Sp(m,1,1) = -(2*k_c*(dx*dy/dz))/dz;
end

% East South Top
Sc(m,1,l) = bc_mat(1,1)*(dz*dy/dx)/dx-bc_mat(1,3)*(2*k_c*(dz*dy/dx))/dx+...
    bc_mat(4,1)*(dx*dz/dy)/dy-bc_mat(4,3)*(2*k_c*(dx*dz/dy))/dy+...
    bc_mat(5,1)*(dx*dy/dz)/dz-bc_mat(5,3)*(2*k_c*(dx*dy/dz))/dz;
if bc_mat(1,3) ~= 0 && bc_mat(4,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(m,1,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(4,3) ~= 0
    Sp(m,1,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(4,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(m,1,l) = -(2*k_c*(dx*dz/dy))/dy...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0 && bc_mat(5,3) ~= 0
    Sp(m,1,l) = -(2*k_c*(dz*dy/dx))/dx...
        -(2*k_c*(dx*dy/dz))/dz;
elseif bc_mat(1,3) ~= 0
    Sp(m,1,l) = -(2*k_c*(dz*dy/dx))/dx;
elseif bc_mat(4,3) ~= 0
    Sp(m,1,l) = -(2*k_c*(dx*dz/dy))/dy;
elseif bc_mat(5,3) ~= 0
    Sp(m,1,l) = -(2*k_c*(dx*dy/dz))/dz;
end

% Assign Values to Coefficient Matrices (Interior Nodes)
for i = 2:1:m-1
    for j = 2:1:n-1
        for k = 2:1:l-1
            An(i,j,k) = k_c*(dx*dz/dy)/dy;
            As(i,j,k) = k_c*(dx*dz/dy)/dy;
            Ae(i,j,k) = k_c*(dy*dz/dx)/dx;
            Aw(i,j,k) = k_c*(dy*dz/dx)/dx;
            At(i,j,k) = k_c*(dx*dy/dz)/dz;
            Ab(i,j,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(i,j,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,j,k) = Sc(i,j,k)*dx*dy*dz+Apn(i,j,k)*Tn(i,j,k);
            Ap(i,j,k) = An(i,j,k)+As(i,j,k)+Aw(i,j,k)+Ae(i,j,k)+...
                At(i,j,k)+Ab(i,j,k)+Apn(i,j,k)-Sp(i,j,k)*dx*dy*dz;
        end
    end
end

% Assign Values to Coefficient Matrices (Side Nodes)
% West Side
for j = 2:1:n-1
    for k = 2:1:l-1
            An(1,j,k) = k_c*(dx*dz/dy)/dy;
            As(1,j,k) = k_c*(dx*dz/dy)/dy;
            Ae(1,j,k) = k_c*(dy*dz/dx)/dx;
            Aw(1,j,k) = 0;
            At(1,j,k) = k_c*(dx*dy/dz)/dz;
            Ab(1,j,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(1,j,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(1,j,k) = Sc(1,j,k)*dx*dy*dz+Apn(1,j,k)*Tn(1,j,k);
            Ap(1,j,k) = An(1,j,k)+As(1,j,k)+Aw(1,j,k)+Ae(1,j,k)+...
                At(1,j,k)+Ab(1,j,k)+Apn(1,j,k)-Sp(1,j,k)*dx*dy*dz;
    end
end

% East Side
for j = 2:1:n-1
    for k = 2:1:l-1
            An(m,j,k) = k_c*(dx*dz/dy)/dy;
            As(m,j,k) = k_c*(dx*dz/dy)/dy;
            Ae(m,j,k) = 0;
            Aw(m,j,k) = k_c*(dy*dz/dx)/dx;
            At(m,j,k) = k_c*(dx*dy/dz)/dz;
            Ab(m,j,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(m,j,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(m,j,k) = Sc(m,j,k)*dx*dy*dz+Apn(m,j,k)*Tn(m,j,k);
            Ap(m,j,k) = An(m,j,k)+As(m,j,k)+Aw(m,j,k)+Ae(m,j,k)+...
                At(m,j,k)+Ab(m,j,k)+Apn(m,j,k)-Sp(m,j,k)*dx*dy*dz;
    end
end

% South Side
for i = 2:1:m-1
    for k = 2:1:l-1
            An(i,1,k) = k_c*(dx*dz/dy)/dy;
            As(i,1,k) = 0;
            Ae(i,1,k) = k_c*(dy*dz/dx)/dx;
            Aw(i,1,k) = k_c*(dy*dz/dx)/dx;
            At(i,1,k) = k_c*(dx*dy/dz)/dz;
            Ab(i,1,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(i,1,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,1,k) = Sc(i,1,k)*dx*dy*dz+Apn(i,1,k)*Tn(i,1,k);
            Ap(i,1,k) = An(i,1,k)+As(i,1,k)+Aw(i,1,k)+Ae(i,1,k)+...
                At(i,1,k)+Ab(i,1,k)+Apn(i,1,k)-Sp(i,1,k)*dx*dy*dz;
    end
end

% North Side
for i = 2:1:m-1
    for k = 2:1:l-1
            An(i,n,k) = 0;
            As(i,n,k) = k_c*(dx*dz/dy)/dy;
            Ae(i,n,k) = k_c*(dy*dz/dx)/dx;
            Aw(i,n,k) = k_c*(dy*dz/dx)/dx;
            At(i,n,k) = k_c*(dx*dy/dz)/dz;
            Ab(i,n,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(i,n,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,n,k) = Sc(i,n,k)*dx*dy*dz+Apn(i,n,k)*Tn(i,n,k);
            Ap(i,n,k) = An(i,n,k)+As(i,n,k)+Aw(i,n,k)+Ae(i,n,k)+...
                At(i,n,k)+Ab(i,n,k)+Apn(i,n,k)-Sp(i,n,k)*dx*dy*dz;
    end
end

% Bottom Side
for i = 2:1:m-1
    for j = 2:1:n-1
            An(i,j,1) = k_c*(dx*dz/dy)/dy;
            As(i,j,1) = k_c*(dx*dz/dy)/dy;
            Ae(i,j,1) = k_c*(dy*dz/dx)/dx;
            Aw(i,j,1) = k_c*(dy*dz/dx)/dx;
            At(i,j,1) = k_c*(dx*dy/dz)/dz;
            Ab(i,j,1) = 0;
            
            Apn(i,j,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,j,1) = Sc(i,j,1)*dx*dy*dz+Apn(i,j,1)*Tn(i,j,1);
            Ap(i,j,1) = An(i,j,1)+As(i,j,1)+Aw(i,j,1)+Ae(i,j,1)+...
                At(i,j,1)+Ab(i,j,1)+Apn(i,j,1)-Sp(i,j,1)*dx*dy*dz;
    end
end

% Top Side
for i = 2:1:m-1
    for j = 2:1:n-1
            An(i,j,l) = k_c*(dx*dz/dy)/dy;
            As(i,j,l) = k_c*(dx*dz/dy)/dy;
            Ae(i,j,l) = k_c*(dy*dz/dx)/dx;
            Aw(i,j,l) = k_c*(dy*dz/dx)/dx;
            At(i,j,l) = 0;
            Ab(i,j,l) = k_c*(dx*dy/dz)/dz;
            
            Apn(i,j,l) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,j,l) = Sc(i,j,l)*dx*dy*dz+Apn(i,j,l)*Tn(i,j,l);
            Ap(i,j,l) = An(i,j,l)+As(i,j,l)+Aw(i,j,l)+Ae(i,j,l)+...
                At(i,j,l)+Ab(i,j,l)+Apn(i,j,l)-Sp(i,j,l)*dx*dy*dz;
    end
end

% Assign Values to Coefficient Matrices (Edge Nodes)
% West South
for k = 2:1:l-1
            An(1,1,k) = k_c*(dx*dz/dy)/dy;
            As(1,1,k) = 0;
            Ae(1,1,k) = k_c*(dy*dz/dx)/dx;
            Aw(1,1,k) = 0;
            At(1,1,k) = k_c*(dx*dy/dz)/dz;
            Ab(1,1,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(1,1,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(1,1,k) = Sc(1,1,k)*dx*dy*dz+Apn(1,1,k)*Tn(1,1,k);
            Ap(1,1,k) = An(1,1,k)+As(1,1,k)+Aw(1,1,k)+Ae(1,1,k)+...
                At(1,1,k)+Ab(1,1,k)+Apn(1,1,k)-Sp(1,1,k)*dx*dy*dz;
end

% West North
for k = 2:1:l-1
            An(1,n,k) = 0;
            As(1,n,k) = k_c*(dx*dz/dy)/dy;
            Ae(1,n,k) = k_c*(dy*dz/dx)/dx;
            Aw(1,n,k) = 0;
            At(1,n,k) = k_c*(dx*dy/dz)/dz;
            Ab(1,n,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(1,n,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(1,n,k) = Sc(1,n,k)*dx*dy*dz+Apn(1,n,k)*Tn(1,n,k);
            Ap(1,n,k) = An(1,n,k)+As(1,n,k)+Aw(1,n,k)+Ae(1,n,k)+...
                At(1,n,k)+Ab(1,n,k)+Apn(1,n,k)-Sp(1,n,k)*dx*dy*dz;
end

% North East
for k = 2:1:l-1
            An(m,n,k) = 0;
            As(m,n,k) = k_c*(dx*dz/dy)/dy;
            Ae(m,n,k) = 0;
            Aw(m,n,k) = k_c*(dy*dz/dx)/dx;
            At(m,n,k) = k_c*(dx*dy/dz)/dz;
            Ab(m,n,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(m,n,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(m,n,k) = Sc(m,n,k)*dx*dy*dz+Apn(m,n,k)*Tn(m,n,k);
            Ap(m,n,k) = An(m,n,k)+As(m,n,k)+Aw(m,n,k)+Ae(m,n,k)+...
                At(m,n,k)+Ab(m,n,k)+Apn(m,n,k)-Sp(m,n,k)*dx*dy*dz;
end

% East South
for k = 2:1:l-1
            An(m,1,k) = k_c*(dx*dz/dy)/dy;
            As(m,1,k) = 0;
            Ae(m,1,k) = 0;
            Aw(m,1,k) = k_c*(dy*dz/dx)/dx;
            At(m,1,k) = k_c*(dx*dy/dz)/dz;
            Ab(m,1,k) = k_c*(dx*dy/dz)/dz;
            
            Apn(m,1,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(m,1,k) = Sc(m,1,k)*dx*dy*dz+Apn(m,1,k)*Tn(m,1,k);
            Ap(m,1,k) = An(m,1,k)+As(m,1,k)+Aw(m,1,k)+Ae(m,1,k)+...
                At(m,1,k)+Ab(m,1,k)+Apn(m,1,k)-Sp(m,1,k)*dx*dy*dz;
end

% West Bottom
for j = 2:1:n-1
            An(1,j,1) = k_c*(dx*dz/dy)/dy;
            As(1,j,1) = k_c*(dx*dz/dy)/dy;
            Ae(1,j,1) = k_c*(dy*dz/dx)/dx;
            Aw(1,j,1) = 0;
            At(1,j,1) = k_c*(dx*dy/dz)/dz;
            Ab(1,j,1) = 0;
            
            Apn(1,j,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(1,j,1) = Sc(1,j,1)*dx*dy*dz+Apn(1,j,1)*Tn(1,j,1);
            Ap(1,j,1) = An(1,j,1)+As(1,j,1)+Aw(1,j,1)+Ae(1,j,1)+...
                At(1,j,1)+Ab(1,j,1)+Apn(1,j,1)-Sp(1,j,1)*dx*dy*dz;
end

% North Bottom
for i = 2:1:m-1
            An(i,n,1) = 0;
            As(i,n,1) = k_c*(dx*dz/dy)/dy;
            Ae(i,n,1) = k_c*(dy*dz/dx)/dx;
            Aw(i,n,1) = k_c*(dy*dz/dx)/dx;
            At(i,n,1) = k_c*(dx*dy/dz)/dz;
            Ab(i,n,1) = 0;
            
            Apn(i,n,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,n,1) = Sc(i,n,1)*dx*dy*dz+Apn(i,n,1)*Tn(i,n,1);
            Ap(i,n,1) = An(i,n,1)+As(i,n,1)+Aw(i,n,1)+Ae(i,n,1)+...
                At(i,n,1)+Ab(i,n,1)+Apn(i,n,1)-Sp(i,n,1)*dx*dy*dz;
end

% East Bottom
for j = 2:1:n-1
            An(m,j,1) = k_c*(dx*dz/dy)/dy;
            As(m,j,1) = k_c*(dx*dz/dy)/dy;
            Ae(m,j,1) = 0;
            Aw(m,j,1) = k_c*(dy*dz/dx)/dx;
            At(m,j,1) = k_c*(dx*dy/dz)/dz;
            Ab(m,j,1) = 0;
            
            Apn(m,j,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(m,j,1) = Sc(m,j,1)*dx*dy*dz+Apn(m,j,1)*Tn(m,j,1);
            Ap(m,j,1) = An(m,j,1)+As(m,j,1)+Aw(m,j,1)+Ae(m,j,1)+...
                At(m,j,1)+Ab(m,j,1)+Apn(m,j,1)-Sp(m,j,1)*dx*dy*dz;
end

% South Bottom
for i = 2:1:m-1
            An(i,1,1) = k_c*(dx*dz/dy)/dy;
            As(i,1,1) = 0;
            Ae(i,1,1) = k_c*(dy*dz/dx)/dx;
            Aw(i,1,1) = k_c*(dy*dz/dx)/dx;
            At(i,1,1) = k_c*(dx*dy/dz)/dz;
            Ab(i,1,1) = 0;
            
            Apn(i,1,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,1,1) = Sc(i,1,1)*dx*dy*dz+Apn(i,1,1)*Tn(i,1,1);
            Ap(i,1,1) = An(i,1,1)+As(i,1,1)+Aw(i,1,1)+Ae(i,1,1)+...
                At(i,1,1)+Ab(i,1,1)+Apn(i,1,1)-Sp(i,1,1)*dx*dy*dz;
end

% West Top
for j = 2:1:n-1
            An(1,j,l) = k_c*(dx*dz/dy)/dy;
            As(1,j,l) = k_c*(dx*dz/dy)/dy;
            Ae(1,j,l) = k_c*(dy*dz/dx)/dx;
            Aw(1,j,l) = 0;
            At(1,j,l) = 0;
            Ab(1,j,l) = k_c*(dx*dy/dz)/dz;
            
            Apn(1,j,l) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(1,j,l) = Sc(1,j,l)*dx*dy*dz+Apn(1,j,l)*Tn(1,j,l);
            Ap(1,j,l) = An(1,j,l)+As(1,j,l)+Aw(1,j,l)+Ae(1,j,l)+...
                At(1,j,l)+Ab(1,j,l)+Apn(1,j,l)-Sp(1,j,l)*dx*dy*dz;
end

% North Top
for i = 2:1:m-1
            An(i,n,l) = 0;
            As(i,n,l) = k_c*(dx*dz/dy)/dy;
            Ae(i,n,l) = k_c*(dy*dz/dx)/dx;
            Aw(i,n,l) = k_c*(dy*dz/dx)/dx;
            At(i,n,l) = 0;
            Ab(i,n,l) = k_c*(dx*dy/dz)/dz;
            
            Apn(i,n,l) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,n,l) = Sc(i,n,l)*dx*dy*dz+Apn(i,n,l)*Tn(i,n,l);
            Ap(i,n,l) = An(i,n,l)+As(i,n,l)+Aw(i,n,l)+Ae(i,n,l)+...
                At(i,n,l)+Ab(i,n,l)+Apn(i,n,l)-Sp(i,n,l)*dx*dy*dz;
end

% East Top
for j = 2:1:n-1
            An(m,j,l) = k_c*(dx*dz/dy)/dy;
            As(m,j,l) = k_c*(dx*dz/dy)/dy;
            Ae(m,j,l) = 0;
            Aw(m,j,l) = k_c*(dy*dz/dx)/dx;
            At(m,j,l) = 0;
            Ab(m,j,l) = k_c*(dx*dy/dz)/dz;
            
            Apn(m,j,l) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(m,j,l) = Sc(m,j,l)*dx*dy*dz+Apn(m,j,l)*Tn(m,j,l);
            Ap(m,j,l) = An(m,j,l)+As(m,j,l)+Aw(m,j,l)+Ae(m,j,l)+...
                At(m,j,l)+Ab(m,j,l)+Apn(m,j,l)-Sp(m,j,l)*dx*dy*dz;
end

% South Top
for i = 2:1:m-1
            An(i,1,l) = k_c*(dx*dz/dy)/dy;
            As(i,1,l) = 0;
            Ae(i,1,l) = k_c*(dy*dz/dx)/dx;
            Aw(i,1,l) = k_c*(dy*dz/dx)/dx;
            At(i,1,l) = 0;
            Ab(i,1,l) = k_c*(dx*dy/dz)/dz;
            
            Apn(i,1,l) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
            bp(i,1,l) = Sc(i,1,l)*dx*dy*dz+Apn(i,1,l)*Tn(i,1,l);
            Ap(i,1,l) = An(i,1,l)+As(i,1,l)+Aw(i,1,l)+Ae(i,1,l)+...
                At(i,1,l)+Ab(i,1,l)+Apn(i,1,l)-Sp(i,1,l)*dx*dy*dz;
end

% Assign Values to Coefficient Matrices (Corner Nodes)
% West South Bottom
An(1,1,1) = k_c*(dx*dz/dy)/dy;
As(1,1,1) = 0;
Ae(1,1,1) = k_c*(dy*dz/dx)/dx;
Aw(1,1,1) = 0;
At(1,1,1) = k_c*(dx*dy/dz)/dz;
Ab(1,1,1) = 0;

Apn(1,1,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(1,1,1) = Sc(1,1,1)*dx*dy*dz+Apn(1,1,1)*Tn(1,1,1);
Ap(1,1,1) = An(1,1,1)+As(1,1,1)+Aw(1,1,1)+Ae(1,1,1)+...
    At(1,1,1)+Ab(1,1,1)+Apn(1,1,1)-Sp(1,1,1)*dx*dy*dz;

% West North Bottom
An(1,n,1) = 0;
As(1,n,1) = k_c*(dx*dz/dy)/dy;
Ae(1,n,1) = k_c*(dy*dz/dx)/dx;
Aw(1,n,1) = 0;
At(1,n,1) = k_c*(dx*dy/dz)/dz;
Ab(1,n,1) = 0;

Apn(1,n,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(1,n,1) = Sc(1,n,1)*dx*dy*dz+Apn(1,n,1)*Tn(1,n,1);
Ap(1,n,1) = An(1,n,1)+As(1,n,1)+Aw(1,n,1)+Ae(1,n,1)+...
    At(1,n,1)+Ab(1,n,1)+Apn(1,n,1)-Sp(1,n,1)*dx*dy*dz;

% North East Bottom
An(m,n,1) = 0;
As(m,n,1) = k_c*(dx*dz/dy)/dy;
Ae(m,n,1) = 0;
Aw(m,n,1) = k_c*(dy*dz/dx)/dx;
At(m,n,1) = k_c*(dx*dy/dz)/dz;
Ab(m,n,1) = 0;

Apn(m,n,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(m,n,1) = Sc(m,n,1)*dx*dy*dz+Apn(m,n,1)*Tn(m,n,1);
Ap(m,n,1) = An(m,n,1)+As(m,n,1)+Aw(m,n,1)+Ae(m,n,1)+...
    At(m,n,1)+Ab(m,n,1)+Apn(m,n,1)-Sp(m,n,1)*dx*dy*dz;

% East South Bottom
An(m,1,1) = k_c*(dx*dz/dy)/dy;
As(m,1,1) = 0;
Ae(m,1,1) = 0;
Aw(m,1,1) = k_c*(dy*dz/dx)/dx;
At(m,1,1) = k_c*(dx*dy/dz)/dz;
Ab(m,1,1) = 0;

Apn(m,1,1) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(m,1,1) = Sc(m,1,1)*dx*dy*dz+Apn(m,1,1)*Tn(m,1,1);
Ap(m,1,1) = An(m,1,1)+As(m,1,1)+Aw(m,1,1)+Ae(m,1,1)+...
    At(m,1,1)+Ab(m,1,1)+Apn(m,1,1)-Sp(m,1,1)*dx*dy*dz;

% West South Top
An(1,1,k) = k_c*(dx*dz/dy)/dy;
As(1,1,k) = 0;
Ae(1,1,k) = k_c*(dy*dz/dx)/dx;
Aw(1,1,k) = 0;
At(1,1,k) = 0;
Ab(1,1,k) = k_c*(dx*dy/dz)/dz;

Apn(1,1,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(1,1,k) = Sc(1,1,k)*dx*dy*dz+Apn(1,1,k)*Tn(1,1,k);
Ap(1,1,k) = An(1,1,k)+As(1,1,k)+Aw(1,1,k)+Ae(1,1,k)+...
    At(1,1,k)+Ab(1,1,k)+Apn(1,1,k)-Sp(1,1,k)*dx*dy*dz;

% West North Top
An(1,n,k) = 0;
As(1,n,k) = k_c*(dx*dz/dy)/dy;
Ae(1,n,k) = k_c*(dy*dz/dx)/dx;
Aw(1,n,k) = 0;
At(1,n,k) = 0;
Ab(1,n,k) = k_c*(dx*dy/dz)/dz;

Apn(1,n,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(1,n,k) = Sc(1,n,k)*dx*dy*dz+Apn(1,n,k)*Tn(1,n,k);
Ap(1,n,k) = An(1,n,k)+As(1,n,k)+Aw(1,n,k)+Ae(1,n,k)+...
    At(1,n,k)+Ab(1,n,k)+Apn(1,n,k)-Sp(1,n,k)*dx*dy*dz;

% North East Top
An(m,n,k) = 0;
As(m,n,k) = k_c*(dx*dz/dy)/dy;
Ae(m,n,k) = 0;
Aw(m,n,k) = k_c*(dy*dz/dx)/dx;
At(m,n,k) = 0;
Ab(m,n,k) = k_c*(dx*dy/dz)/dz;

Apn(m,n,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(m,n,k) = Sc(m,n,k)*dx*dy*dz+Apn(m,n,k)*Tn(m,n,k);
Ap(m,n,k) = An(m,n,k)+As(m,n,k)+Aw(m,n,k)+Ae(m,n,k)+...
    At(m,n,k)+Ab(m,n,k)+Apn(m,n,k)-Sp(m,n,k)*dx*dy*dz;

% East South Top
An(m,1,k) = k_c*(dx*dz/dy)/dy;
As(m,1,k) = 0;
Ae(m,1,k) = 0;
Aw(m,1,k) = k_c*(dy*dz/dx)/dx;
At(m,1,k) = 0;
Ab(m,1,k) = k_c*(dx*dy/dz)/dz;

Apn(m,1,k) = 0; %rho*c_sp*dx*dy*dz/dt; Using SS
bp(m,1,k) = Sc(m,1,k)*dx*dy*dz+Apn(m,1,k)*Tn(m,1,k);
Ap(m,1,k) = An(m,1,k)+As(m,1,k)+Aw(m,1,k)+Ae(m,1,k)+...
    At(m,1,k)+Ab(m,1,k)+Apn(m,1,k)-Sp(m,1,k)*dx*dy*dz;

return
