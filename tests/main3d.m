% Conversion Algorithm for Mukherjee TDMA

% Clear all previous data
clear variables
close all
clc

% Set L, W, H
grid_L = .1;         % Actual dimension of grid x
grid_H = .1;         % Actual dimension of grid y
grid_W = .1;         % Actual dimension of grid z

% Set Material Conditions
k_c = 385;            % Termal Conductivity

% Set Geometry Conditions
m = 30;     % x-dim
n = 30;     % y-dim
l = 30;     % z-dim

% Create Temperature Matrix
T = zeros(m,n,l);   % Current Time Step 
Tn = zeros(m,n,l);  % New Time Step // Off for Steady State

% Non-dimensionalize
dx = grid_L/m;      % Actual dimension of grid x
dy = grid_H/n;      % Actual dimension of grid y
dz = grid_W/l;      % Actual dimension of grid z

% Define Boundary Conditions
% North BC
bc_qc_n = 0;
bc_qh_n = 0;
bc_T_n = 0;

% South BC
bc_qc_s = 0;
bc_qh_s = 0;
bc_T_s = 0;

% East BC
bc_qc_e = 0;
bc_qh_e = 0;
bc_T_e = 0;

% West BC
bc_qc_w = 500*1e7;
bc_qh_w = 0;
bc_T_w = 0;

% Top BC
bc_qc_t = 0;
bc_qh_t = 0;
bc_T_t = 100;

% Bottom BC
bc_qc_b = 0;
bc_qh_b = 0;
bc_T_b = 0;

% Create Matrix of Boundary Conditions
bc_mat = [bc_qc_e,bc_qh_e,bc_T_e;
    bc_qc_w,bc_qh_w,bc_T_w;
    bc_qc_n,bc_qh_n,bc_T_n;
    bc_qc_s,bc_qh_s,bc_T_s;
    bc_qc_t,bc_qh_t,bc_T_t;
    bc_qc_b,bc_qh_b,bc_T_b];

% Generate Coefficient Matrices
[Ap,Ae,Aw,An,As,At,Ab,bp] = make_dense_coeffs(m,n,l,...
    dx,dy,dz,k_c,bc_mat,Tn);

% Generate Sparse Matrix
[A,d] = make_sparse_matrix(Ap,Ae,Aw,An,As,At,Ab,bp);

x = d/A;

for k = 1:1:l
    for i = 1:1:m
        for j = 1:1:n
            T(i,j,k) = x(j+(i-1)*n+(k-1)*n*m);
        end
    end
end

plane1 = squeeze(T(:,:,15));
contourf(plane1)
