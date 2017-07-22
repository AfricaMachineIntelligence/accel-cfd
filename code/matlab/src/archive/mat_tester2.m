% Test Different Solution Algorithms

% Convert Output CSC into a Spare Matrix
fprintf('Matrix Conversion\n');
tic
A = sparse(row,col,val);
toc
fprintf('Matrix Conversion Done\n\n');

% Clear Data
clear row col val

% Start L, U factorization for Preconditioning
d_val = spdiags(A,0);
d_row = 1:1:length(d_val);
d_mat = sparse(d_row,d_row,d_val);

[L U] = ilu(A);

% Set Tolerance for Iterative Solvers
e_tol = 1e-6;

% Run Biconjugate Gradient Solvers
fprintf('Biconjugate gradients\n');
tic
[x, bc_f, bc_rr, bc_itr, bc_resvec] = bicg(A,d',e_tol,1000);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bc_rr,bc_itr);
fprintf(text);

fprintf('Biconjugate gradients with PC\n');
tic
[x, bcp_f, bcp_rr, bcp_itr, bcp_resvec] = bicg(A,d',e_tol,1000,L,U);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bcp_rr,bcp_itr);
fprintf(text);

% Run Biconjugate Gradient Stabilized Solvers
fprintf('Biconjugate gradients stabilized method\n');
tic
[x, bg_f, bg_rr, bg_itr, bg_resvec] = bicgstab(A,d',e_tol,1000);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bg_rr,bg_itr);
fprintf(text);

fprintf('Biconjugate gradients stabilized method with PC\n');
tic
[x, bgp_f, bgp_rr, bgp_itr, bgp_resvec] = bicgstab(A,d',e_tol,1000,L,U);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bgp_rr,bgp_itr);
fprintf(text);

% Run GMRES Solvers
fprintf('Generalized minimum residual method (with restarts)\n');
tic
[x, gm_f, gm_rr, gm_itr, gm_resvec] = gmres(A,d',100,e_tol,100);
toc
output_text = 'Res: %1.3e\nItr: [%d %d]\n\n';
text = sprintf(output_text,gm_rr, gm_itr);
fprintf(text);

fprintf('Generalized minimum residual method (with restarts) with PC\n');
tic
[x, gmp_f, gmp_rr, gmp_itr, gmp_resvec] = gmres(A,d',100,e_tol,100,L,U);
toc
output_text = 'Res: %1.3e\nItr: [%d %d]\n\n';
text = sprintf(output_text,gmp_rr, gmp_itr);
fprintf(text);

% Create Arrays for Plotting Residule Vectors
bc_t = 0:1:(length(bc_resvec)-1);
bcp_t = 0:1:(length(bcp_resvec)-1);
bg_t = 0:1:(length(bg_resvec)-1);
bgp_t = 0:1:(length(bgp_resvec)-1);
gm_t = 0:1:(length(gm_resvec)-1);
gmp_t = 0:1:(length(gmp_resvec)-1);

% Calculate Norm of A
nA = norm(A,'fro');

% Plot Residule vs. Iteration
semilogy(bc_t,bc_resvec/nA)
hold on
semilogy(bcp_t,bcp_resvec/nA)
semilogy(bg_t,bg_resvec/nA)
semilogy(bgp_t,bgp_resvec/nA)
semilogy(gm_t,gm_resvec/nA)
semilogy(gmp_t,gmp_resvec/nA)
plot([0 length(bc_resvec)],[1e-6 1e-6],'k')

grid minor
xlabel 'Iteration'
ylabel 'Relative Residual'
legend 'BiCG' 'BiCG w/ LU' 'BiCGStab' 'BiCGStab w/ LU' 'GMRES' 'GMRES w/ LU' 
