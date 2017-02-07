% Test Different Solution Algorithms

A = sparse(row,col,val);

[L U] = ilu(A);

fprintf('Biconjugate gradients stabilized method\n');
tic
[x, bg_f, bg_rr, bg_itr, bg_resvec] = bicgstab(A,d',1e-6,1000);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bg_rr,bg_itr);
fprintf(text);

fprintf('Biconjugate gradients stabilized method with LU Preconditioning\n');
tic
[x, bgp_f, bgp_rr, bgp_itr, bgp_resvec] = bicgstab(A,d',1e-6,1000,L,U);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bg_rr,bg_itr);
fprintf(text);

fprintf('Generalized minimum residual method (with restarts)\n');
tic
[x, gm_f, gm_rr, gm_itr, gm_resvec] = gmres(A,d',100,1e-6);
toc
output_text = 'Res: %1.3e\nItr: [%d %d]\n\n';
text = sprintf(output_text,gm_rr, gm_itr);
fprintf(text);

fprintf('Generalized minimum residual method (with restarts) with PC\n');
tic
[x, gm_f, gm_rr, gm_itr, gmp_resvec] = gmres(A,d',100,1e-6,100,L,U);
toc
output_text = 'Res: %1.3e\nItr: [%d %d]\n\n';
text = sprintf(output_text,gm_rr, gm_itr);
fprintf(text);


norm_A = norm(A,'fro');

bg_t = 0:1:(length(bg_resvec)-1);
bgp_t = 0:1:(length(bgp_resvec)-1);
gm_t = 0:1:(length(gm_resvec)-1);
gmp_t = 0:1:(length(gmp_resvec)-1);

semilogy(bg_t,bg_resvec/norm_A)
hold on
semilogy(bgp_t,bgp_resvec/norm_A)
semilogy(gm_t,gm_resvec/norm_A)
semilogy(gmp_t,gmp_resvec/norm_A)
plot([0 length(bg_resvec)],[1e-6 1e-6])

grid minor
xlabel 'Iterations'
ylabel 'Relative Residual'

titleSpec = 'Convergence - 3D Grid (%dx%dx%d)';
chart_title = sprintf(titleSpec,m,n,l);
title(chart_title)
legend 'BiCGStab' 'LU-BiCGStab' 'GMRES' 'LU-GMRES'
