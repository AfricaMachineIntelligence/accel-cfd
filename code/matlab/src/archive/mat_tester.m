% Test Different Solution Algorithms
clc

%A = csr2dense(val,col,row);
fprintf('Matrix Conversion\n');
tic
A = csr2csc(val,col,row);
toc
fprintf('Matrix Conversion Done\n\n');

M = speye(length(A))*(1/.6);

fprintf('Biconjugate gradients stabilized method\n');
tic
[x, bg_f, bg_rr, bg_itr] = bicgstab(A,d',1e-6,1000);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bg_rr,bg_itr);
fprintf(text);

fprintf('Biconjugate gradients stabilized method with PC\n');
tic
[x, bg_f, bg_rr, bg_itr] = bicgstab(A,d',1e-6,1000,M);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bg_rr,bg_itr);
fprintf(text);

fprintf('Generalized minimum residual method (with restarts)\n');
tic
[x, gm_f, gm_rr, gm_itr] = gmres(A,d',100,1e-6);
toc
output_text = 'Res: %1.3e\nItr: [%d %d]\n\n';
text = sprintf(output_text,gm_rr, gm_itr);
fprintf(text);

fprintf('Generalized minimum residual method (with restarts) with PC\n');
tic
[x, gm_f, gm_rr, gm_itr] = gmres(A,d',100,1e-6,100,M);
toc
output_text = 'Res: %1.3e\nItr: [%d %d]\n\n';
text = sprintf(output_text,gm_rr, gm_itr);
fprintf(text);

fprintf('Biconjugate gradients stabilized (1) method\n');
tic
[x, bgl_f, bgl_rr, bgl_itr] = bicgstabl(A,d',1e-6,1000);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bgl_rr,bgl_itr);
fprintf(text);

fprintf('Biconjugate gradients stabilized (1) method with PC\n');
tic
[x, bgl_f, bgl_rr, bgl_itr] = bicgstabl(A,d',1e-6,1000,M);
toc
output_text = 'Res: %1.3e\nItr: %d\n\n';
text = sprintf(output_text,bgl_rr,bgl_itr);
fprintf(text);

fprintf('MATLAB mldivide\n');
tic
x = d/A;
toc
fprintf('\n\n')
