close all
clear all
clc

% Set the size of the geometric grid
N = 3;

% Determine Number of Variables
i_nodes = (N-2)^2;
s_nodes = 4*(N-2);
c_nodes = 4;

% Calculate Total A matrix Variables
i_vars = i_nodes*5;
s_vars = s_nodes*4;
c_vars = c_nodes*3;

t_vars = i_vars+s_vars+c_vars;

n_perc = t_vars/(N^2)^2;

line1_text = 'Total number of A matrix (2D) vars for %3.0f is : %3.0f.\n';
line2_text = 'Sparsity is: %0.5f\n';

fprintf(line1_text,N,t_vars)
fprintf(line2_text,n_perc)
