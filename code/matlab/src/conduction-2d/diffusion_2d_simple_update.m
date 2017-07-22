% --- 2D Conduction ---
%
% Example 7.2 of Versteeg 2007

% Clear all previous data
close all
clear global
clc

% Define Initial and Boundary Conditions
qw = 500*1e3;   % Heat Flux, West kW/m^2
Tn = 100;   % Constant Temperature, North °C
k = 385;   % Thermal Conductivity W/m.K

alpha = 0.5;

itr = 0;

% Define Grid Size
X = 0.3;
Y = 0.4;
Z = 0.01;

m = 120;
n = 160;

dx = X/m;
dy = Y/n;

Aw = dy*Z;
Ae = dy*Z;
As = dx*Z;
An = dx*Z;

% Define Coefficients for Interior Nodes
Su_i = 0;
aw_i = k*Aw/dx;
ae_i = k*Ae/dx;
as_i = k*As/dy;
an_i = k*An/dy;
ap_i = aw_i+ae_i+as_i+an_i;

% Define Coefficients for WS Node
Su_ws = qw*Aw;
Sp_ws = 0;
aw_ws = 0;
ae_ws = k*Ae/dx;
as_ws = 0;
an_ws = k*An/dy;
ap_ws = aw_ws+ae_ws+as_ws+an_ws-Sp_ws;

% Define Coefficients for W Wall
Su_w = qw*Aw;
Sp_w = 0;
aw_w = 0;
ae_w = k*Ae/dx;
as_w = k*As/dy;
an_w = k*An/dy;
ap_w = aw_w+ae_w+as_w+an_w-Sp_w;

% Define Coefficients for WN Node
Su_wn = qw*Aw+(2*k*An/dy)*Tn;
Sp_wn = -2*k*An/dy;
aw_wn = 0;
ae_wn = k*Ae/dx;
as_wn = k*As/dy;
an_wn = 0;
ap_wn = aw_wn+ae_wn+as_wn+an_wn-Sp_wn;

% Define Coefficients for N Wall
Su_n = (2*k*An/dy)*Tn;
Sp_n = -2*k*An/dy;
aw_n = k*Aw/dx;
ae_n = k*Ae/dx;
as_n = k*As/dy;
an_n = 0;
ap_n = aw_n+ae_n+as_n+an_n-Sp_n;

% Define Coefficients for NE Node
Su_ne = (2*k*An/dy)*Tn;
Sp_ne = -2*k*An/dy;
aw_ne = k*Aw/dx;
ae_ne = 0;
as_ne = k*As/dy;
an_ne = 0;
ap_ne = aw_ne+ae_ne+as_ne+an_ne-Sp_ne;

% Define Coefficients for E Wall
Su_e = 0;
Sp_e = 0;
aw_e = k*Aw/dx;
ae_e = 0;
as_e = k*As/dy;
an_e = k*An/dy;
ap_e = aw_e+ae_e+as_e+an_e-Sp_e;

% Define Coefficients for ES Node
Su_es = 0;
Sp_es = 0;
aw_es = k*Aw/dx;
ae_es = 0;
as_es = 0;
an_es = k*An/dy;
ap_es = aw_es+ae_es+as_es+an_es-Sp_es;

% Define Coefficients for S Wall
Su_s = 0;
Sp_s = 0;
aw_s = k*Aw/dx;
ae_s = k*Ae/dx;
as_s = 0;
an_s = k*An/dy;
ap_s = aw_s+ae_s+as_s+an_s-Sp_s;

% Set up Temperature Matrix
T = zeros(n,m);
err_m = zeros(n,m)+1;
a = zeros(1,n);
b = zeros(1,n);
c = zeros(1,n);
d = zeros(1,n);

err_m(1,:) = 0;
err_m(n,:) = 0;
err_m(:,1) = 0;
err_m(:,m) = 0;

T_temp = T;
kappa=1;

% Start TDMA Solution Loops - North to South
err = 1;

tic 
while err > 1e-5
%   West i-Loop (west nodes)
    i = 1;    
    j = 1; 
    
    a(j) = 0;
    b(j) = ap_ws;
    c(j) = -an_ws;
    d(j) = ae_ws*T(j,i+1)+Su_ws;

    for j = 2:n-1
        a(j) = -as_w;
        b(j) = ap_w;
        c(j) = -an_w;
        d(j) = ae_w*T(j,i+1)+Su_w;
    end
    
    a(n) = -as_wn;
    b(n) = ap_wn;
    d(n) = ae_wn*T(n,i+1)+Su_wn;
    
    x = TDMAsolver(a,b,c,d);
    
    for k = 1:length(x)
        T(k,i) = x(k);
    end
    
%   Middle i-Loop (interior nodes)    
    for i = 2:m-1
            
        j = 1; 
    
        a(j) = 0;
        b(j) = ap_s;
        c(j) = -an_s;
        d(j) = aw_s*T(j,i-1)+ae_s*T(j,i+1)+Su_s;

        for j = 2:n-1
            a(j) = -as_i;
            b(j) = ap_i;
            c(j) = -an_i;
            d(j) = aw_i*T(j,i-1)+ae_i*T(j,i+1)+Su_i;
        end

        a(n) = -as_n;
        b(n) = ap_n;
        d(n) = aw_n*T(n,i-1)+ae_n*T(n,i+1)+Su_n;

        x = TDMAsolver(a,b,c,d);

        for k = 1:length(x)
            T(k,i) = x(k);
        end
    
    end
    
%   East i-Loop (east nodes)
    i = m;    
    j = 1; 
    
    a(j) = 0;
    b(j) = ap_es;
    c(j) = -an_es;
    d(j) = aw_es*T(j,i-1)+Su_es;

    for j = 2:n-1
        a(j) = -as_e;
        b(j) = ap_e;
        c(j) = -an_e;
        d(j) = aw_e*T(j,i-1)+Su_e;
    end
    
    a(n) = -as_ne;
    b(n) = ap_ne;
    d(n) = aw_ne*T(n,i-1)+Su_ne;
    
    x = TDMAsolver(a,b,c,d);
    
    for k = 1:length(x)
        T(k,i) = x(k);
    end
    
%   Convergence Check
    itr = itr+1;
    for i = 2:n-1
        for j = 2:m-1
            S_aT_nb = T(i,j-1)*aw_i+T(i,j+1)*aw_i+T(i-1,j)*as_i+T(i+1,j)*an_i;
            b = Su_i;
            err_m(i,j) = T(i,j)-T_temp(i,j)-alpha*((S_aT_nb+b)/ap_i-T_temp(i,j));
            tot_err(kappa) = err_m(i,j);
            kappa = kappa +1;
        end
    end
    err = max(max(abs(err_m)));
    T_temp = T;
end
toc

contourf(T,20)
colorbar
grid minor
colormap(jet)
title('2D Conduction 60x80 nodes')
i = 1:1:length(err);
