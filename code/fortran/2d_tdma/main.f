      program main

c     Written by: Matt Blomquist
c     Date: 2016-09-20

c     This program solves a 2D Conduction Problem using a line-by-line TDMA technique to 
c     solve the disretized equations.

      implicit none

c     Define Variables Used
      real :: X, Y, Z, m, n, dx, dy
      real :: Aw, Ae, As, An
      real :: err

c     Define Coefficients Used      
      real :: Su_i, aw_i, ae_i, as_i, an_i, ap_i
      real :: Su_ws, Sp_ws, aw_ws, ae_ws, as_ws, an_ws, ap_ws 
      real :: Su_w, Sp_w, aw_w, ae_w, as_w, an_w, ap_w 
      real :: Su_wn, Sp_wn, aw_wn, ae_wn, as_wn, an_wn, ap_wn 
      real :: Su_n, Sp_n, aw_n, ae_n, as_n, an_n, ap_n 
      real :: Su_ne, Sp_ne, aw_ne, ae_ne, as_ne, an_ne, ap_ne 
      real :: Su_e, Sp_e, aw_e, ae_e, as_e, an_e, ap_e
      real :: Su_es, Sp_es, aw_es, ae_es, as_es, an_es, ap_es
      real :: Su_s, Sp_s, aw_s, ae_s, as_s, an_s, ap_s

c     Define 
      real :: qw, Tn, k

c     Define Empty Arrays
      real, dimension(n,m) :: T, T_temp          
      real, dimension(n) :: a, b, c, d

c     Define Loop Constants
      integer :: i, j, itr

c     Define Initial and Boundary Conditions   
c     Heat Flux, West kW/m^2
      qw = 500000

c     Constant Temperature, North ∞C
      Tn = 100   

c     Thermal Conductivity W/m.K
      k = 385   

      itr = 0

c     Define Grid Size
      X = 0.3
      Y = 0.4
      Z = 0.01

      m = 30
      n = 40

      dx = X/m
      dy = Y/n

      Aw = dy*Z
      Ae = dy*Z
      As = dx*Z
      An = dx*Z

c     Define Coefficients for Interior Nodes
      Su_i = 0
      aw_i = k*Aw/dx
      ae_i = k*Ae/dx
      as_i = k*As/dy
      an_i = k*An/dy
      ap_i = aw_i+ae_i+as_i+an_i

c     Define Coefficients for WS Node
      Su_ws = qw*Aw
      Sp_ws = 0  
      aw_ws = 0   
      ae_ws = k*Ae/dx
      as_ws = 0   
      an_ws = k*An/dy
      ap_ws = aw_ws+ae_ws+as_ws+an_ws-Sp_ws

c     Define Coefficients for W Wall
      Su_w = qw*Aw
      Sp_w = 0
      aw_w = 0
      ae_w = k*Ae/dx
      as_w = k*As/dy
      an_w = k*An/dy
      ap_w = aw_w+ae_w+as_w+an_w-Sp_w

c     Define Coefficients for WN Node
      Su_wn = qw*Aw+(2*k*An/dy)*Tn
      Sp_wn = -2*k*An/dy
      aw_wn = 0   
      ae_wn = k*Ae/dx
      as_wn = k*As/dy
      an_wn = 0   
      ap_wn = aw_wn+ae_wn+as_wn+an_wn-Sp_wn

c     Define Coefficients for N Wall
      Su_n = (2*k*An/dy)*Tn
      Sp_n = -2*k*An/dy
      aw_n = k*Aw/dx
      ae_n = k*Ae/dx
      as_n = k*As/dy
      an_n = 0    
      ap_n = aw_n+ae_n+as_n+an_n-Sp_n

c     Define Coefficients for NE Node
      Su_ne = (2*k*An/dy)*Tn
      Sp_ne = -2*k*An/dy
      aw_ne = k*Aw/dx
      ae_ne = 0   
      as_ne = k*As/dy
      an_ne = 0   
      ap_ne = aw_ne+ae_ne+as_ne+an_ne-Sp_ne

c     Define Coefficients for E Wall
      Su_e = 0    
      Sp_e = 0    
      aw_e = k*Aw/dx
      ae_e = 0    
      as_e = k*As/dy
      an_e = k*An/dy
      ap_e = aw_e+ae_e+as_e+an_e-Sp_e

c     Define Coefficients for ES Node
      Su_es = 0   
      Sp_es = 0   
      aw_es = k*Aw/dx
      ae_es = 0   
      as_es = 0   
      an_es = k*An/dy
      ap_es = aw_es+ae_es+as_es+an_es-Sp_es

c     Define Coefficients for S Wall
      Su_s = 0    
      Sp_s = 0    
      aw_s = k*Aw/dx
      ae_s = k*Ae/dx
      as_s = 0    
      an_s = k*An/dy
      ap_s = aw_s+ae_s+as_s+an_s-Sp_s

c     Set up Temperature Matrix
      T = 0.0   
      a = 0.0      
      b = 0.0      
      c = 0.0      
      d = 0.0

      T_temp = T

c     Start TDMA Solution Loops - North to South
      err = 1     

c     tic 
      while err > .001
c           West i-Loop (west nodes)
            i = 1    
            j = 1 
    
            a(j) = 0
            b(j) = ap_ws
            c(j) = -an_ws
            d(j) = ae_ws*T(j,i+1)+Su_ws

            for j = 2:n-1
                  a(j) = -as_w
                  b(j) = ap_w
                  c(j) = -an_w
                  d(j) = ae_w*T(j,i+1)+Su_w
            end
    
            a(n) = -as_wn
            b(n) = ap_wn
            d(n) = ae_wn*T(n,i+1)+Su_wn
    
            x = TDMAsolver(a,b,c,d)
    
            for k = 1:length(x)
                  T(k,i) = x(k)    
            end
    
c           Middle i-Loop (interior nodes)    
            for i = 2:m-1
            
                  j = 1  
    
                  a(j) = 0    
                  b(j) = ap_s
                  c(j) = -an_s
                  d(j) = aw_s*T(j,i-1)+ae_s*T(j,i+1)+Su_s

                  for j = 2:n-1
                        a(j) = -as_i
                        b(j) = ap_i
                        c(j) = -an_i
                        d(j) = aw_i*T(j,i-1)+ae_i*T(j,i+1)+Su_i
                  end

                  a(n) = -as_n
                  b(n) = ap_n
                  d(n) = aw_n*T(n,i-1)+ae_n*T(n,i+1)+Su_n

                  x = TDMAsolver(a,b,c,d) 

                  for k = 1:length(x)
                        T(k,i) = x(k)     
                  end
    
            end
    
c           East i-Loop (east nodes)
            i = m    
            j = 1  
    
            a(j) = 0    
            b(j) = ap_es
            c(j) = -an_es
            d(j) = aw_es*T(j,i-1)+Su_es

            for j = 2:n-1
                  a(j) = -as_e
                  b(j) = ap_e
                  c(j) = -an_e
                  d(j) = aw_e*T(j,i-1)+Su_e
            end
    
            a(n) = -as_ne
            b(n) = ap_ne
            d(n) = aw_ne*T(n,i-1)+Su_ne
    
            x = TDMAsolver(a,b,c,d) 
    
            for k = 1:length(x)
                  T(k,i) = x(k)     
            end
    
            itr = itr+1
            err = mean(mean(abs(T-T_temp)))     
            T_temp = T
      end

c     toc

c     Write Data to a File

      end program main