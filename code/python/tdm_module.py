import numpy as np

def TDMAsolver(a, b, c, d):
    # a, b, and c are the main diagonals of a tri-diagonal matrix with a representing the lower
    # diagonal and c representing the upper diagonal. d is the right vector.
    # this function returns x, which is the solution to Ax = d
    n = len(b)
    x = np.zeros(n)

    c[0] = c[0] / b[0]
    d[0] = d[0] / b[0] 

    for i in range(1, n-1, 1):
        temp = b[i] - a[i] * c[i-1]
        c[i] = c[i] / temp
        d[i] = (d[i] - a[i] * d[i-1]) / temp

    d[n-1] = (d[n-1] - a[n-1] * d[n-2]) / (b[n-1] - a[n-1] * c[n-2])

    x[n-1] = d[n-1]

    for k in range(n-2, -1, -1):
        x[k] = d[k] - c[k] * x[k+1]
    #print(x)
    return x 

def cond_2d(A,T,err,r,n,m):
    # This function performs a TDMA for a 2D grid.
    
    e = np.zeros((n, m))+1.
    e[:,0] = 0.
    e[:,m-1] = 0.
    e[0,:] = 0.
    e[n-1,:] = 0.

    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    d = np.zeros(n)
    x = np.zeros(n)
    
    T_temp = T

    err_max = 1.
    while (err_max > err):

    
        # West i-Loop (west nodes)
        i = 0;    
        j = 0; 
    
        a[j] = 0
        b[j] = A[6,8]
        c[j] = -A[5,8]
        d[j] = A[3,8]*T[j,i+1]+A[0,8]

        for j in range(1, n-1, 1):
            a[j] = -A[4,1]
            b[j] = A[6,1]
            c[j] = -A[5,1]
            d[j] = A[3,1]*T[j,i+1]+A[0,1]
        
    
        a[n-1] = -A[4,2]
        b[n-1] = A[6,2]
        d[n-1] = A[3,2]*T[n-1,i]+A[0,2]
    
        x = TDMAsolver(a,b,c,d)
    
        for k in range(0, len(x)):
            T[k,i] = x[k]
        
    
        # Middle i-Loop (interior nodes)    
        for i in range(1, m-2):
            
            j = 0
    
            a[j] = 0
            b[j] = A[6,7]
            c[j] = -A[5,7]
            d[j] = A[2,7]*T[j,i-1]+A[0,7]

            for j in range(1, n-1, 1):
                a[j] = -A[4,0]
                b[j] = A[6,0]
                c[j] = A[5,0]
                d[j] = A[2,0]*T[j,i-1]+A[3,0]*T[j,i+1]+A[0,0]
            
            a[n-1] = -A[4,3]
            b[n-1] = A[6,3]
            d[n-1] = A[2,3]*T[n-1,i-1]+A[3,3]*T[n-1,i+1]+A[0,3]

            x = TDMAsolver(a,b,c,d)
            for k in range(0, len(x)):
                T[k,i] = x[k]
            
    
        # East i-Loop (east nodes)
        i = m-1    
        j = 0 
    
        a[j] = 0
        b[j] = A[6,6]
        c[j] = -A[5,6]
        d[j] = A[2,6]*T[j,i-1]+A[0,6]

        for j in range(1, n-1, 1):
            a[j] = -A[4,5]
            b[j] = A[6,5]
            c[j] = -A[5,5]
            d[j] = A[2,5]*T[j,i-1]+A[0,5]
            
        a[n-1] = -A[4,4]
        b[n-1] = A[0,4]
        d[n-1] = A[2,4]*T[n-1,i-1]+A[0,4]
    
        x = TDMAsolver(a,b,c,d)
    
        for k in range(0, len(x)):
            T[k,i] = x[k]
        
   
        # Convergence Check
        for i in range(1, n-1):
            for j in range(1, m-1):
                S_aT_nb = T[i,j-1]*A[2,0]+T[i,j+1]*A[3,0]+T[i-1,j]*A[4,0]+T[i+1,j]*A[5,0]
                e[i,j] = T[i,j] - T_temp[i,j] - r * ((S_aT_nb+A[0,0]) / A[6,0] - T_temp[i,j])
                
        err_max = np.amax(np.absolute(e))
        #print(err_max)
        #print(np.argmax(np.absolute(e)))
        T_temp = T
        #print(e)

    return T