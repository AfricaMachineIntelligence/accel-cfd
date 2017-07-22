module cfdSolve

implicit none

contains

! This subroutine calculates the number of non-zero values in the 
! matrix A.
subroutine calc_var_nums(t_vars,m,n,l)
      
implicit none
      
integer, intent(in) :: m,n,l
integer, intent(out) :: t_vars
      
integer :: i_nodes, l_nodes, s_nodes, c_nodes
integer :: i_vars, s_vars, l_vars, c_vars
      
! Determine Number of Nodes
i_nodes = (m-2)*(n-2)*(l-2)
l_nodes = 4*((m-2)+(n-2)+(l-2))
s_nodes = 2*((m-2)*(n-2)+(n-2)*(l-2)+(l-2)*(m-2))
c_nodes = 8

! Calculate Total A matrix Variables
i_vars = i_nodes*7
s_vars = s_nodes*6
l_vars = l_nodes*5
c_vars = c_nodes*4

t_vars = i_vars+l_vars+s_vars+c_vars
      
return t_vars

end subroutine calc_var_nums

! This subroutine calculates the index values for the CSR format.
subroutine calc_index_3d(indx,i,j,k,m,n,l)
      
integer, intent(in) :: i,j,k,m,n,l
integer, intent(out) :: indx
      
integer :: dm, dn, dl, 
integer :: b_clc, b_lsl, b_sis
integer :: plane1, plane2
      
! Create Abbriviation for dnterior lengths
dm = m-2
dn = n-2
dl = l-2

! Create Abbrivation for blocks
b_clc = 2*4+dm*5
b_lsl = 2*5+dm*6
b_sis = 2*6+dm*7
plane1 = 2*b_clc+dn*b_lsl
plane2 = 2*b_lsl+dn*b_sis

! Run Index Check and Calculate
if (k == 1) then
if (j == 1) then
if (i == 1) then
indx = 1  
else
indx = 5*(i-2)+4+2
end if
else if (j < n) then
if (i == 1) then
indx = 2+(j-2)*b_lsl+b_clc
else
indx = 6*(i-2)+8+(j-2)*b_lsl+b_clc
end if
else
if (i == 1) then
indx = 2+dn*b_lsl+b_clc
else
indx = 5*(i-2)+7+dn*b_lsl+b_clc
end if
end if
    
else if (k < l) then
if (j == 1) then
if (i == 1) then
indx = 2+(k-2)*plane2+plane1  
else
indx = 6*(i-2)+8+(k-2)*plane2+plane1
end if
else if (j < n) then
if (i == 1) then
indx = 3+(j-2)*b_sis+b_lsl+(k-2)*plane2+plane1
else
indx = 7*(i-2)+10+(j-2)*b_sis+b_lsl+(k-2)*plane2+plane1
end if
else
if (i == 1) then
indx = 3+dn*b_sis+b_lsl+(k-2)*plane2+plane1
else
indx = 6*(i-2)+9+dn*b_sis+b_lsl+(k-2)*plane2+plane1
end if
end if
    
else
if (j == 1) then
if (i == 1) then
indx = 2+dl*plane2+plane1  
else
indx = 5*(i-2)+7+dl*plane2+plane1
end if
else if (j < n) then
if (i == 1) then
indx = 3+(j-2)*b_lsl+b_clc+dl*plane2+plane1
else
indx = 6*(i-2)+9+(j-2)*b_lsl+b_clc+dl*plane2+plane1
end if
else
if (i == 1) then
indx = 3+dn*b_lsl+b_clc+dl*plane2+plane1
else
indx = 5*(i-2)+8+dn*b_lsl+b_clc+dl*plane2+plane1
end if
end if
end if
      
end subroutine calc_index_3d

! This subroutine calculates the sparse matrix from the coefficients
subroutine make_sparse_matrix(val,col,row,d,Ap,Ae,Aw,An,As,At,Ab,bp,m,n,l)
      
integer, intent(in) :: m, n, l
real, intent(in), dimension(m,n,l) :: Ap, Ae, Aw, An, As, At, Ab, bp
real, intent(out), dimension(m*n*l) :: row, d
real, intent(out), dimension() :: val, col
      

! Determine Number of Nodes
i_nodes = (m-2)*(n-2)*(l-2)
l_nodes = 4*((m-2)+(n-2)+(l-2))
s_nodes = 2*((m-2)*(n-2)+(n-2)*(l-2)+(l-2)*(m-2))
c_nodes = 8

! Calculate Total A matrix Variables
i_vars = i_nodes*7
s_vars = s_nodes*6
l_vars = l_nodes*5
c_vars = c_nodes*4

t_vars = i_vars+l_vars+s_vars+c_vars

val = zeros(1,t_vars)
row = zeros(1,t_vars)
col = zeros(1,t_vars)

d = zeros(1,m*n*l)

indx = zeros(m,n,l)

! Calculate indx Values
do i = 1:1:m
do j = 1:1:n
do k = 1:1:l
call calc_index_3d(index,i,j,k,m,n,l)
indx(i,j,k) = index                  
end do
end do
end do

! -------------------------------------------------------------------------
! Fill Interior Nodes
do k = 2:1:l-1
do j = 2:1:n-1
do i = 2:1:m-1
            
            
! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)
            
! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)
            
! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 
            
! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
            
! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)
            
! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)
            
! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
            
end do
end do
end do
! -------------------------------------------------------------------------
! Fill West Wall (i = 1)
do k = 2:1:l-1
do j = 2:1:n-1
        
i = 1
        
! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-1) = -As(i,j,k) 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
        
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
        
end do
end do
! -------------------------------------------------------------------------
! Fill East Wall (i = m)
do k = 2:1:l-1
do j = 2:1:n-1
            
i = m
        
! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+1) = -An(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
            
end do
end do
! -------------------------------------------------------------------------
! Fill South Wall (j = 1)
do k = 2:1:l-1
do i = 2:1:m-1
        
j = 1
        
! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)

end do
end do
! -------------------------------------------------------------------------
! Fill North Wall (j = n)
do k = 2:1:l-1
do i = 2:1:m-1
        
j = n

! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)

end do
end do
! -------------------------------------------------------------------------
! Fill Bottom Wall (k = 1)
do j = 2:1:n-1
do i = 2:1:m-1

k = 1
            
! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)

end do
end do
! -------------------------------------------------------------------------
! Fill Top Wall (k = l)
do j = 2:1:n-1
do i = 2:1:m-1

k = l

! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)

end do
end do
! -------------------------------------------------------------------------
! Fill in Corner (W,S,B)
i = 1
j = 1
k = 1

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Corner (E,S,B)
i = m
j = 1
k = 1

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+1) = -An(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Corner (W,N,B)
i = 1
j = n
k = 1

! South Node
val(indx(i,j,k)-1) = -As(i,j,k) 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Corner (E,N,B)
i = m
j = n
k = 1

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+1) = -At(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Corner (W,S,T)
i = 1
j = 1
k = l

! Bottom Node
val(indx(i,j,k)-1) = -Ab(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-2)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Corner (E,S,T)
i = m
j = 1
k = l

! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+1) = -An(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Corner (W,N,T)
i = 1
j = n
k = l

! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-1) = -As(i,j,k) 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Corner (E,N,T)
i = m
j = n
k = l

! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
! -------------------------------------------------------------------------
! Fill in Line (W2E - S,B)
do i = 2:1:m-1
j = 1
k = 1

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
end do
! -------------------------------------------------------------------------
! Fill in Line (W2E - N,B)
do i = 2:1:m-1
j = n
k = 1

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)
    
! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
end do
! -------------------------------------------------------------------------
! Fill in Line (W2E - S,T)
do i = 2:1:m-1
j = 1
k = l

! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
end do
! -------------------------------------------------------------------------
! Fill in Line (W2E - N,T)
do i = 2:1:m-1
j = n
k = l
    
! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
end do
! -------------------------------------------------------------------------
! Fill in Line (S2N - W,B)
do j = 2:1:n-1
i = 1
k = 1

! South Node
val(indx(i,j,k)-1) = -As(i,j,k) 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------
! Fill in Line (S2N - E,B)
do j = 2:1:n-1
i = m
k = 1

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+1) = -An(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------
! Fill in Line (S2N - W,T)
do j = 2:1:n-1
i = 1
k = l

! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-1) = -As(i,j,k) 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------
! Fill in Line (S2N - E,T)
do j = 2:1:n-1
i = m
k = l
    
! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+1) = -An(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------
! Fill in Line (B2T - W,S)
do k = 2:1:l-1
i = 1
j = 1

! Bottom Node
val(indx(i,j,k)-1) = -Ab(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-2)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+2) = -An(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+3) = -At(i,j,k)
row(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+3) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------
! Fill in Line (B2T - E,S)
do k = 2:1:l-1
i = m
j = 1
    
! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! North Node
val(indx(i,j,k)+1) = -An(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------
! Fill in Line (B2T - W,N)
do k = 2:1:l-1
i = 1
j = n
    
! Bottom Node
val(indx(i,j,k)-2) = -Ab(i,j,k)
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-1) = -As(i,j,k) 
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-1) = i+m*(j-2)+m*n*(k-1)

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! East Node
val(indx(i,j,k)+1) = -Ae(i,j,k) 
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = (i+1)+m*(j-1)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+2) = -At(i,j,k)
row(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+2) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------
! Fill in Line (B2T - E,N)
do k = 2:1:l-1
i = m
j = n
    
! Bottom Node
val(indx(i,j,k)-3) = -Ab(i,j,k)
row(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-3) = i+m*(j-1)+m*n*(k-2)

! South Node
val(indx(i,j,k)-2) = -As(i,j,k) 
row(indx(i,j,k)-2) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)-2) = i+m*(j-2)+m*n*(k-1)

! West Node
val(indx(i,j,k)-1) = -Aw(i,j,k)
row(indx(i,j,k)-1) = i+m*(j-1)+m*n*(k-1) 
col(indx(i,j,k)-1) = (i-1)+m*(j-1)+m*n*(k-1) 

! Center Node (p)
val(indx(i,j,k)) = Ap(i,j,k)
row(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)) = i+m*(j-1)+m*n*(k-1)

! Top Node
val(indx(i,j,k)+1) = -At(i,j,k)
row(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k-1)
col(indx(i,j,k)+1) = i+m*(j-1)+m*n*(k)
            
! Solution Value
d(i+m*(j-1)+m*n*(k-1)) = bp(i,j,k)
    
end do
! -------------------------------------------------------------------------

end subroutine make_sparse_matrix

! This subroutine solves the sparse linear system using the PARDISO routine.
subroutine solvePardiso()

end subroutine solvePardiso

subroutine solveBicg()

end subroutine solveBicg

subroutine solveBicgstab()

end subroutine solveBicgstab

subroutine solveGmres()

end subroutine solveGmres

end module cfdSolve
