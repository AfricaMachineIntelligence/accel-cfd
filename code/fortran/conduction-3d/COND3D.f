     program COND3D
! This program solves a 3D conduction problem base on input.
! ---------------------------------------------------------------
! This program solves a 3D conduction problem with BC as follows:
! Wast Wall = 
! East Wall = 
! North Wall = 
! South Wall = 
! Top Wall = 
! Bottom Wall = 
! ----------------------------------------------------------------
! The following constants are also implemented
! Nodes in m, n, and l are 100, 50, 25
      
      implicit none

! Initialize Parameters 
	  
      real(8) :: xdim, ydim, zdim
      real(8) :: hcon, e, r
      real(8) :: t1,t2
	  
      integer :: m, n, l
      integer :: i, j, k

      character(len=8) :: west_type, north_type, east_type
      character(len=8) :: south_type, top_type, bottom_type

      real(8) :: west_val, north_val, east_val, south_val
      real(8) :: top_val, bottom_val

      real(8) :: dx, dy, dz
      real(8) :: Aw, An, Ae, As, At, Ab
      real(8), dimension(27,9) :: CMAT
      real(8), dimension(6) :: BC

      character(len=4), dimension(6) :: BC_TYPE

      real(8), dimension(300,200,100) :: T

!	-----------------------------------------------------------------
!	DEFINE PARAMETERS

!	Grid Dimensions (x,y,z) // (m)
      xdim = 0.0100
      ydim = 0.0050
      zdim = 0.0025

!	Heat Transfer Coefficient // (W/m2/K)
      hcon = 385 

!	Solution Parameteres
!	Number of Nodes in Gris (m > x, n > y, l > z) // (no dim)
      m = 300
      n = 200
      l = 100

!     Error Threshold for Convergence // (no dim)
      e = 0.005

!	Relaxation Factor for Convergence // (no dim)
      r = 0.0

!	Bondary Conditions // Type = 'COND' or 'INSU' or 'TEMP'
      west_type = 'COND'
      west_val = 50000

      north_type = 'TEMP'
      north_val = 100

      east_type = 'INSU'
      east_val = 0
	  
      south_type = 'INSU'
      south_val = 0

      top_type = 'INSU'
      top_val = 0

      bottom_type = 'INSU'
      bottom_val = 0
			 
!	END DEFINE PARAMETERS
!	-----------------------------------------------------------------

!	-----------------------------------------------------------------
!	DEFINE NON-DIM VALUES
      dx = xdim/m
      dy = ydim/n
      dz = zdim/l

!	Create Boundary Condition Array for Type
      BC_TYPE(1) = west_type
      BC_TYPE(2) = north_type
      BC_TYPE(3) = east_type
      BC_TYPE(4) = south_type
      BC_TYPE(5) = top_type
      BC_TYPE(6) = bottom_type
	   
!	Create Boundary Condition Array for Subroutine Feedthrough
      BC(1) = west_val
      BC(2) = north_val
      BC(3) = east_val
      BC(4) = south_val
      BC(5) = top_val
      BC(6) = bottom_val

!	Initialize Coefficients Matrix
      CMAT = 0
      T = 0

!	Create Coefficients Matrix
      call create_q(CMAT, dx, dy, dz, hcon, BC_TYPE, BC)

!     Start CPU Time Count
      call CPU_TIME(t1)
      print *, 'Start TDMA:',t1

!	Call TDMA Loop to update T
      call TDMA3D(T, CMAT, m, n, l, e, r)

!	End CPU Time Count
      call CPU_TIME(t2)
      print *, 'End TDMA:',t2
      print *, 'Total Duration of TDMA Loop:',t2-t1

!     Write Results of T to file
      open(unit=7, file='temp_large.txt')
      
      do k = 1,l
            do j = 1,n
                  do i = 1,m
                        write(7,*), T(i,j,k)
                  end do
            end do
      end do
      
      close(7)
      
      end program COND3D  

!	-----------------------------------------------------------------
!	-----------------------------------------------------------------
!	-----------------------------------------------------------------
!	-----------------------------------------------------------------
!	-----------------------------------------------------------------

      subroutine TDMA3D(T,Q,m,n,l,e,r)
      
      implicit none
      
! Define Input Variables that we need
      integer, intent(in) :: m, n, l
      real(8), dimension(m,n,l), intent(out) :: T
      real(8), dimension(27,9), intent(in) :: Q
      real(8), intent(in) :: e, r
      
! Define Variables that are used just in this subroutine
      real(8), dimension(m) :: a, b, c, d, x
      real(8), dimension(m,n,l) :: E_MAT, T_TEMP
      real(8) :: E_MAX
      integer :: i, j, k, itr
            
! Initialize New Variables
      a = 0
      b = 0
      c = 0
      d = 0
      x = 0
      E_MAT = 1
      E_MAX = 1
      
! Set T_TEMP to T
      T_TEMP = T
      
! Create an Iteration Count Variable 
      itr = 1
      
! This is the main T matrix update and convergence loop
      do while (E_MAX > e)

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! Start >> Perform 2D T update with Z = 1 (k=1)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! Start >> Bottom Plane Sweep (Z-Min)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

! Initialize i,j,k
      i = 1
      j = 1
      k = 1

! ---------------------------------------------------------------------
! Start >> Bottom-South Edge Solution
! ---------------------------------------------------------------------

! South-West-Bottom Corner
      a(i) =  0
      b(i) = Q(27,9)
      c(i) = -Q(27,5)
      d(i) = Q(27,4)*T(i,j+1,k)+Q(27,7)*T(i,j,k+1)+Q(27,1)

      !	Bottom-South Edge 
      do i = 2, m-1      
            a(i) = -Q(19,3)
            b(i) = Q(19,9)
            c(i) = -Q(19,5)
            d(i) = Q(19,4)*T(i,j+1,k)+Q(19,7)*T(i,j,k+1)+Q(19,1)       
      end do

! Set i to m // Last node in line.
      i = m

! East-South-Bottom Corner
      a(i) = -Q(25,3)
      b(i) = Q(25,9)
      c(i) = 0
      d(i) = Q(25,4)*T(i,j+1,k)+Q(25,7)*T(i,j,k+1)+Q(25,1)

! ---------------------------------------------------------------------
! End >> Bottom-South Edge Solution
! ---------------------------------------------------------------------

! Solve 1D Line TDMA 
      call tdma_solver(c, b, a, d, m, x)
          
! Update T with 1D Line Solution
      T(1:m,j,k) = x

! ---------------------------------------------------------------------
! Start >> Bottom Plane Solution
! ---------------------------------------------------------------------

!	Start South to North Sweep 
      do j = 2,n-1
      
! Initialize i         
            i = 1

!	Bottom-West Edge
            a(i) = 0
            b(1) = Q(16,9)
            c(i) = -Q(16,5)
            d(i) = Q(16,4)*T(i,j+1,k)+Q(16,7)*T(i,j,k+1)+Q(16,6) &
                  *T(i,j-1,k)+Q(16,1)

! Bottom-Plane
            do i = 2, m-1
                  a(i) = -Q(7,3)
                  b(i) = Q(7,9)
                  c(i) = -Q(7,5)
                  d(i) = Q(7,4)*T(i,j+1,k)+Q(7,6)*T(i,j-1,k) &
                        +Q(7,7)*T(i,j,k+1)+Q(7,1)
            end do

! Set i to m // Last node in line.
            i = m

! Bottom-East Edge
            a(m) = -Q(18,3)
            b(m) = Q(18,9)
            c(m) = 0
            d(m) = Q(18,4)*T(i,j+1,k)+Q(18,6)*T(i,j-1,k)+Q(18,7) &
                  *T(i,j,k+1)+Q(18,1)

! ---------------------------------------------------------------------
! End >> Bottom Plane Solution
! ---------------------------------------------------------------------

! Solve 1D Line TDMA 
            call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
            T(1:m,j,k) = x

! Stop South to North Sweep
            end do

! ---------------------------------------------------------------------
! Start >> Bottom-North Edge Solution
! ---------------------------------------------------------------------

! Initialize i and set j to n
            i = 1
            j = n

! Bottom-North-West Corner
            a(i) = 0
            b(i) = Q(21,9)
            c(i) = -Q(21,5)
            d(i) = Q(21,6)*T(i,j-1,k)+Q(21,7)*T(i,j,k+1)+Q(21,1)

! Bottom-North Edge
            do i = 2, m-1
                  a(i) = -Q(17,3)
                  b(i) = Q(17,9)
                  c(i) = -Q(17,5)
                  d(i) = Q(17,6)*T(i,j-1,k)+Q(17,7)*T(i,j,k+1)+Q(17,1)
            end do

! Set i to m // Last node in line.
            i = m

! Bottom-North-East Corner
            a(i) = -Q(23,3)
            b(i) = Q(23,9)
            c(i) = 0
            d(i) = Q(23,6)*T(i,j-1,k)+Q(23,7)*T(i,j,k+1)+Q(23,1)

! ---------------------------------------------------------------------
! End >> Bottom-North Edge Solution
! ---------------------------------------------------------------------

! Solve 1D Line TDMA 
            call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
            T(1:m,j,k) = x

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! End >> Bottom Plane Sweep (Z-Min)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! *********************************************************************
! *********************************************************************
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! Start >> Bottom to Top Sweep (Z-Direction / Interior Nodes)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

! Start Bottom to Top Loop
            do k = 2,l-1

! ---------------------------------------------------------------------
! Start >> South Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

!	Initialize i and j
                  i = 1
                  j = 1

!	West-South Edge
                  a(i) = 0
                  b(i) = Q(11,9)
                  c(i) = -Q(11,5)
                  d(i) = Q(11,4)*T(i,j+1,k)+Q(11,7)*T(i,j,k+1)+Q(11,8) &
                  *T(i,j,k-1)+Q(11,1)

!	South Wall
                  do i = 2, m-1
                        a(i) = -Q(5,3)
                        b(i) = Q(5,9)
                        c(i) = -Q(5,5)
                        d(i) = Q(5,4)*T(i,j+1,k)+Q(5,7)*T(i,j,k+1) &
                              +Q(5,8)*T(i,j,k-1)+Q(5,1)
                  end do

! Set i to m // Last node in line.
                  i = m

! South-East Corner
                  a(m) = -Q(10,3)
                  b(m) = Q(10,9)
                  c(m) = 0
                  d(m) = Q(10,4)*T(i,j+1,k)+Q(10,7)*T(i,j,k+1)+Q(10,8) &
                        *T(i,j,k-1)+Q(10,1)

! Solve 1D Line TDMA 
                  call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
                  T(1:m,j,k) = x

! ---------------------------------------------------------------------
! End >> South Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> South to North Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! Start South to North Loop
                  do j = 2, n-1

!	Initialize i
                        i = 1

!	West Wall
                        a(i) = 0
                        b(i) = Q(2,9)
                        c(i) = -Q(2,5)
                        d(i) = Q(2,4)*T(i,j+1,k)+Q(2,6)*T(i,j-1,k) &
                              +Q(2,7)*T(i,j,k+1)+Q(2,8)*T(i,j,k-1) &
                              +Q(2,1)

! Interior Nodes 
                        do i = 2, m-1 
                              a(i) = -Q(1,3)
                              b(i) = Q(1,9)
                              c(i) = -Q(1,5)
                              d(i) = Q(1,4)*T(i,j+1,k)+Q(1,6) &
                                    *T(i,j-1,k)+Q(1,7)*T(i,j,k+1) &
                                    +Q(1,8)*T(i,j,k-1)+Q(1,1)
                        end do

! Set i to m // Last node in line.
                        i = m

! East Wall
                        a(m) = -Q(4,3)
                        b(m) = Q(4,9)
                        c(m) = 0
                        d(m) = Q(3,4)*T(i,j+1,k)+Q(3,6)*T(i,j-1,k) &
                              +Q(3,7)*T(i,j,k+1)+Q(3,8)*T(i,j,k-1) &
                              +Q(3,1)

! Solve 1D Line TDMA 
                        call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
                        T(1:m,j,k) = x

! End South to North Loop
                  end do

! ---------------------------------------------------------------------
! End >> South to North Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> North Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! Initialize i and j
                  i = 1
                  j = n

! West-North Edge
                  a(i) = 0
                  b(i) = Q(8,9)
                  c(i) = -Q(8,5)
                  d(i) = Q(8,6)*T(i,j-1,k)+Q(8,7)*T(i,j,k+1)+Q(8,8) &
                        *T(i,j,k-1)+Q(8,1)

! North Wall
                  do i = 2, m-1
                        a(i) = -Q(3,3)
                        b(i) = Q(3,9)
                        c(i) = -Q(3,5)
                        d(i) = Q(5,6)*T(i,j-1,k)+Q(5,7)*T(i,j,k+1) &
                              +Q(5,8)*T(i,j,k-1)+Q(5,1)
                  end do

! Set i to m // Last node in line.
                  i = m

!	North-East Edge
                  a(m) = -Q(9,3)
                  b(m) = Q(9,9)
                  c(m) = 0
                  d(m) = Q(9,6)*T(i,j-1,k)+Q(9,7)*T(i,j,k+1)+Q(9,8) &
                        *T(i,j,k-1)+Q(9,1)

! Solve 1D Line TDMA 
                  call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
                  T(1:m,j,k) = x

! ---------------------------------------------------------------------
! End >> North Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! End k loop
            end do

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! End >> Bottom to Top Sweep (Z-Direction / Interior Nodes)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! *********************************************************************
! *********************************************************************
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! Start >> Top Plane Sweep (Z-Max)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> Top-South Edge (Z-Max)
! ---------------------------------------------------------------------

! Initialize i, j, and k
            i = 1
            j = 1
            k = l

! Top-South-West Corner
            a(i) = 0
            b(i) = Q(26,9)
            c(i) = -Q(26,5)
            d(i) = Q(26,4)*T(i,j+1,k)+Q(26,8)*T(i,j,k-1)+Q(26,1)

! Top-South Edge 
            do i = 2, m-1
                  a(i) = -Q(15,3)
                  b(i) = Q(15,9) 
                  c(i) = -Q(15,5)
                  d(i) = Q(15,4)*T(i,j+1,k)+Q(15,8)*T(i,j,k-1)+Q(15,1)
            end do
! Set i to m // Last node in line.
            i = m

!	Top-East-South Corner
            a(m) = -Q(24,3)
            b(m) = Q(24,9)
            c(m) = 0
            d(m) = Q(24,4)*T(i,j+1,k)+Q(24,8)*T(i,j,k-1)+Q(24,1)

! Solve 1D Line TDMA 
            call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
            T(1:m,j,k) = x

! ---------------------------------------------------------------------
! End >> Top-South Edge (Z-Max)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> Top Plane (Z-Max)
! ---------------------------------------------------------------------

! Start South to North Loop
            do j = 2, n-1

! Initialize i
                  i = 1
              
! Top-West Wall
                  a(1) =  0
                  b(1) = Q(12,9)
                  c(1) = -Q(12,5)
                  d(1) = Q(12,4)*T(i,j+1,k)+Q(12,6)*T(i,j-1,k) &
                        +Q(12,8)*T(i,j,k-1)+Q(12,1)

! Top Wall
                  do i = 2, m-1
                        a(i) = -Q(6,3)
                        b(i) = Q(6,9)
                        c(i) = -Q(6,5)
                        d(i) = Q(6,4)*T(i,j+1,k)+Q(6,6)*T(i,j-1,k) &
                              +Q(6,8)*T(i,j,k-1)+Q(6,1)
                  end do

! Set i to m // Last node in line.
                  i = m

! Top-East Edge
                  a(m) = -Q(14,3)
                  b(m) = Q(14,9)
                  c(m) = 0
                  d(m) = Q(14,4)*T(i,j+1,k)+Q(14,6)*T(i,j-1,k) &
                        +Q(14,8)*T(i,j,k-1)+Q(14,1)

! Solve 1D Line TDMA 
                  call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
                  T(1:m,j,k) = x

! End j loop (North to South)
            end do

! ---------------------------------------------------------------------
! End >> Top Plane (Z-Max)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> North-Top Edge (Z-Max)
! ---------------------------------------------------------------------


! Initialize i and j
            i = 1
            j = n

! West-North-Top Corner
            a(i) =  0
            b(i) = Q(20,9)
            c(i) = -Q(20,5)
            d(i) = Q(20,6)*T(i,j-1,k)+Q(20,8)*T(i,j,k-1)+Q(20,1)

! North Top Edge
            do i = 2, m-1
                  a(i) = -Q(13,3)
                  b(i) = Q(13,9) 
                  c(i) = -Q(13,5)
                  d(i) = Q(13,6)*T(i,j-1,k)+Q(13,8)*T(i,j,k-1)+Q(13,1)
            end do

! Set i to m // Last node in line.
            i = m

! East-North-Top Corner
            a(m) = -Q(22,3)
            b(m) = Q(22,9)
            c(m) = 0
            d(m) = Q(22,6)*T(i,j-1,k)+Q(22,8)*T(i,j,k-1)+Q(22,1)

! Solve 1D Line TDMA 
            call tdma_solver(a, b, c, d, m, x)

! Update T with 1D Line Solution
            T(1:m,j,k) = x

! ---------------------------------------------------------------------
! End >> North-Top Edge (Z-Max)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! End >> Perform 2D T update with Z = 1 (k=1)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

! /////////////////////////////////////////////////////////////////////
! /////////////////////////////////////////////////////////////////////
! /////////////////////////////////////////////////////////////////////
! /////////////////////////////////////////////////////////////////////
! /////////////////////////////////////////////////////////////////////

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! Start >> Perform Convergence Check 
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

! *********************************************************************
! Start >> Bottom Plane Sweep
! *********************************************************************

! ---------------------------------------------------------------------
! Start >> South-Bottom Sweep (Z-Min)
! ---------------------------------------------------------------------

! Initialize i and j
            i = 1
            j = 1
            k = 1

! West-South-Bottom Corner
            E_MAT(1,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(27,5) &
                        *T(i+1,j,k)+Q(27,4)*T(i,j+1,k)+Q(27,7) &
                        *T(i,j,k+1))/Q(27,9)-T_TEMP(i,j,k))

! South Bottom Edge
            do i = 2,m-1
                  E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(19,5) &
                              *T(i+1,j,k)+Q(19,4)*T(i,j+1,k)+Q(19,7) &
                              *T(i,j,k+1)+Q(19,3)*T(i-1,j,k))/Q(19,9) &
                              -T_TEMP(i,j,k))
            end do

! Set i to m // Last node in line.
            i = m

! South-East-Bottom 
            E_MAT(m,j,k) =T(i,j,k)-T_TEMP(i,j,k)-r*((Q(25,3) &
                        *T(i-1,j,k)+Q(25,4)*T(i,j+1,k)+Q(25,7) &
                        *T(i,j,k+1))/Q(25,9)-T_TEMP(i,j,k))

! ---------------------------------------------------------------------
! End >> South-Bottom Plane Sweep (Z-Min)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> Bottom Plane Sweep (Z-Min)
! ---------------------------------------------------------------------

! Start South to North Loop
            do j = 2,n-1

! Initialize i
                  i = 1

! West-Bottom Edge
                  E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(16,5) &
                              *T(i+1,j,k)+Q(16,4)*T(i,j+1,k)+Q(16,7) &
                              *T(i,j,k+1)+Q(16,6)*T(i,j-1,k))/Q(16,9) & 
                              -T_TEMP(i,j,k))

! Bottom Plane
                  do i = 2,m-1
                        E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r &
                                    *((Q(7,3)*T(i-1,j,k)+Q(7,4) &
                                    *T(i,j+1,k)+Q(7,7)*T(i,j,k+1) &
                                    +Q(7,5)*T(i+1,j,k))/Q(7,9) &
                                    -T_TEMP(i,j,k))
                  end do

! Set i to m // Last node in line.
                  i = m

! East-bottom Edge
                  E_MAT(m,j,1) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(18,3) &
                              *T(i-1,j,k)+Q(18,4)*T(i,j+1,k)+Q(18,7) &
                              *T(i,j,k+1)+Q(18,6)*T(i,j-1,k))/Q(18,9) &
                              -T_TEMP(i,j,k))

! Stop South to North Loop
            end do

! ---------------------------------------------------------------------
! End >> Bottom Plane Sweep (Z-Min)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> Bottom-North Edge Sweep (Z-Min)
! ---------------------------------------------------------------------

! Initialize i and j
            i = 1
            j = n

! West-North-Bottom Corner 
            E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(21,5) &
                        *T(i+1,j,k)+Q(21,6)*T(i,j-1,k)+Q(21,7) &
                        *T(i,j,k+1))/Q(21,9)-T_TEMP(i,j,k))

! North Bottom Edge
            do i = 2,m-1
                  E_MAT(i,j,1) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(17,5) &
                              *T(i+1,j,k)+Q(17,6)*T(i,j-1,k)+Q(17,7) &
                              *T(i,j,k+1)+Q(19,3)*T(i-1,j,k))/Q(17,9) &
                              -T_TEMP(i,j,k))
            end do


! Set i to m // Last node in line.
            i = m

! North-East-Bottom Corner
            E_MAT(m,n,1) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(23,3) &
                        *T(i-1,j,k)+Q(23,6)*T(i,j-1,k)+Q(23,7) &
                        *T(i,j,k+1))/Q(23,9)-T_TEMP(i,j,k))

! ---------------------------------------------------------------------
! End >> Bottom-North Edge Sweep (Z-Min)
! ---------------------------------------------------------------------

! *********************************************************************
! End >> Bottom Plane Sweep ||| Start >> Interior Z-Plane Sweep 
! *********************************************************************

! Start k loop (Bottom to Top)
            do k = 2,l-1
            
! ---------------------------------------------------------------------
! Start >> South Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! Initialize i and j  
                  i = 1
                  j = 1
          
! West-South Edge
                  E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(11,5) &
                              *T(i+1,j,k)+Q(11,4)*T(i,j+1,k)+Q(11,7) &
                              *T(i,j,k+1)+Q(11,8)*T(i,j,k-1))/Q(11,9) &
                              -T_TEMP(i,j,k))

! South Plane
                  do i = 2,m-1
                        E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r &
                                    *((Q(5,5)*T(i+1,j,k)+Q(5,4) &
                                    *T(i,j+1,k)+Q(5,7)*T(i,j,k+1) &
                                    +Q(5,8)*T(i,j,k-1)+Q(5,3) &
                                    *T(i-1,j,k))/Q(5,9)-T_TEMP(i,j,k))
                  end do
                  
! Set i to m // Last node in line.
                  i = m
                  
! East-South Edge
                  E_MAT(m,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(10,5) &
                              *T(i,j+1,k)+Q(10,7)*T(i-1,j,k)+Q(10,3) &
			                *T(i,j,k+1)+Q(10,8)*T(i,j,k-1))/Q(10,9) &
                              -T_TEMP(i,j,k))

! ---------------------------------------------------------------------
! End >> South Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> Interior Node Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! Start j loop (South to North Sweep)
                  do j = 2,n-1
                        
! Initialize i
                        i = 1 
                  
! West Plane
                        E_MAT(1,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r &
                                    *((Q(2,5)*T(i+1,j,k)+Q(2,4) &
                                    *T(i,j+1,k)+Q(2,7)*T(i,j,k+1) &
                                    +Q(2,8)*T(i,j,k-1)+Q(2,6) &
                                    *T(i,j-1,k))/Q(2,9)-T_TEMP(i,j,k))

! Interior Nodes
                        do i = 2,m-1
                              E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)- &
                                          r*((Q(1,5)*T(i+1,j,k)+Q(1,4) &
                                          *T(i,j+1,k)+Q(1,7)*T(i,j,k+1) &
                                          + Q(1,8)*T(i,j,k-1)+Q(1,6) &
                                          *T(i,j-1,k)+Q(1,3) &
                                          *T(i-1,j,k))/Q(1,9) &
                                          -T_TEMP(i,j,k))
                        end do

! Set i to m // Last node in line.
                        i = m
            
! East Plane 
                        E_MAT(m,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r &
                                    *((Q(4,4)*T(i,j+1,k)+ Q(4,7) &
                                    * T(i,j,k+1)+Q(4,8)*T(i,j,k-1) &
                                    +Q(4,6)*T(i,j-1,k)+Q(4,3) &
                                    *T(i-1,j,k))/Q(4,9)-T_TEMP(i,j,k))

                  end do

! ---------------------------------------------------------------------
! End >> Interior Node Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------
                  
! ---------------------------------------------------------------------
! Start >> North Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! Initialize i and j
                  i = 1
                  j = n
              
! West-North Edge
                  E_MAT(i,j,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(12,4) &
                              *T(i+1,j,k)+Q(12,7)*T(i,j,k+1)+Q(12,8) &
                              *T(i,j,k-1)+Q(12,6)*T(i,j-1,k))/Q(12,9) &
                              - T_TEMP(i,j,k))

! North Plane
                  do i = 2,m-1
                        E_MAT(i,n,k) = T(i,j,k)-T_TEMP(i,j,k)-r &
                              *((Q(3,4)*T(i+1,j,k)+Q(3,3)*T(i,j,k+1) &
                              +Q(3,8)*T(i,j,k-1)+Q(3,6)*T(i,j-1,k) &
                              +Q(3,3)*T(i-1,j,k))/Q(3,9)-T_TEMP(i,j,k))
                  end do

! Set i to m // Last node in line.
                  i = m
                  
! North-East Edge
                  E_MAT(m,n,k) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(9,4) &
                              *T(i,j+1,k)+Q(9,7)*T(i,j,k+1) &
                              +Q(9,8)*T(i,j,k-1)+Q(9,6) &
                              *T(i,j-1,k)+Q(9,3)*T(i-1,j,k)) &
                              /Q(9,9)-T_TEMP(i,j,k))

! ---------------------------------------------------------------------
! End >> North Plane Sweep (Multiple Z Planes)
! ---------------------------------------------------------------------

! End k-loop 
            end do
      
! *********************************************************************
! End >> Interior Z-Plane Sweep ||| Start >> Top Plane Sweep
! *********************************************************************
            
! ---------------------------------------------------------------------
! Start >> Top South Edge Sweep (Z Max)
! ---------------------------------------------------------------------
       
! Initialize i, j, and k
            i = 1
            j = 1
            k = l
          
! Top-West-South Corner
            E_MAT(1,1,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(20,5) &
			          *T(i+1,j,k)+Q(20,4)*T(i,j+1,k)+Q(20,8) &
			          *T(i,j,k-1))/Q(20,9)-T_TEMP(i,j,k))

! Top South Edge
            do i = 2,m-1
                  E_MAT(i,1,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(15,5) &
                              *T(i+1,j,k)+Q(15,4)*T(i,j+1,k)+Q(15,8) &
                              *T(i,j,k-1)+Q(15,3)*T(i-1,j,k)) &
                              /Q(15,9)-T_TEMP(i,j,k))
            end do
            
! Set i to m // Last node in line.
            i = m
            
! Top-South-East
            E_MAT(m,j,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(24,3) &
                        *T(i-1,j,k)+Q(24,4)*T(i,j+1,k)+Q(24,8) &
                        *T(i,j,k-1))/Q(24,9)-T_TEMP(i,j,k))

! ---------------------------------------------------------------------
! End >> Top South Edge Sweep (Z Max)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> Top Plane Sweep (Z Max)
! ---------------------------------------------------------------------
                      
! Start j loop (South to North)
            do j = 2,n-1
                  
! Initialize i
                  i = 1
                  
! Top-West Edge
                  E_MAT(i,j,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(12,5) &
                              * T(i+1,j,k)+Q(12,4)*T(i,j+1,k)+Q(12,8) &
                              *T(i,j,k-1)+Q(12,6)*T(i,j-1,k))/Q(12,9) &
                              -T_TEMP(i,j,k))

! Top Plane       
                  do i = 2,m-1
                        E_MAT(i,j,l) = T(i,j,k)-T_TEMP(i,j,k)-r &
                                    *((Q(6,3)*T(i-1,j,k)+Q(6,4) &
                                    *T(i,j+1,k)+Q(6,8)*T(i,j,k-1) &
                                    +Q(6,5)*T(i+1,j,k))/Q(6,9) &
                                    -T_TEMP(i,j,k))
                  end do

! Set i to m // Last node in line.
                  i = m
                  
! Top East Edge 
                  E_MAT(m,j,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(14,3) &
                              *T(i-1,j,k)+Q(14,4)*T(i,j+1,k)+Q(14,8) &
                              *T(i,j,k-1)+Q(14,6)*T(i,j-1,k))/Q(14,9) &
                              -T_TEMP(i,j,k))
                              
! End j loop (South to North)
            end do
            
! ---------------------------------------------------------------------
! End >> Top Plane Sweep (Z Max)
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Start >> Top North Edge Sweep (Z Max)
! ---------------------------------------------------------------------

! Initialize i, and j
            i = 1
            j = n

! Top-West-North Corner
            E_MAT(1,n,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(20,5) &
                        *T(i+1,j,k)+Q(20,6)*T(i,j-1,k)+Q(20,8) &
                        *T(i,j,k-1))/Q(20,9)-T_TEMP(i,j,k))

! Top North Edge 
            do i = 2,m-1
                  E_MAT(i,n,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(13,5) &
                              *T(i+1,j,k)+Q(13,6)*T(i,j-1,k)+Q(13,8) &
                              *T(i,j,k-1)+Q(13,3)*T(i-1,j,k))/Q(13,9) &
                              -T_TEMP(i,j,k))
            end do

! Set i to m // Last node in line.
            i = m
            
! Top-North-East Corner 
            E_MAT(m,n,l) = T(i,j,k)-T_TEMP(i,j,k)-r*((Q(22,3) &
                        *T(i-1,j,k)+Q(22,6)*T(i,j-1,k)+Q(22,8) &
                        *T(i,j,k-1))/Q(22,9)-T_TEMP(i,j,k))
                        
! ---------------------------------------------------------------------
! Start >> Top North Edge Sweep (Z Max)
! ---------------------------------------------------------------------

! *********************************************************************
! End >> Top Plane Sweep
! *********************************************************************
                        
! Error Update is complete now...
          
! FIND MAX E_MAT
            E_MAX = SUM(ABS(E_MAT))/(m*n*l)
		  
! Print E_MAT MAX / Print Threshold
            print *, 'TDMA Error:', E_MAX, 'Threshold:', E, 'Run:', itr
            itr = itr + 1
          
! Now that we are finished creating the new error matrix, T(1:m,j,k) = _TEMP to T.
            T_TEMP = T
          
      end do
!     Write Results of E to file
      open(unit=3, file='error_large.txt')
      
      do k = 1,l
            do j = 1,n
                  do i = 1,m
                        write(3,*), E_MAT(i,j,k)
                  end do
            end do
      end do
      
      close(3)
      return
      end subroutine TDMA3D
                      
!	-----------------------------------------------------------------
!	-----------------------------------------------------------------
!	-----------------------------------------------------------------
!	-----------------------------------------------------------------      
!     This is the 1D Line Solution of the TDMA
      subroutine tdma_solver (a, b, c, d, n, x)
    
      implicit none
      
      integer, intent(in) :: n
      real(8), intent(in), dimension(n) :: a, c
      real(8), intent(inout), dimension(n) :: b, d
      real(8), intent(out) :: x(n)
      
      integer :: i
      real(8) :: q
      
      do i = 2,n
      q = a(i)/b(i - 1)
      b(i) = b(i) - c(i - 1)*q
      d(i) = d(i) - d(i - 1)*q
      end do
      
      q = d(n)/b(n)
      x(n) = q
      
      do i = n - 1,1,-1
      q = (d(i) - c(i)*q)/b(i)
      x(i) = q
      end do
            
      return
      end subroutine tdma_solver
      

      subroutine create_Q (Q, dx, dy, dz, k, BC_TYPE, BC)
!	-----------------------------------------------------------------


!	Define Parameters
	real(8), dimension(27,9), intent(out) :: Q
	real(8), dimension(6), intent(in) :: BC
	real(8), intent(in) :: dx, dy, dz, k 
	real(8) :: Aw, An, Ae, As, At, Ab

	character(len=4), dimension(6), intent(in) :: BC_TYPE

!	Set Matrix to 0 for initialization
	Q = 0

!	Create Intermidiate Constants
	Aw = dy*dz/dx
	An = dx*dz/dy
	Ae = dy*dz/dx
	As = dx*dz/dy
	At = dx*dy/dz
	Ab = dx*dy/dz

!     Interior Nodes = 1 // Index = 1
	Q(1,1) = 0
	Q(1,2) = 0
	Q(1,3) = k*Aw/dx
	Q(1,4) = k*An/dy
	Q(1,5) = k*Ae/dx
	Q(1,6) = k*As/dy
	Q(1,7) = k*At/dz
	Q(1,8) = k*Ab/dz
	Q(1,9) = Q(1,3)+Q(1,4)+Q(1,5)+Q(1,6)+Q(1,7)+Q(1,8)-Q(1,2)

!	Side Nodes = 6
!	West // Index = 2
      
	if (BC_TYPE(1) .eq. 'COND') then
          Q(2,1) = Aw*BC(1)
		Q(2,2) = 0
	elseif (BC_TYPE(1) .eq. 'TEMP') then
		Q(2,1) = 2*k*Aw/dx*BC(1)
		Q(2,2) = -2*k*Aw/dx
	else 
		Q(2,1) = 0
		Q(2,2) = 0
	endif
	Q(2,3) = 0
	Q(2,4) = k*An/dy
	Q(2,5) = k*Ae/dx
	Q(2,6) = k*As/dy
	Q(2,7) = k*At/dz
	Q(2,8) = k*Ab/dz
	Q(2,9) = Q(2,3)+Q(2,4)+Q(2,5)+Q(2,6)+Q(2,7)+Q(2,8)-Q(2,2)

!	North // Index = 3

	if (BC_TYPE(2) .eq. 'COND') then
		Q(3,1) = An*BC(2)
		Q(3,2) = 0
	elseif (BC_TYPE(2) .eq. 'TEMP') then
		Q(3,1) = 2*k*An/dy*BC(2)
		Q(3,2) = -2*k*An/dy
	else 
		Q(3,1) = 0
		Q(3,2) = 0
	endif

	Q(3,3) = k*Aw/dx
	Q(3,4) = 0
	Q(3,5) = k*Ae/dx
	Q(3,6) = k*As/dy
	Q(3,7) = k*At/dz
	Q(3,8) = k*Ab/dz
	Q(3,9) = Q(3,3)+Q(3,4)+Q(3,5)+Q(3,6)+Q(3,7)+Q(3,8)-Q(3,2)

!	East // Index = 4

	if (BC_TYPE(3) .eq. 'COND') then
		Q(4,1) = Ae*BC(3)
		Q(4,2) = 0
	elseif (BC_TYPE(3) .eq. 'TEMP') then
		Q(4,1) = 2*k*Ae/dx*BC(3)
		Q(4,2) = -2*k*Ae/dx
	else 
		Q(4,1) = 0
		Q(4,2) = 0
	endif

	Q(4,3) = k*Aw/dx
	Q(4,4) = k*An/dy
	Q(4,5) = 0
	Q(4,6) = k*As/dy
	Q(4,7) = k*At/dz
	Q(4,8) = k*Ab/dz
	Q(4,9) = Q(4,3)+Q(4,4)+Q(4,5)+Q(4,6)+Q(4,7)+Q(4,8)-Q(4,2)

!	South // Index = 5

	if (BC_TYPE(4) .eq. 'COND') then
		Q(5,1) = As*BC(4)
		Q(5,2) = 0
	elseif (BC_TYPE(4) .eq. 'TEMP') then
		Q(5,1) = 2*k*As/dy*BC(4)
		Q(5,2) = -2*k*As/dy
	else 
		Q(5,1) = 0
		Q(5,2) = 0
	endif

	Q(5,3) = k*Aw/dx
	Q(5,4) = k*An/dy
	Q(5,5) = k*Ae/dx
	Q(5,6) = 0
	Q(5,7) = k*At/dz
	Q(5,8) = k*Ab/dz
	Q(5,9) = Q(5,3)+Q(5,4)+Q(5,5)+Q(5,6)+Q(5,7)+Q(5,8)-Q(5,2)

!	Top // Index = 6

	if (BC_TYPE(5) .eq. 'COND') then
		Q(6,1) = At*BC(5)
		Q(6,2) = 0
	elseif (BC_TYPE(5) .eq. 'TEMP') then
		Q(6,1) = 2*k*At/dz*BC(5)
		Q(6,2) = -2*k*At/dz
	else 
		Q(6,1) = 0
		Q(6,2) = 0
	endif

	Q(6,3) = k*Aw/dx
	Q(6,4) = k*An/dy
	Q(6,5) = k*Ae/dx
	Q(6,6) = k*As/dy
	Q(6,7) = 0
	Q(6,8) = k*Ab/dz
	Q(6,9) = Q(6,3)+Q(6,4)+Q(6,5)+Q(6,6)+Q(6,7)+Q(6,8)-Q(6,2)

!	Bottom // Index = 7

	if (BC_TYPE(6) .eq. 'COND') then
		Q(7,1) = Ab*BC(6)
		Q(7,2) = 0
	elseif (BC_TYPE(6) .eq. 'TEMP') then
		Q(7,1) = 2*k*Ab/dz*BC(6)
		Q(7,2) = -2*k*Ab/dz
	else 
		Q(7,1) = 0
		Q(7,2) = 0
	endif

	Q(7,3) = k*Aw/dx
	Q(7,4) = k*An/dy
	Q(7,5) = k*Ae/dx
	Q(7,6) = k*As/dy
	Q(7,7) = k*At/dz
	Q(7,8) = 0
	Q(7,9) = Q(7,3)+Q(7,4)+Q(7,5)+Q(7,6)+Q(7,7)+Q(7,8)-Q(7,2)

!	Edge Nodes = 12
!	West-North // Index = 8
	Q(8,1) = Q(2,1)+Q(3,1)
	Q(8,2) = Q(2,2)+Q(3,2)
	Q(8,3) = 0
	Q(8,4) = 0
	Q(8,5) = k*Ae/dx
	Q(8,6) = k*As/dy
	Q(8,7) = k*At/dz
	Q(8,8) = k*Ab/dz
	Q(8,9) = Q(8,3)+Q(8,4)+Q(8,5)+Q(8,6)+Q(8,7)+Q(8,8)-Q(8,2)

!	North-East // Index = 9
	Q(9,1) = Q(3,1)+Q(4,1)
	Q(9,2) = Q(3,2)+Q(4,2)
	Q(9,3) = k*Aw/dx
	Q(9,4) = 0
	Q(9,5) = 0
	Q(9,6) = k*As/dy
	Q(9,7) = k*At/dz
	Q(9,8) = k*Ab/dz
	Q(9,9) = Q(9,3)+Q(9,4)+Q(9,5)+Q(9,6)+Q(9,7)+Q(9,8)-Q(9,2)

!	East-South // Index = 10
	Q(10,1) = Q(4,1)+Q(5,1)
	Q(10,2) = Q(4,2)+Q(5,2)
	Q(10,3) = k*Aw/dx
	Q(10,4) = k*An/dy
	Q(10,5) = 0
	Q(10,6) = 0
	Q(10,7) = k*At/dz
	Q(10,8) = k*Ab/dz
	Q(10,9) = Q(10,3)+Q(10,4)+Q(10,5)+Q(10,6)+Q(10,7)+Q(10,8)-Q(10,2)

!	South-West // Index = 11
	Q(11,1) = Q(2,1)+Q(5,1)
	Q(11,2) = Q(2,2)+Q(5,2)
	Q(11,3) = 0
	Q(11,4) = k*An/dy
	Q(11,5) = k*Ae/dx
	Q(11,6) = 0
	Q(11,7) = k*At/dz
	Q(11,8) = k*Ab/dz
	Q(11,9) = Q(11,3)+Q(11,4)+Q(11,5)+Q(11,6)+Q(11,7)+Q(11,8)-Q(11,2)

!	Top-West // Index = 12
	Q(12,1) = Q(2,1)+Q(6,1)
	Q(12,2) = Q(2,2)+Q(6,2)
	Q(12,3) = 0
	Q(12,4) = k*An/dy
	Q(12,5) = k*Ae/dx
	Q(12,6) = k*As/dy
	Q(12,7) = 0
	Q(12,8) = k*Ab/dz
	Q(12,9) = Q(12,3)+Q(12,4)+Q(12,5)+Q(12,6)+Q(12,7)+Q(12,8)-Q(12,2)

!	Top-North // Index = 13
	Q(13,1) = Q(3,1)+Q(6,1)
	Q(13,2) = Q(3,2)+Q(6,2)
	Q(13,3) = k*Aw/dx
	Q(13,4) = 0
	Q(13,5) = k*Ae/dx
	Q(13,6) = k*As/dy
	Q(13,7) = 0
	Q(13,8) = k*Ab/dz
	Q(13,9) = Q(13,3)+Q(13,4)+Q(13,5)+Q(13,6)+Q(13,7)+Q(13,8)-Q(13,2)

!	Top-East // Index = 14
	Q(14,1) = Q(4,1)+Q(6,1)
	Q(14,2) = Q(4,2)+Q(6,2)
	Q(14,3) = k*Aw/dx
	Q(14,4) = k*An/dy
	Q(14,5) = 0
	Q(14,6) = k*As/dy
	Q(14,7) = 0
	Q(14,8) = k*Ab/dz
	Q(14,9) = Q(14,3)+Q(14,4)+Q(14,5)+Q(14,6)+Q(14,7)+Q(14,8)-Q(14,2)

!	Top-South // Index = 15
	Q(15,1) = Q(5,1)+Q(6,1)
	Q(15,2) = Q(5,2)+Q(6,2)
	Q(15,3) = k*Aw/dx
	Q(15,4) = k*An/dy
	Q(15,5) = k*Ae/dx
	Q(15,6) = 0
	Q(15,7) = 0
	Q(15,8) = k*Ab/dz
	Q(15,9) = Q(15,3)+Q(15,4)+Q(15,5)+Q(15,6)+Q(15,7)+Q(15,8)-Q(15,2)

!	Bottom-West // Index = 16
	Q(16,1) = Q(2,1)+Q(7,1)
	Q(16,2) = Q(2,2)+Q(7,2)
	Q(16,3) = 0
	Q(16,4) = k*An/dy
	Q(16,5) = k*Ae/dx
	Q(16,6) = k*As/dy
	Q(16,7) = k*At/dz
	Q(16,8) = 0
	Q(16,9) = Q(16,3)+Q(16,4)+Q(16,5)+Q(16,6)+Q(16,7)+Q(16,8)-Q(16,2)

!	Bottom-North // Index = 17
	Q(17,1) = Q(3,1)+Q(7,1)
	Q(17,2) = Q(3,2)+Q(7,2)
	Q(17,3) = k*Aw/dx
	Q(17,4) = 0
	Q(17,5) = k*Ae/dx
	Q(17,6) = k*As/dy
	Q(17,7) = k*At/dz
	Q(17,8) = 0
	Q(17,9) = Q(17,3)+Q(17,4)+Q(17,5)+Q(17,6)+Q(17,7)+Q(17,8)-Q(17,2)

!	Bottom-East // Index = 18
	Q(18,1) = Q(4,1)+Q(7,1)
	Q(18,2) = Q(4,2)+Q(7,2)
	Q(18,3) = k*Aw/dx
	Q(18,4) = k*An/dy
	Q(18,5) = 0
	Q(18,6) = k*As/dy
	Q(18,7) = k*At/dz
	Q(18,8) = 0
	Q(18,9) = Q(18,3)+Q(18,4)+Q(18,5)+Q(18,6)+Q(18,7)+Q(18,8)-Q(18,2)

!	Bottom-South // Index = 19
	Q(19,1) = Q(5,1)+Q(7,1)
	Q(19,2) = Q(5,2)+Q(7,2)
	Q(19,3) = k*Aw/dx
	Q(19,4) = k*An/dy
	Q(19,5) = k*Ae/dx
	Q(19,6) = 0
	Q(19,7) = k*At/dz
	Q(19,8) = 0
	Q(19,9) = Q(19,3)+Q(19,4)+Q(19,5)+Q(19,6)+Q(19,7)+Q(19,8)-Q(19,2)

!	Corner Nodes = 8
!	West-North-Top // Index = 20
	Q(20,1) = Q(2,1)+Q(3,1)+Q(6,1)
	Q(20,2) = Q(2,2)+Q(3,2)+Q(6,2)
	Q(20,3) = 0
	Q(20,4) = 0
	Q(20,5) = k*Ae/dx
	Q(20,6) = k*As/dy
	Q(20,7) = 0
	Q(20,8) = k*Ab/dz
	Q(20,9) = Q(20,3)+Q(20,4)+Q(20,5)+Q(20,6)+Q(20,7)+Q(20,8)-Q(20,2)

!	West-North-Bottom // Index = 21
	Q(21,1) = Q(2,1)+Q(3,1)+Q(7,1)
	Q(21,2) = Q(2,2)+Q(3,2)+Q(7,2)
	Q(21,3) = 0
	Q(21,4) = 0
	Q(21,5) = k*Ae/dx
	Q(21,6) = k*As/dy
	Q(21,7) = k*At/dz
	Q(21,8) = 0
	Q(21,9) = Q(21,3)+Q(21,4)+Q(21,5)+Q(21,6)+Q(21,7)+Q(21,8)-Q(21,2)

!	North-East-Top // Index = 22
	Q(22,1) = Q(3,1)+Q(4,1)+Q(6,1)
	Q(22,2) = Q(3,2)+Q(4,2)+Q(6,2)
	Q(22,3) = k*Aw/dx
	Q(22,4) = 0
	Q(22,5) = 0
	Q(22,6) = k*As/dy
	Q(22,7) = 0
	Q(22,8) = k*Ab/dz
	Q(22,9) = Q(22,3)+Q(22,4)+Q(22,5)+Q(22,6)+Q(22,7)+Q(22,8)-Q(22,2)

!	North-East-Bottom // Index = 23
	Q(23,1) = Q(3,1)+Q(4,1)+Q(7,1)
	Q(23,2) = Q(3,2)+Q(4,2)+Q(7,2)
	Q(23,3) = k*Aw/dx
	Q(23,4) = 0
	Q(23,5) = 0
	Q(23,6) = k*As/dy
	Q(23,7) = k*At/dz
	Q(23,8) = 0
	Q(23,9) = Q(23,3)+Q(23,4)+Q(23,5)+Q(23,6)+Q(23,7)+Q(23,8)-Q(23,2)

!	East-South-Top // Index = 24
	Q(24,1) = Q(4,1)+Q(5,1)+Q(6,1)
	Q(24,2) = Q(4,2)+Q(5,2)+Q(6,2)
	Q(24,3) = k*Aw/dx
	Q(24,4) = k*An/dy
	Q(24,5) = 0
	Q(24,6) = 0
	Q(24,7) = 0
	Q(24,8) = k*Ab/dz
	Q(24,9) = Q(24,3)+Q(24,4)+Q(24,5)+Q(24,6)+Q(24,7)+Q(24,8)-Q(24,2)

!	East-South-Bottom // Index = 25
	Q(25,1) = Q(4,1)+Q(5,1)+Q(7,1)
	Q(25,2) = Q(4,2)+Q(5,2)+Q(7,2)
	Q(25,3) = k*Aw/dx
	Q(25,4) = k*An/dy
	Q(25,5) = 0
	Q(25,6) = 0
	Q(25,7) = k*At/dz
	Q(25,8) = 0
	Q(25,9) = Q(25,3)+Q(25,4)+Q(25,5)+Q(25,6)+Q(25,7)+Q(25,8)-Q(25,2)

!	South-West-Top // Index = 26
	Q(26,1) = Q(2,1)+Q(5,1)+Q(6,1)
	Q(26,2) = Q(2,2)+Q(5,2)+Q(6,2)
	Q(26,3) = 0
	Q(26,4) = k*An/dy
	Q(26,5) = k*Ae/dx
	Q(26,6) = 0
	Q(26,7) = 0
	Q(26,8) = k*Ab/dz
	Q(26,9) = Q(26,3)+Q(26,4)+Q(26,5)+Q(26,6)+Q(26,7)+Q(26,8)-Q(26,2)

!	South-West-Bottom // Index = 27
	Q(27,1) = Q(2,1)+Q(5,1)+Q(7,1)
	Q(27,2) = Q(2,2)+Q(5,2)+Q(7,2)
	Q(27,3) = 0
	Q(27,4) = k*An/dy
	Q(27,5) = k*Ae/dx
	Q(27,6) = 0
	Q(27,7) = k*At/dz
	Q(27,8) = 0
	Q(27,9) = Q(27,3)+Q(27,4)+Q(27,5)+Q(27,6)+Q(27,7)+Q(27,8)-Q(27,2)
      
!	END CREATE Q
!	-----------------------------------------------------------------

      return
	end subroutine create_Q
