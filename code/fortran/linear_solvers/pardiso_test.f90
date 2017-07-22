!  pardiso_test.f90 
!
!  FUNCTIONS:
!  pardiso_test - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: pardiso_test
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
    !include 'mkl_pardiso.f90'
    
    program pardiso_test
        
    implicit none
    
    include 'mkl.fi'
    
    real(8), dimension(5,5) :: A
    real(8), dimension(5) :: b, x
    real(8), dimension(13) :: v, c, r
    
    integer :: i, j, m, n, info
    integer, dimension(8) :: job
    integer, dimension(13) :: ja, ia
    real(8), dimension(13) :: Acsr
    
    m = 5
    n = 5
    
    i = 1
    j = 1
    
    A(i,j) = 1
    A(i,j+1) = 1
    
    do i = 2, 4
          do j = 2, 4
                A(i,j-1) = -1
                A(i,j) = 1
                A(i,j+1) = -1
          end do
    end do 
    
    i = 5
    j = 5
    
    A(i,j) = 1
    A(i,j-1) = -1
    b = 2
    
    job(1) = 0
    job(2) = 1
    job(3) = 0
    job(4) = 2
    job(5) = 25
    job(6) = 1
    Acsr = 0
    ja = 0
    ia = 0
    
!   This function converts dense matrix to CSR format.
    call mkl_ddnscsr(job,m,n,A,m,Acsr,ja,ia,info)
    
    print *, 'Acsr:', Acsr
    
    ! Variables

    ! Body of pardiso_test
    print *, 'Hello World'

    end program pardiso_test
