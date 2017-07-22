      program mpi_tdma
      
      implicit none

      include 'include/mpif.h'
      integer comm, rank, numproc, ierror, k, n           
      real, dimension(10000,10000) :: A1,A2,A3
      real, dimension(10000) :: B1,B2,B3,x1,x2,x3,x1_t,x2_t,x3_t

      real :: t1,t2
      
      n = 10000
      call cpu_time(t1)
      print *,"CPU Time Start:",t1         
      
      
      call MPI_INIT(ierror)

      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

        if(rank == 0) then
c           This is where the first TDMA starts.
c           Read in A1 Matrix Values        
            open(12, file="data/A1_matrix.txt")
            read(12,*) A1
            close(12)
c            A1 = transpose(A1)
                         
c           Read in B1 Matrix Values        
            open(12, file="data/B1_matrix.txt")
            read(12,*) B1
            close(12)
            
c           Read in x1 Matrix Values as Temp        
            open(12, file="data/x1_matrix.txt")
            read(12,*) x1_t
            close(12)

            call tdm_solve(A1,B1,x1,n)
            call cpu_time(t2)
            print *,"CPU Time End A1:",t2
c           This is where the first TDMA ends.

c           This is where we check x1 from Fortran with x1_t from MATLAB
c            print *,"This is x1 ",x1
c            print *,"This is x1_t ",x1_t
c           If true, print success and the answer!       
        endif
        
        if(rank == 1) then
c           This is where the second TDMA starts.
c           Read in A2 Matrix Values        
            open(12, file="data/A2_matrix.txt")
            read(12,*) A2
            close(12)
c            A2 = transpose(A2) 
            
c           Read in B2 Matrix Values        
            open(12, file="data/B2_matrix.txt")
            read(12,*) B2
            close(12)
            
c           Read in x2 Matrix Values as Temp        
            open(12, file="data/x2_matrix.txt")
            read(12,*) x2_t
            close(12)

            call tdm_solve(A2,B2,x2,n)
            call cpu_time(t2)
            print *,"CPU Time End A2:",t2
c           This is where the second TDMA ends. 

c           This is where we check x2 from Fortran with x2_t from MATLAB
c            print *,"This is x2 ",x2
c            print *,"This is x2_t ",x2_t
c           If true, print success and the answer! 
            
        endif
        
        if(rank == 2) then
c           This is where the TDMA starts.
c           Read in A3 Matrix Values        
            open(12, file="data/A3_matrix.txt")
            read(12,*) A3
            close(12)
c            A3 = transpose(A3) 
            
c           Read in B3 Matrix Values        
            open(12, file="data/B3_matrix.txt")
            read(12,*) B3
            close(12)
            
c           Read in x3 Matrix Values as Temp        
            open(12, file="data/x3_matrix.txt")
            read(12,*) x3_t
            close(12)

            call tdm_solve(A3,B3,x3,n)
            call cpu_time(t2)
            print *,"CPU Time End A3:",t2
c           This is where the TDMA ends. 

c           This is where we check x3 from Fortran with x3_t from MATLAB
c            print *,"This is x3 ",x3
c            print *,"This is x3_t ",x3_t
c           If true, print success and the answer! 

        endif   

      call MPI_FINALIZE(ierror)
      
      call cpu_time(t2)
      print *,"The program took this many seconds ",t2-t1

      end program mpi_tdma
      
      subroutine tdm_solve(A,B,x,n)
      
        implicit none
        
        integer,intent(in) :: n
        real, dimension(n,n) :: A
        real, dimension(n), intent(in) :: B
        real, dimension(n), intent(out) :: x
        real, dimension(n) :: at,bt,ct,dt,cp,dp
        
        real(8) :: m
        integer :: i
        
        do i = 1,n-1
            at(i) = A(i,i)
            bt(i) = A(i,i+1)
            ct(i) = A(i+1,i)
            dt(i) = B(i)
        end do
        
        at(n) = A(n,n)
        dt(n) = B(n)
        
        ! initialize c-prime and d-prime
        cp(1) = ct(1)/bt(1)
        dp(1) = dt(1)/bt(1)
        ! solve for vectors c-prime and d-prime
        do i = 2,n
           m = bt(i)-cp(i-1)*at(i)
           cp(i) = ct(i)/m
           dp(i) = (dt(i)-dp(i-1)*at(i))/m
        end do
        ! initialize x
         x(n) = dp(n)
         
        ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
        
        
      end subroutine tdm_solve