      program mpi_tdma
      
      implicit none

      include 'include/mpif.h'
      integer comm, rank, numproc, ierror, k
      
      real, dimension(10,10) :: A1,A2,A3
      real, dimension(10) :: B1,B2,B3,x1,x2,x3,x1_t,x2_t,x3_t

      call MPI_INIT(ierror)

      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

        if(rank == 0) then
        
c           Read in A1 Matrix Values        
            open(12, file="data/A1_matrix.txt")
            read(12,*) A1
            close(12)
            
            
c           Read in B1 Matrix Values        
            open(12, file="data/B1_matrix.txt")
            read(12,*) B1
            close(12)
            
c           Read in x1 Matrix Values as Temp        
            open(12, file="data/x1_matrix.txt")
            read(12,*) x1_t
            close(12)   
            A1 = transpose(A1)     
            
c           This is where the first TDMA starts.
            call tdm_solve(A1,B1,x1)

c           This is where the first TDMA ends.

c           This is where we check x1 from Fortran with x1_t from MATLAB
            print *,"This is x1 ",x1
            print *,"This is x1_t ",x1_t
c           If true, print success and the answer!       
        endif
        
        if(rank == 1) then
c           Read in A2 Matrix Values        
            open(12, file="data/A2_matrix.txt")
            read(12,*) A2
            close(12)
            
c           Read in B2 Matrix Values        
            open(12, file="data/B2_matrix.txt")
            read(12,*) B2
            close(12)
            
c           Read in x2 Matrix Values as Temp        
            open(12, file="data/x2_matrix.txt")
            read(12,*) x2_t
            close(12)

c           This is where the second TDMA starts.
            call tdm_solve(A2,B2,x2)

c           This is where the second TDMA ends. 

c           This is where we check x2 from Fortran with x2_t from MATLAB
c            print *,"This is x2 ",x2
c            print *,"This is x2_t ",x2_t
c           If true, print success and the answer! 
            
        endif
        
        if(rank == 2) then
c           Read in A3 Matrix Values        
            open(12, file="data/A3_matrix.txt")
            read(12,*) A3
            close(12)
            
c           Read in B3 Matrix Values        
            open(12, file="data/B3_matrix.txt")
            read(12,*) B3
            close(12)
            
c           Read in x3 Matrix Values as Temp        
            open(12, file="data/x3_matrix.txt")
            read(12,*) x3_t
            close(12)

c           This is where the TDMA starts.
            call tdm_solve(A3,B3,x3)

c           This is where the TDMA ends. 

c           This is where we check x3 from Fortran with x3_t from MATLAB
c            print *,"This is x3 ",x3
c            print *,"This is x3_t ",x3_t
c           If true, print success and the answer! 

        endif   

      call MPI_FINALIZE(ierror)

      end program mpi_tdma
      
      subroutine tdm_solve(A,B,x)
      
        implicit none
        
        real, dimension(10,10) :: A
        real, dimension(10) :: B,x
        real :: c
        
        integer :: i,j
        
        A(1,2) = A(1,2)/A(1,1)
        B(1) = B(1)/A(1,1)
        
        do i = 2,10-1
            c = A(i,i)-A(i-1,i)* A(i,i-1)
            A(i+1,i)= A(i+1,i)/c;
            B(i)=(B(i)-A(i-1,i)*B(i-1))/c
        end do
            
        i = 10
        
        x(i) = (B(i)-A(i-1,i)*B(i-1))/(A(i,i)-A(i-1,i)*A(i,i-1))
        
        do i = 1,9
            x(i) = -A(10-i,10-i)*x(10-i)+B(10-i)
        end do
        
      end subroutine tdm_solve