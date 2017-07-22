      program addvec

      include 'include/mpif.h'
      integer comm, rank, numproc, ierror
      
      real, dimension(10) :: A, B, C
      
      A = (/ 1,2,3,4,5,6,7,8,9,10 /)
      B = (/ 2,4,6,8,10,12,14,16,18,20 /)

      call MPI_INIT(ierror)

      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

        if(rank.eq.0) then
            C = 2*A
            print *,"I'm processor ",rank," I multiplied A by 2",C
        endif
        
        if(rank.eq.1) then
            C = A+B
            print *,"I'm processor ",rank," I multiplied A by 2",C
        endif
        
        if(rank.eq.2) then
            C = 2*B
            print *,"I'm processor ",rank," I multiplied B by 2",C
        endif   

      call MPI_FINALIZE(ierror)

      end program addvec