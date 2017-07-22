      program helloworld

      include 'include/mpif.h'
      integer comm, rank, numproc, ierror

      call MPI_INIT(ierror)

      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

c     print *,"Hello World from Processor ",rank," of ",numproc

        if(rank.eq.0) then
            print *,"I'm processor ",rank
        endif
        
        if(rank.eq.1) then
            print *,"I'm processor ",rank
        endif
        
        if(rank.eq.2) then
            print *,"I'm processor ",rank
        endif   

      call MPI_FINALIZE(ierror)

      end program helloworld