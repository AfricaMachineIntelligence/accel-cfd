       program make_matrix
          
          implicit none

          real, dimension(10,10,10) :: A
          integer i, j, k
          
          !DIR$ PARALLEL 
          do i = 1,10
            do j = 1,10
                do k = 1,10
                    A(i,j,k) = 1
                enddo
            enddo
          enddo

          print *, A

       end program make_matrix