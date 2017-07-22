      program main

c     Written by: Matt Blomquist
c     Date: 2016-08-09

c     This program solves a 3D Conduction Problem using a line-by-line TDMA technique to 
c     solve the disretized equations.

      implicit none


c     Define Variables Used
      real :: t_start,t_end,m_val

      integer :: x_length,y_length,z_length 
      integer :: i,j,k

c     Generate Matrix
      real, dimension(300,200,100) :: COND

c     Define Matrix Size
      x_length = 300
      y_length = 200
      z_length = 100      

c     Set Values to 0
      do i = 1,x_length
            do j = 0,y_length
                  do k = 0,z_length
                        COND(i,j,k) = 0.0
                  enddo
            enddo
      enddo      

c     Define Boundary Conditions


c     Discretize equations


c     Start Clock / Run 3D TDMA
      call cpu_time(t_start)


c     End 3D TDMA / Stop Clock
      call cpu_time(t_end)


c     Print 3D TDMA Duration
      print *,'Total TDMA Sweep duration:',t_end-t_start

      end program main