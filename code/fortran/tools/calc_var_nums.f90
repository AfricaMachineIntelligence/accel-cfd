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