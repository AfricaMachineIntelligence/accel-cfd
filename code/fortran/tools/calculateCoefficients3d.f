! This program is written in free-form

subroutine calculateCoefficients3d(m, n, l, c, r)

    integer, intent(in) :: m, n, l
    integer, intent(out) :: c, r

    integer :: nodes_c, nodes_l, nodes_w, nodes_i
    integer :: coeff_c, coeff_l, coeff_w, coeff_i

    implicit none

!   Calculate the number of nodes in a 3D grid
    nodes_c = 8
    nodes_l = (m - 2) * 4 + (n - 2) * 4 + (l - 2) * 4
    nodes_w = ((m - 2)*(n - 2)) * 2 + ((n - 2)*(l - 2)) * 2 + ((l - 2)*(m - 2)) * 2
    nodes_i = (m - 2)*(n - 2)*(l - 2)
      
!   Calculate the number of coefficients in the for each node type
    coeff_c = nodes_c * 4
    coeff_l = nodes_l * 5
    coeff_w = nodes_w * 6
    coeff_i = nodes_i * 7

!   Calculate total coefficients and rows
    c = coeff_c + coeff_l + coeff_w + coeff_i
    r = m * n * l

    return

end subroutine