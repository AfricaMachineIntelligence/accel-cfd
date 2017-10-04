! This program is written in free-form

subroutine solver3d()

    integer, intent(in) :: m,n,l
    integer, intent(in) :: c, r
     
    real(dp), dimension(m, n, l), intent(in) :: Ab, As, Aw, Ap, Ae, An, At
    real(dp), dimension() :: values, columns, 
    real(dp), dimension() :: rowIndex

    implicit none

!   Calculate the number of coefficients for vectors to hold
    call calculateCoefficients3d(m,n,l,c,r)

!   Populate arrays
    call solverSetup3d(m, n, l, c, r, Ab, As, Aw, Ap, Ae, An, At, values, columns, rowIndex)

!   Run solver based on set value
!   1. Intel MKL Paradiso
!   2. BiCGSTAB(ell)

!   Print Solver Start
    print(*), "System Size (m,n,l) \n Tolerance Set to: tol \n Max Iterations Set to: maxit \n Solving..."

    if (set.eq.1) then

!     Load include files
      include "mkl.fi"
      include "mkl_paradiso.f90"

!     Call Paradiso function
      call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)

    else 

!     Call BiCGSTAB Algorithm
      call bistbl (l, n, x, b, mv, solve, mxmv, work, ldw, rwork, ldrw, iwork, info)

    end if

!   Print Results and Return
    print(*), "Total Solution Time: calc_time\n Final Resisdual: min_res \n Number of Iterations: itr \n"

    return

end subroutine