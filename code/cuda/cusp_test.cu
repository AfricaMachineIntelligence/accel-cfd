// Original File from: https://cusplibrary.github.io/group__krylov__methods.html#ga23cfa8325966505d6580151f91525887
// Modified to verify proper compilation

#include <cusp/csr_matrix.h>
#include <cusp/monitor.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/gallery/poisson.h>
int main(void)
{
    // create an empty sparse matrix structure (CSR format)
    cusp::csr_matrix<int, float, cusp::device_memory> A;
    // initialize matrix
    cusp::gallery::poisson5pt(A, 10, 10);
    // allocate storage for solution (x) and right hand side (b)
    cusp::array1d<float, cusp::device_memory> x(A.num_rows, 0);
    cusp::array1d<float, cusp::device_memory> b(A.num_rows, 1);
    // set stopping criteria:
    //  iteration_limit    = 100
    //  relative_tolerance = 1e-6
    //  absolute_tolerance = 0
    //  verbose            = true
    cusp::monitor<float> monitor(b, 100, 1e-6, 0, true);
    // set preconditioner (identity)
    cusp::identity_operator<float, cusp::device_memory> M(A.num_rows, A.num_rows);
    // solve the linear system A x = b
    cusp::krylov::bicgstab(A, x, b, monitor, M);
    return 0;
}
