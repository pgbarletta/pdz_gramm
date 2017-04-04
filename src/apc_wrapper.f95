subroutine apc_wrapper(m, n, mtx_a, mtx_b, indices)
implicit none

integer*8, intent(in) :: m, n
real*8, dimension (m, n), intent(in) :: mtx_a, mtx_b
integer*8, intent(out) :: indices(n)


indices(1) = 1
indices(2) = 44
indices(n-1) = 1
indices(n) = 2



end subroutine
