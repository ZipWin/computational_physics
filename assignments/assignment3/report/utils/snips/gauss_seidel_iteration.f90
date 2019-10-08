subroutine gauss_seidel_iteration(co_matrix, ndim, solution, omega)
    implicit none
    integer(8), intent(in) :: ndim
    real(8), optional, intent(in) :: omega
    real(8), dimension(ndim, ndim+1), intent(in) :: co_matrix
    type(result), dimension(ndim), intent(out) :: solution

    real(8), dimension(ndim, ndim+1) :: A
    real(8) :: lambda = 1.0d0
    real(8), dimension(ndim) :: X, X_tmp, Error
    integer(8) :: iter, i

    A = co_matrix
    call pivot_exchange(A, ndim, ndim+1)

    if (present(omega)) then
        lambda = omega
    end if

    X = 0.0d0
    do iter = 1, 10000
        X_tmp = X
        do i = 1, ndim
            X(i) = -(dot_product(A(i, :ndim), X) - A(i, i) * X(i) - A(i, ndim+1)) / A(i, i)
        end do
        X = (1 - lambda) * X_tmp + lambda * X
        Error = abs(X - X_tmp)

        if (maxval(Error) < EPSILON) then
            exit
        end if

    end do
    print "('Iter: ',i4)", iter
    solution%value = X
    solution%error = Error

    return
end subroutine