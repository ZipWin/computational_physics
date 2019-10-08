module utils
    implicit none
    type result
        real(8) :: value
        real(8) :: error
    end type

    contains
    subroutine pivot_exchange(matrix, ndim1, ndim2)
        implicit none
        integer(8), intent(in) :: ndim1, ndim2
        real(8), dimension(ndim1, ndim2) :: matrix

        integer :: col, i
        real(8), dimension(ndim2) :: tmp_row

        do col = 1, ndim2
            do i = col, ndim1
                if (matrix(i, col) > matrix(col, col)) then
                    tmp_row = matrix(col, :)
                    matrix(col, :) = matrix(i, :)
                    matrix(i, :) = tmp_row
                end if
            end do
        end do

        return
    end subroutine

end module


module linear_eqs_solver
    use utils
    implicit none
    real(8), parameter :: EPSILON = 1e-3

    contains
    subroutine gauss_elimination(co_matrix, ndim, solution)
        implicit none
        integer(8), intent(in) :: ndim
        real(8), dimension(ndim, ndim+1), intent(in) :: co_matrix
        type(result), dimension(ndim), intent(out) :: solution

        real(8), dimension(ndim, ndim+1) :: A
        real(8), dimension(ndim) :: tmp_row, X
        integer(8) :: t, i

        A = co_matrix
        do t = 1, ndim
            ! pivot exchange
            do i = t, ndim
                if (abs(A(i, t)) > abs(A(t, t))) then
                    tmp_row = A(t, :)
                    A(t, :) = A(i, :)
                    A(i, :) = tmp_row
                end if
            end do
            ! Gauss elimination
            do i = t+1, ndim
                A(i, :) = A(i, :) - A(i, t) / A(t, t) * A(t, :)
            end do
        end do

        print "('The Gauss-eliminated extented coefficients matrix is')"
        print "(10f7.2)", (A(i, :), i=1,9)

        X = 0.0d0
        do t = ndim, 1, -1
            X(t) = (A(t, ndim+1) - dot_product(A(t, t+1:ndim), X(t+1:ndim))) / A(t, t)
        end do

        solution%value = X
        solution%error = 0.0d0

        return
    end subroutine

    subroutine doolittle(co_matrix, ndim, solution)
        implicit none
        integer(8), intent(in) :: ndim
        real(8), dimension(ndim, ndim+1), intent(in) :: co_matrix
        type(result), dimension(ndim), intent(out) :: solution

        real(8), dimension(ndim, ndim+1) :: A
        real(8), dimension(ndim, ndim) :: L, U
        real(8), dimension(ndim) :: X, Y
        integer(8) :: i, j, k

        A = co_matrix
        L = 0.0d0; U = 0.0d0;
        
        call pivot_exchange(A, ndim, ndim+1)

        do k = 1, ndim
            do j = k, ndim
                U(k, j) = A(k, j) - dot_product(L(k, :k-1), U(:k-1, j))
            end do
            L(k, k) = 1.0d0
            do i = k+1, ndim
                L(i, k) = (A(i, k) - dot_product(L(i, :k-1), U(:k-1, k))) / U(k, k)
            end do
        end do

        do i = 1, ndim
            Y(i) = A(i, ndim+1) - dot_product(L(i, :i-1), Y(:i-1))
        end do
        do i = ndim, 1, -1
            X(i) = (Y(i) - dot_product(U(i, i+1:), X(i+1:))) / U(i, i)
        end do

        solution%value = X
        solution%error = 0.0d0

        return
    end subroutine


    subroutine jacobi_iteration(co_matrix, ndim, solution)
        implicit none
        integer(8), intent(in) :: ndim
        real(8), dimension(ndim, ndim+1), intent(in) :: co_matrix
        type(result), dimension(ndim), intent(out) :: solution

        real(8), dimension(ndim, ndim+1) :: A
        real(8), dimension(ndim) :: X, X_tmp, Error
        integer(8) :: iter, i

        A = co_matrix
        call pivot_exchange(A, ndim, ndim+1)

        X = 0.0d0
        do iter = 1, 10000
            do i = 1, ndim
                X_tmp(i) = -(dot_product(A(i, :ndim), X) - A(i, i) * X(i) - A(i, ndim+1)) / A(i, i)
            end do
            Error = abs(X - X_tmp)
            X = X_tmp

            if (maxval(Error) < EPSILON) then
                exit
            end if

        end do
        print *, iter
        solution%value = X
        solution%error = Error

        return
    end subroutine


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

end module


program main
    use utils
    use linear_eqs_solver
    implicit none
    integer(8), parameter :: ndim = 9
    real(8), dimension(ndim, ndim+1) :: A
    type(result), dimension(ndim) :: solution
    integer(8) :: i

    A(1, :9) = (/31.0d0, -13.0d0, 0.0d0, 0.0d0, 0.0d0, -10.0d0, 0.0d0, 0.0d0, 0.0d0/)
    A(2, :9) = (/-13.0d0, 35.0d0, -9.0d0, 0.0d0, -11.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
    A(3, :9) = (/0.0d0, -9.0d0, 31.0d0, -10.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
    A(4, :9) = (/0.0d0, 0.0d0, -10.0d0, 79.0d0, -30.0d0, 0.0d0, 0.0d0, 0.0d0, -9.0d0/)
    A(5, :9) = (/0.0d0, 0.0d0, 0.0d0, -30.0d0, 57.0d0, -7.0d0, 0.0d0, -5.0d0, 0.0d0/)
    A(6, :9) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, -7.0d0, 47.0d0, -30.0d0, 0.0d0, 0.0d0/)
    A(7, :9) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, -30.0d0, 41.0d0, 0.0d0, 0.0d0/)
    A(8, :9) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, -5.0d0, 0.0d0, 0.0d0, 27.0d0, -2.0d0/)
    A(9, :9) = (/0.0d0, 0.0d0, 0.0d0, -9.0d0, 0.0d0, 0.0d0, 0.0d0, -2.0d0, 29.0d0/)
    A(:, 10) = (/-15.0d0, 27.0d0, -23.0d0, 0.0d0, -20.0d0, 12.0d0, -7.0d0, 7.0d0, 10.0d0/)

    print "('The extented coefficients matrix is')"
    print "(10f7.2)", (A(i, :), i=1,9)

    call gauss_seidel_iteration(A, ndim, solution)

    print "('The solution vector is')"
    print "(9f8.4)", (solution(i)%value, i=1,9)
    print "('The errors are')"
    print "(9f8.4)", (solution(i)%error, i=1,9)

    print "('Checking result...')"
    print "(9f8.2)", (dot_product(A(i, :9), solution%value), i=1,9)

end program