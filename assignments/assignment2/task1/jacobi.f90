program main
    implicit none
    integer, parameter :: MAX_ROOTS_NUM = 100
    real(8), external :: phi1, phi2
    real(8), parameter :: step_size = 0.1d0, epsilon1 = 1e-5, epsilon2 = 1e-5
    real(8), dimension(MAX_ROOTS_NUM) :: roots, errors
    integer, dimension(MAX_ROOTS_NUM) :: iters
    real(8) :: left, right
    character(50) :: file_name = "history.dat"
    integer :: iter, root_num, i
    logical :: stat

    left = -10.0d0
    right = 10.0d0

    print "(a46)", "The form of iteration function phi: (3x)^{1/3}"
    call multi_jacobi(phi1, left, right, step_size, epsilon1, epsilon2, file_name, root_num, roots, errors, iters)
    print "(a35, i3)", "The number of roots that is found: ", root_num
    do i = 1, root_num
        print *, " "
        print "(a8, f16.6)", "root: ", roots(i)
        print "(a8, f16.6)", "error: ", errors(i)
        print "(a8, i16)", "iter: ", iters(i)
    end do

    print *, " "
    print "(a41)", "The form of iteration function phi: x^3/3"
    call multi_jacobi(phi2, left, right, step_size, epsilon1, epsilon2, file_name, root_num, roots, errors, iters)
    print "(a35, i3)", "The number of roots that is found: ", root_num
    do i = 1, root_num
        print *, " "
        print "(a8, f16.6)", "root: ", roots(i)
        print "(a8, f16.6)", "error: ", errors(i)
        print "(a8, i16)", "iter: ", iters(i)
    end do

end program


!------------------------------------------------------------------------------
! This is the function for iteration x = phi(x)
! To find different roots, different phi should be used
! The form of function phi is chosen manually
!
! Argments:
!          x: a real(8) number
!------------------------------------------------------------------------------
real(8) function phi1(x) result(phi)
    implicit none
    real(8), intent(in) :: x

    phi = sign(abs(3*x)**(1.0d0/3.0d0), x)

    return
end function

real(8) function phi2(x) result(phi)
    implicit none
    real(8), intent(in) :: x

    phi = x**3 / 3

    return
end function


real(8) function numerical_derivative(func, x) result(dfdx)
!------------------------------------------------------------------------------
! This function computes the numerical derivative of func at x. 
! It returns the results in real(8).
!
! Argments:
!          func: a real(8) function
!          x: a real(8) number
!------------------------------------------------------------------------------
    implicit none
    real(8), external :: func
    real(8), intent(in) :: x

    real(8) :: h
    h = 1e-5

    dfdx = (func(x+h) - func(x-h)) / (2*h)

    return
end function


subroutine convergence_check(phi, a, b, stat)
!------------------------------------------------------------------------------
! This subroutine check whether the iteration x = phi(x) converges in [a, b].
! It return the results in stat. If converges, stat is true; if not, false.
!
! Argments:
!          phi: a real(8) function
!          a: real(8), the left boundary of the interval
!          b: real(8), the right boundary of the interval
!          stat: logical, save the state of convergence
!------------------------------------------------------------------------------
    implicit none
    real(8), external :: phi, numerical_derivative
    real(8), intent(in) :: a, b
    logical, intent(out) :: stat

    integer(8), parameter :: EVAL_NUM = 1000
    real(8), dimension(EVAL_NUM) :: xs, phis, dphis
    real(8) :: eval_step
    integer :: i

    eval_step = abs(a - b) / EVAL_NUM

    do i = 1, EVAL_NUM
        xs(i) = a + eval_step * i
        phis(i) = phi(xs(i))
        dphis(i) = numerical_derivative(phi, xs(i))
    end do

    stat = .false.
    if (maxval(phis) <= maxval((/a, b/)) .and. minval(phis) >= minval((/a, b/))) then
        if (maxval(abs(dphis)) < 1) then
            stat = .true.
        end if
    end if

    return
end subroutine


subroutine jacobi(phi, x_init, epsilon1, epsilon2, file_name, root, error, iter, stat)
!------------------------------------------------------------------------------
! This subroutine finds a root for a given Jacobi iteration x = phi(x).
! 
! Argments:
!          phi: a real(8) function
!          x_init: real(8), the initial point
!          epsilon1: real(8), the expected error between x_n and x_{n-1}
!          epsilon2: real(8), the expected error between f(x_n) and 0
!          file_name: a string, the iteration history will be save into a file named as file_name
!          root: save the found root
!          error: save the error
!          iter: save the number of times of iteration
!          stat: if converges, true; else, false
!------------------------------------------------------------------------------
    implicit none
    real(8), external :: phi
    real(8), intent(in) :: x_init, epsilon1, epsilon2
    character(50), intent(in) :: file_name
    real(8), intent(out) :: root, error
    integer, intent(out) :: iter
    logical, intent(out) :: stat

    integer, parameter :: MAX_ITER_NUM = 1000000
    real(8), dimension(MAX_ITER_NUM) :: history
    real(8) :: x, f
    f(x) = phi(x) - x
    integer :: i
    logical :: stop_condition

    x = x_init
    i = 1
    history(i) = x

    stop_condition = .false.
    do while (.not. stop_condition)
        x = phi(x)
        i = i + 1
        history(i) = x
        stop_condition = (abs(x - history(i-1)) <= epsilon1)    &
                         .and. (abs(f(x)) <= epsilon2)          &
                         .and. (i <= max_iter_num)
    end do

    iter = i
    root = history(i)
    error = abs(history(i) - history(i-1))

    if (i == max_iter_num) then
        stat = .false.
    else
        stat = .true.
    end if

    open(file=file_name, unit=10, action="write")
    do i = 1, iter
        write (10, "(i8, a4, f8.6)") i, "   ", history(i)
    end do
    close(unit=10)

    return
end subroutine


subroutine multi_jacobi(phi, left, right, step_size, epsilon1, epsilon2, file_name, root_num, roots, errors, iters)
!------------------------------------------------------------------------------
! This subroutine search all the possible roots in [left, right] using Jacobi.
! It divides the interval into several adjoint small intervals, whose length is step_size.
! 
! Argments:
!          phi: a real(8) function
!          left: real(8), the left boundary of the interval
!          right: real(8), the right boundary of the interval
!          step_size: real(8), the small interval length
!          epsilon1: real(8), the expected error between x_n and x_{n-1}
!          epsilon2: real(8), the expected error between f(x_n) and 0
!          file_name: a string, the iteration history will be save into a file named as file_name
!          root_num: the number of roots that has been found
!          roots: save the found roots
!          errors: save the errors
!          iters: save the numbers of times of iteration
!------------------------------------------------------------------------------
    implicit none
    real(8), external :: phi
    real(8), intent(in) :: left, right, step_size, epsilon1, epsilon2
    character(50), intent(in) :: file_name
    integer, intent(out) :: root_num

    integer, parameter :: MAX_ROOTS_NUM = 100
    real(8), dimension(MAX_ROOTS_NUM) :: roots, errors
    integer, dimension(MAX_ROOTS_NUM) :: iters
    real(8) :: root, error, a, b, x_init
    integer :: iter, i
    logical :: can_converge, is_converge

    i = 1
    a = left
    b = left + step_size
    do while (b <= right)
        call convergence_check(phi, a, b, can_converge)
        if (can_converge) then
            call random_seed()
            call random_number(x_init)
            x_init = a + x_init * (a - b)
            call jacobi(phi, x_init, epsilon1, epsilon2, file_name, root, error, iter, is_converge)
            if (is_converge) then
                if (i .ne. 1) then
                    if (abs(root - roots(i-1)) >= epsilon1) then
                        roots(i) = root
                        errors(i) = error
                        iters(i) = iter
                        i = i + 1
                    end if
                else
                    roots(i) = root
                    errors(i) = error
                    iters(i) = iter
                    i = i + 1
                end if
            end if
        end if
        a = a + step_size
        b = b + step_size
    end do
    root_num = i - 1

    return
end subroutine