module utils
    implicit none
    integer, parameter :: EVAL_NUM = 1000

    type result
        real(8) :: value
        real(8) :: error
    end type

    contains
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
        real(8), external :: phi
        real(8), intent(in) :: a, b
        logical, intent(out) :: stat

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

end module


module roots_seeker
    use utils

    implicit none
    integer, parameter :: MAX_ITER_NUM = 1000000
    integer, parameter :: MAX_ROOTS_NUM = 100

    contains
    !------------------------------------------------------------------------------
    ! This following subroutines finds only one root in a given interval or initial point
    ! 
    ! Argments:
    !          f: real(8), the original function
    !          phi: a real(8) function (not need for some methods)
    !          left_init: the lower limit of the interval
    !          right_init: the greater limit of the interval
    !          x_init: real(8), the initial point
    !          epsilon1: real(8), the expected error between x_n and x_{n-1}
    !          epsilon2: real(8), the expected error between f(x_n) and 0
    !          file_name: a string, the iteration history will be save into a file named as file_name
    !          root: save the found root
    !          error: save the error
    !          iter: save the number of times of iteration
    !          stat: if converges, true; else, false
    !------------------------------------------------------------------------------
    subroutine bisection(f, left_init, right_init, epsilon1, epsilon2, file_name, root, iter, stat)
        implicit none
        real(8), external :: f
        real(8), intent(in) :: left_init, right_init, epsilon1, epsilon2
        character(50), intent(in) :: file_name
        type(result), intent(out) :: root
        integer, intent(out) :: iter
        logical, intent(out) :: stat

        real(8), dimension(MAX_ITER_NUM) :: history
        real(8) :: left, middle, right, lv, mv, rv
        integer :: i
        logical :: stop_condition

        left = left_init
        right = right_init
        middle = (left + right) / 2
        i = 1
        history(i) = middle
        stat = .true.

        stop_condition = .false.
        do while (.not. stop_condition)
            lv = f(left)
            mv = f(middle)
            rv = f(right)

            if (abs(lv - 0.0d0) <= epsilon2) then
                right = left + 2 * epsilon1
            else if (abs(rv - 0.0d0) <= epsilon2) then
                left = right - 2 * epsilon1
            end if

            if (lv * rv < 0) then
                if (lv * mv < 0) then
                    right = middle
                else if (mv * rv < 0) then
                    left = middle
                end if
            else
                stat = .false.
                exit
            end if

            middle = (left + right) / 2
            i = i + 1
            history(i) = middle

            stop_condition = ((abs(right - left) <= epsilon1)       &
                             .and. (abs(f(middle)) <= epsilon2))    &
                             .or. (i >= max_iter_num)

        end do

        root%value = middle
        root%error = right - left
        iter = i

        if (i == MAX_ITER_NUM) then
            stat = .false.
        end if

        open(file=file_name, unit=10, action="write")
        do i = 1, iter
            write (10, "(i8, a4, f8.6)") i, "   ", history(i)
        end do
        close(unit=10)

        return
    end subroutine

    subroutine jacobi(f, phi, x_init, epsilon1, epsilon2, file_name, root, iter, stat)
        implicit none
        real(8), external :: phi, f
        real(8), intent(in) :: x_init, epsilon1, epsilon2
        character(50), intent(in) :: file_name
        type(result), intent(out) :: root
        integer, intent(out) :: iter
        logical, intent(out) :: stat

        real(8), dimension(MAX_ITER_NUM) :: history
        real(8) :: x
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
            stop_condition = ((x - history(i-1) <= epsilon1)       &
                             .and. (abs(f(x)) <= epsilon2))    &
                             .or. (i >= max_iter_num)
        end do

        iter = i
        root%value = history(i)
        root%error = abs(history(i) - history(i-1))

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

    subroutine newton_downhill(f, x_init, epsilon1, epsilon2, file_name, root, iter, stat)
        implicit none
        real(8), external :: f
        real(8), intent(in) :: x_init, epsilon1, epsilon2
        character(50), intent(in) :: file_name
        type(result), intent(out) :: root
        integer, intent(out) :: iter
        logical, intent(out) :: stat

        real(8), dimension(MAX_ITER_NUM) :: history
        real(8) :: x, tmp, lambda
        integer :: i
        logical :: stop_condition

        x = x_init
        i = 1
        history(i) = x

        stop_condition = .false.
        do while (.not. stop_condition)
            lambda = 1.0d0 ! reset the value of lambda
            do while (.true.)
                tmp = x - lambda * f(x) / numerical_derivative(f, x)
                if (abs(f(tmp)) < abs(f(x))) then
                    exit
                else
                    lambda = lambda / 2
                end if
            end do
            x = tmp
            i = i + 1
            history(i) = x
            stop_condition = ((x - history(i-1) <= epsilon1)       &
                             .and. (abs(f(x)) <= epsilon2))    &
                             .or. (i >= max_iter_num)
        end do

        iter = i
        root%value = history(i)
        root%error = abs(history(i) - history(i-1))

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

    subroutine post_accelerating(f, phi, x_init, epsilon1, epsilon2, file_name, root, iter, stat)
        implicit none
        real(8), external :: f, phi
        real(8), intent(in) :: x_init, epsilon1, epsilon2
        character(50), intent(in) :: file_name
        type(result), intent(out) :: root
        integer, intent(out) :: iter
        logical, intent(out) :: stat

        real(8), dimension(MAX_ITER_NUM) :: history
        real(8) :: x, L
        integer :: i
        logical :: stop_condition

        x = x_init
        i = 1
        history(i) = x

        stop_condition = .false.
        do while (.not. stop_condition)
            L = numerical_derivative(phi, x)
            x = (phi(x) - L * x) / (1.0d0 - L)
            i = i + 1
            history(i) = x
            stop_condition = ((x - history(i-1) <= epsilon1)       &
                             .and. (abs(f(x)) <= epsilon2))    &
                             .or. (i >= max_iter_num)
        end do

        iter = i
        root%value = history(i)
        root%error = abs(history(i) - history(i-1))

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

    subroutine atiken(f, phi, x_init, epsilon1, epsilon2, file_name, root, iter, stat)
        implicit none
        real(8), external :: f, phi
        real(8), intent(in) :: x_init, epsilon1, epsilon2
        character(50), intent(in) :: file_name
        type(result), intent(out) :: root
        integer, intent(out) :: iter
        logical, intent(out) :: stat

        real(8), dimension(MAX_ITER_NUM) :: history
        real(8) :: x
        integer :: i
        logical :: stop_condition

        x = x_init
        i = 1
        history(i) = x

        stop_condition = .false.
        do while (.not. stop_condition)
            x = phi(phi(x)) - (phi(phi(x)) - phi(x))**2 / (phi(phi(x)) - 2*phi(x) + x)
            i = i + 1
            history(i) = x
            stop_condition = ((x - history(i-1) <= epsilon1)       &
                             .and. (abs(f(x)) <= epsilon2))    &
                             .or. (i >= max_iter_num)
        end do

        iter = i
        root%value = history(i)
        root%error = abs(history(i) - history(i-1))

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

    subroutine multi_seeker(root_num, roots, seeker, left, right, step_size, epsilon1, epsilon2, file_name, iters, f, phi)
    !------------------------------------------------------------------------------
    ! This subroutine search all the possible roots in [left, right].
    ! It divides the interval into several adjoint small intervals, whose length is step_size.
    ! 
    ! Argments:
    !          seeker: the roots solver to use, including bisection, Jacobi, Newton-downhill, post acceleration, Aitken
    !          f: real(8), the original function
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
        real(8), external :: f
        real(8), external, optional :: phi
        character(20), intent(in) :: seeker
        real(8), intent(in) :: left, right, step_size, epsilon1, epsilon2
        character(50), intent(in) :: file_name
        integer, intent(out) :: root_num

        type(result), dimension(MAX_ROOTS_NUM) :: roots
        integer, dimension(MAX_ROOTS_NUM) :: iters
        type(result) :: root
        real(8) :: a, b, x_init
        integer :: iter, i, idx
        logical :: can_converge, is_converge, is_found

        if (trim(seeker) == "bisection") then
            ! Bisection
            i = 1
            a = left
            b = left + step_size
            do while (b <= right)
                call bisection(f, a, b, epsilon1, epsilon2, file_name, root, iter, is_converge)
                if (is_converge) then
                    is_found = .false.
                    do idx = 1, i-1
                        if (abs(root%value - roots(idx)%value) < 2*epsilon1) then
                            is_found = .true.
                            exit
                        end if
                    end do
                    if (.not. is_found) then
                        roots(i) = root
                        iters(i) = iter
                        i = i + 1
                    else if (root%error <= roots(idx)%error) then
                        roots(idx) = root
                        iters(idx) = iter
                    end if
                end if
                a = a + step_size
                b = b + step_size
            end do
            root_num = i - 1
        end if

        if (trim(seeker) == "downhill") then
            ! Newton downhill
            i = 1
            a = left
            b = left + step_size
            do while (b <= right)
                call random_seed()
                call random_number(x_init)
                x_init = a + x_init * (a - b)
                call newton_downhill(f, x_init, epsilon1, epsilon2, file_name, root, iter, is_converge)
                if (is_converge) then
                    is_found = .false.
                    do idx = 1, i-1
                        if (abs(root%value - roots(idx)%value) < 2*epsilon1) then
                            is_found = .true.
                            exit
                        end if
                    end do
                    if (.not. is_found) then
                        roots(i) = root
                        iters(i) = iter
                        i = i + 1
                    else if (root%error <= roots(idx)%error) then
                        roots(idx) = root
                        iters(idx) = iter
                    end if
                end if
                a = a + step_size
                b = b + step_size
            end do
            root_num = i - 1
        end if

        if (trim(seeker) == "jacobi" .or. trim(seeker) == "post" .or. trim(seeker) == "atiken") then
            i = 1
            a = left
            b = left + step_size
            do while (b <= right)
                call convergence_check(phi, a, b, can_converge)
                if (can_converge) then
                    call random_seed()
                    call random_number(x_init)
                    x_init = a + x_init * (a - b)
                    if (trim(seeker) == "jacobi") then
                        call jacobi(f, phi, x_init, epsilon1, epsilon2, file_name, root, iter, is_converge)
                    else if (trim(seeker) == "post") then
                        call post_accelerating(f, phi, x_init, epsilon1, epsilon2, file_name, root, iter, is_converge)
                    else if (trim(seeker) == "atiken") then
                        call atiken(f, phi, x_init, epsilon1, epsilon2, file_name, root, iter, is_converge)
                    end if
                    if (is_converge) then
                        is_found = .false.
                        do idx = 1, i-1
                            if (abs(root%value - roots(idx)%value) < 2*epsilon1) then
                                is_found = .true.
                                exit
                            end if
                        end do
                        if (.not. is_found) then
                            roots(i) = root
                            iters(i) = iter
                            i = i + 1
                        else if (root%error <= roots(idx)%error) then
                            roots(idx) = root
                            iters(idx) = iter
                        end if
                    end if
                end if
                a = a + step_size
                b = b + step_size
            end do
            root_num = i - 1
        end if

        return
    end subroutine

end module


program main
    use utils
    use roots_seeker

    implicit none
    real(8), external :: f, phi1, phi2
    real(8), parameter :: step_size = 0.1d0, epsilon1 = 1e-5, epsilon2 = 1e-5
    type(result), dimension(MAX_ROOTS_NUM) :: roots
    integer, dimension(MAX_ROOTS_NUM) :: iters
    real(8) :: left, right
    character(50) :: file_name = "history.dat"
    character(20) :: seeker = "post"
    integer :: iter, root_num, i
    logical :: stat

    left = -10.0d0
    right = 10.0d0

    print *, "method: "//trim(seeker)
    print *, "phi: (3*x)^{1/3}"
    call multi_seeker(root_num, roots, seeker, left, right, step_size, epsilon1, epsilon2, file_name, iters, f, phi1)
    print "(a35, i3)", "The number of roots that is found: ", root_num
    do i = 1, root_num
        print "(a8, f32.16)", "root: ", roots(i)%value
        print "(a8, f32.16)", "error: ", roots(i)%error
        print "(a8, i32)", "iter: ", iters(i)
        print *, " "
    end do

    print *, "phi: x^3/3"
    call multi_seeker(root_num, roots, seeker, left, right, step_size, epsilon1, epsilon2, file_name, iters, f, phi2)
    print "(a35, i3)", "The number of roots that is found: ", root_num
    do i = 1, root_num
        print "(a8, f32.16)", "root: ", roots(i)%value
        print "(a8, f32.16)", "error: ", roots(i)%error
        print "(a8, i32)", "iter: ", iters(i)
        print *, " "
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
real(8) function f(x)
    implicit none
    real(8), intent(in) :: x

    f = sign(abs(x)**3/3, x) - x

    return
end function

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