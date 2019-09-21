program main
    implicit none
    real(8), external :: f
    real(8), parameter :: epsilon = 1e-5
    real(8) :: left, middle, right
    real(8) :: lv, mv, rv
    integer, parameter :: max_num = 100
    real(8), dimension(max_num) :: roots, errors
    integer :: idx = 1, roots_num
    
    left = -10.0d0
    right = 10.0d0

    call multi_bisection(left, right, epsilon, max_num, roots, errors, roots_num)

    print *, roots(1:roots_num)
    print *, errors(1:roots_num)

end program


real(8) function f(x)
    implicit none
    real(8), intent(in) :: x

    f = x * (x + 5) * (x - 3) * (x + 2.1)

    return
end function


subroutine multi_bisection(left_init, right_init, epsilon, max_num, roots, errors, roots_num)
    real(8), intent(in) :: left_init, right_init, epsilon
    integer, intent(in) :: max_num
    real(8), dimension(max_num), intent(out) :: roots, errors
    integer, intent(out) :: roots_num
    
    real(8) :: step_size, left, right, root, error
    integer :: i, idx
    logical :: stat

    step_size = (right_init - left_init) / max_num

    idx = 1
    do i = 1, max_num
        left = left_init + (i-1) * step_size
        right = left + step_size
        call bisection(left, right, epsilon, root, error, stat)
        if (stat) then
            roots(idx) = root
            errors(idx) = error
            idx = idx + 1
        end if
    end do

    roots_num = idx - 1
    
    return
end subroutine


subroutine bisection(left_init, right_init, epsilon, root, error, stat)
    implicit none
    real(8), intent(in) :: left_init, right_init, epsilon
    real(8), intent(out) :: root, error
    logical, intent(out) :: stat
    
    real(8), external :: f
    real(8) :: left, middle, right, lv, mv, rv

    left = left_init
    right = right_init

    do while (.true.)
        middle = (left + right) / 2
        lv = f(left)
        mv = f(middle)
        rv = f(right)

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

        error = right - left
        if (error <= epsilon) then
            root = middle
            stat = .true.
            exit
        end if
    end do

    return
end subroutine