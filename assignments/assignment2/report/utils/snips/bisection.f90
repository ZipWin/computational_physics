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

        stop_condition = (abs(right - left) <= epsilon1)    &
                         .and. (abs(f(middle)) <= epsilon2)          &
                         .and. (i <= max_iter_num)

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