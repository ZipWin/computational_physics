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
    character(20) :: seeker = "downhill"
    integer :: iter, root_num, i
    logical :: stat

    left = -10.0d0
    right = 10.0d0

    print *, "method: "//trim(seeker)
    call multi_seeker(root_num, roots, seeker, left, right, step_size, epsilon1, epsilon2, file_name, iters, f)
    print "(a35, i3)", "The number of roots that is found: ", root_num
    do i = 1, root_num
        print "(a8, f32.16)", "root: ", roots(i)%value
        print "(a8, f32.16)", "error: ", roots(i)%error
        print "(a8, i32)", "iter: ", iters(i)
        print *, " "
    end do

end program