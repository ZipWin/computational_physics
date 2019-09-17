module basic_io

    contains

    function input_matrix(ndim1, ndim2) result(matrix)
        ! Input matrix manually
        implicit none
        integer, intent(in) :: ndim1, ndim2
        integer, parameter :: ndim = 100
        integer :: i
        real*8, dimension(ndim, ndim) :: matrix

        do i = 1, ndim1
            print "(a17,i2,a8)", "Please input the ", i, 'th row: '
            read *, matrix(i, :ndim2)
        end do

        return
    end function  input_matrix

    subroutine write_matrix(matrix, file_name)
        ! Save matrix to file
        implicit none
        real*8, intent(in) :: matrix(:, :)
        character(10), intent(in) :: file_name
        integer :: i
        character(10) :: row_len

        open(file=file_name, unit=10, action="write")
        write(row_len, "(i10)") size(matrix, dim=2)
        do i = 1, size(matrix, dim=1)
            write(unit=10, fmt="("//row_len//"f8.3)") matrix(i, :)
        end do
        close(unit=10)

        return
    end subroutine write_matrix

    subroutine read_matrix(matrix, file_name)
        ! Read matrix from file
        implicit none
        real*8, intent(out) :: matrix(:, :)
        character(10), intent(in) :: file_name
        integer :: i
        character(10) :: row_len

        open(file=file_name, unit=10, action="read")
        write(row_len, "(i10)") size(matrix, dim=2)
        do i = 1, size(matrix, dim=1)
            read(unit=10, fmt="("//row_len//"f8.3)") matrix(i, :)
        end do
        close(unit=10)

        return
    end subroutine read_matrix

end module basic_io


program io
    use basic_io

    implicit none
    integer :: ndim1, ndim2
    integer, parameter :: ndim = 100
    real*8, dimension(ndim, ndim) :: matrix, matrix2
    integer :: i
    character(10), parameter :: file_name = "matrix.txt"
    ! Save the size of the 2nd dimension of matrix for format print
    character(10) :: row_len

    print "(a37)", "Please input the size of the matrix: "
    read *, ndim1, ndim2

    matrix = input_matrix(ndim1, ndim2)

    print "(15a)", "The input matrix is: "
    ! Translate the size of 2nd dimension to string
    write(row_len, "(i10)") size(matrix(:ndim1, :ndim2), dim=2)
    do i = 1, size(matrix(:ndim1, :ndim2), dim=1)
        print "("//row_len//"f8.3)", matrix(i, :ndim2)
    end do

    call write_matrix(matrix(:ndim1, :ndim2), file_name)
    print "(26a, 10a)", "Matrix have been saved in ", file_name

    call read_matrix(matrix2(:ndim1, :ndim2), file_name)
    print "(17a, 10a)", "Read matrix from ", file_name
    do i = 1, size(matrix2(:ndim1, :ndim2), dim=1)
        print "("//row_len//"f8.3)", matrix2(i, :ndim2)
    end do

end program io