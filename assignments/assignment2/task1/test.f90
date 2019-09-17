program main
    implicit none
    character(20) :: s1, s2
    integer :: i

    s1 = "hello.world!"

    i = 1
    do while (i <= 50)
        if (s1(i:i) == '.') exit
        i = i + 1
    end do

    s2 = s1(:i)

    print *, s1
    print *, s2

end program