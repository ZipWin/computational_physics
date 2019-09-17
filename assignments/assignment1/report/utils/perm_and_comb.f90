program perm_and_comb
    implicit none
    real(8) :: n, m, P, C
    n = 12.0
    m = 8.0
    
    call permutation(n, m, P)
    call combination(n, m, C)

    print "(a4, f4.1, a6, f4.1)", "n = ", n, "  m = ", m
    print "(a12, f16.1)", "Permutation:", P
    print "(a12, f16.1)", "Combination:", C

end program


subroutine permutation(n, m, P)
    implicit none
    real(8), intent(in) :: n, m
    real(8), intent(out) :: P
    real(8), external :: mul2
    
    if (m < 0 .or. n < m) then
        print "(a)", "Illegal input!"
        return
    else
        P = mul2(n, n - m + 1)
    end if

    return
end subroutine


subroutine combination(n, m, C)
    implicit none
    real(8), intent(in) :: n, m
    real(8) :: P
    real(8), intent(out) :: C
    real(8), external :: mul2

    if (m < 0 .or. n < m) then
        print "(a)", "Illegal input!"
        return
    else
        if (m >= n/2) then
            C = mul2(n, m+1) / mul2(n-m, 1.0d0)
        else
            C = mul2(n, n-m+1) / mul2(m, 1.0d0)
        end if
    end if

    return
end subroutine


recursive function mul2(n, m) result(r)
    implicit none
    real(8), intent(in) :: n, m
    real(8) :: r
    
    if (n < m) then
        r = 1.0
    else
        r = n * mul2(n-1, m)
    end if
    
    return
end function