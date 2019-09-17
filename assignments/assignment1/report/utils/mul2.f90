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