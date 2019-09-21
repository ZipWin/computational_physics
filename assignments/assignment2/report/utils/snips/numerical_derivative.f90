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