!****************************************************
!   Thomas R. Cameron, Davidson College             
!   Last Modified: 22 February 2018
!****************************************************
!   Module: modified_laguerre contains the paramters 
!   and subroutines associated with computing the 
!   initial estimates for the roots of an univariate 
!   polynomial. 
!****************************************************
!   Contains the following paramaters, functions, and
!   subroutines.
!       eps: machine double precision.
!       
!       laguerre: updates a single root approximation
!       for a univariate polynomial using modified
!       laguerre correction.
!
!       modify: aberth like correction terms for
!       laguerre's method.
!****************************************************
module modified_laguerre
    use initial_estimates
    implicit none
    real(kind=dp), parameter    :: eps = epsilon(0.0D0)
contains
    !************************************
    !               laguerre            *
    !************************************
    ! Updates jth component of array
    ! roots using modified laguerre 
    ! correction term for polynomial of
    ! degree deg, whose coefficients are
    ! stored in p and their moduli is
    ! stored in alpha. 
    ! The backward error in the jth root
    ! approximation is computed and 
    ! stored in berr. If less than eps,
    ! then check is turned false and the
    ! update is halted.
    ! If |roots(j)|>1 then the reversal
    ! polynomial is used to obtain
    ! correction term.
    !************************************
    subroutine laguerre(p, alpha, deg, j, check, roots, berr)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:)
        real(kind=dp), intent(out)      :: berr
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, c, g, h, z, corr1, corr2
        
        ! main
        z = roots(j)
        r = abs(z)
        if(r>1) then
            z = 1/z
            r = 1/r
            a = p(1)
            berr = alpha(1)*(3.8*deg+1)
            do k=2,deg+1
                a = z*a+p(k)
                berr = r*berr+alpha(k)*(3.8*(deg+1-k)+1)
            end do
            berr = abs(a)/berr
            if(berr<eps) then
                check = .false.
                return
            end if
            b = deg*p(1)
            c = deg*(deg-1)*p(1)
            do k=2,deg-1
                b = z*b+(deg-k+1)*p(k)
                c = z*c+(deg-k+1)*(deg-k)*p(k)
            end do
            b = z*b+p(deg)
            b = b/a
            c = c/a
            g = z*(deg-z*b)
            h = z**2*(deg-2*z*b+z**2*(b**2-c))
        else
            a = p(deg+1)
            berr = alpha(deg+1)*(3.8*deg+1)
            do k=deg,1,-1
                a = z*a+p(k)
                berr = r*berr+alpha(k)*(3.8*(k-1)+1)
            end do
            berr = abs(a)/berr
            if(berr<eps) then
                check = .false.
                return
            end if
            b = deg*p(deg+1)
            c = deg*(deg-1)*p(deg+1)
            do k=deg,3,-1
                b = z*b+(k-1)*p(k)
                c = z*c+(k-1)*(k-2)*p(k)
            end do
            b = z*b+p(2)
            b = b/a
            c = c/a
            g = b
            h = b**2-c
        end if
        call modify(deg, j, roots, corr1, corr2)
        g = g-corr1
        h = h-corr2
        z = sqrt((deg-1)*(deg*h-g**2))
        corr1 = g-z
        corr2 = g+z
        if(abs(corr1)>abs(corr2)) then
            roots(j) = roots(j)-deg/corr1
        else
            roots(j) = roots(j)-deg/corr2
        end if
    end subroutine laguerre
    !************************************
    !               modify              *
    !************************************
    ! Computes aberth like correction
    ! terms corr1 and corr2 about the
    ! jth root approximation for a 
    ! polynomial of degree n.
    !************************************
    subroutine modify(n, j, roots, corr1, corr2)
        implicit none
        ! argument variables
        integer, intent(in)             :: n, j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr1, corr2
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: z, zj
        
        ! main
        corr1 = cmplx(0.0D0, 0.0D0, kind=dp)
        corr2 = cmplx(0.0D0, 0.0D0, kind=dp)
        zj = roots(j)
        do i=1,j-1
            z = 1/(zj - roots(i))
            corr1 = corr1 + z
            corr2 = corr2 + z**2
        end do
        do i=j+1,n
            z = 1/(zj - roots(i))
            corr1 = corr1 + z
            corr2 = corr2 + z**2
        end do
    end subroutine modify
    
end module modified_laguerre