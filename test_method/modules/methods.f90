!****************************************************
!   Thomas R. Cameron, Davidson College             
!   Last Modified: 24 February 2018
!****************************************************
!   Module: methods contains the paramters and 
!   subroutines associated with updating the root
!   approximations of an univariate polynomial. 
!****************************************************
!   Contains the following paramaters, functions, and
!   subroutines.
!       eps: machine double precision.
!*****************************************************
module methods
    use initial_estimates
    implicit none
    real(kind=dp), parameter        :: eps = epsilon(0.0D0)
contains
    !************************************
    !               laguerre_seq        *
    !************************************
    subroutine laguerre_seq(p, alpha, ralpha, deg, j, check, roots, berr)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:), ralpha(:)
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
            berr = ralpha(1)
            do k=2,deg+1
                a = z*a+p(k)
                berr = r*berr+ralpha(k)
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
            berr = alpha(deg+1)
            do k=deg,1,-1
                a = z*a+p(k)
                berr = r*berr+alpha(k)
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
        call modify_lseq(j, roots, corr1, corr2)
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
    end subroutine laguerre_seq
    !************************************
    !               laguerre_con        *
    !************************************
    subroutine laguerre_con(p, alpha, ralpha, deg, j, check, roots, berr)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:), ralpha(:)
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
            berr = ralpha(1)
            do k=2,deg+1
                a = z*a+p(k)
                berr = r*berr+ralpha(k)
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
            berr = alpha(deg+1)
            do k=deg,1,-1
                a = z*a+p(k)
                berr = r*berr+alpha(k)
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
        call modify_lcon(deg, j, roots, corr1, corr2)
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
    end subroutine laguerre_con
    !************************************
    !               aberth              *
    !************************************
    subroutine aberth(p, alpha, ralpha, deg, j, check, roots, berr)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:), ralpha(:)
        real(kind=dp), intent(out)      :: berr
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, g, z, corr
        
        ! main
        z = roots(j)
        r = abs(z)
        if(r>1) then
            z = 1/z
            r = 1/r
            a = p(1)
            berr = ralpha(1)
            do k=2,deg+1
                a = z*a+p(k)
                berr = r*berr+ralpha(k)
            end do
            berr = abs(a)/berr
            if(berr<eps) then
                check = .false.
                return
            end if
            b = deg*p(1)
            do k=2,deg
                b = z*b+(deg-k+1)*p(k)
            end do
            g = b/a
            g = 1/(z*(deg-z*g))
        else
            a = p(deg+1)
            berr = alpha(deg+1)
            do k=deg,1,-1
                a = z*a+p(k)
                berr = r*berr+alpha(k)
            end do
            berr = abs(a)/berr
            if(berr<eps) then
                check = .false.
                return
            end if
            b = deg*p(deg+1)
            do k=deg,2,-1
                b = z*b+(k-1)*p(k)
            end do
            g = a/b
        end if
        call modify_aberth(deg, j, roots, corr)
        roots(j) = roots(j) - g/(1-g*corr)
    end subroutine aberth
    !************************************
    !               modify_lseq         *
    !************************************
    subroutine modify_lseq(j, roots, corr1, corr2)
        implicit none
        ! argument variables
        integer, intent(in)             :: j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr1, corr2
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: z, zj
        
        ! main
        corr1 = cmplx(zero, zero, kind=dp)
        corr2 = cmplx(zero, zero, kind=dp)
        zj = roots(j)
        do i=1,j-1
            z = 1/(zj - roots(i))
            corr1 = corr1 + z
            corr2 = corr2 + z**2
        end do
    end subroutine modify_lseq
    !************************************
    !               modify_lcon         *
    !************************************
    subroutine modify_lcon(n, j, roots, corr1, corr2)
        implicit none
        ! argument variables
        integer, intent(in)             :: n, j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr1, corr2
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: z, zj
        
        ! main
        corr1 = cmplx(zero, zero, kind=dp)
        corr2 = cmplx(zero, zero, kind=dp)
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
    end subroutine modify_lcon
    !************************************
    !               modify_aberth       *
    !************************************
    subroutine modify_aberth(n, j, roots, corr)
        implicit none
        ! argument variables
        integer, intent(in)             :: n, j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: z, zj
        
        ! main
        corr = cmplx(zero, zero, kind=dp)
        zj = roots(j)
        do i=1,j-1
            z = 1/(zj - roots(i))
            corr = corr + z
        end do
        do i=j+1,n
            z = 1/(zj - roots(i))
            corr = corr + z
        end do
    end subroutine modify_aberth
end module methods