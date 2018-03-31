!********************************************************************************
!   FPML: Fourth order Parallelizable Modification of Laguerre's method
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 20 March 2018
!********************************************************************************
! MIT License
!
! Copyright (c) 2018 Thomas R. Cameron
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!********************************************************************************
!   This module contains the following paramaters, functions, and 
!   subroutines:
!       dp: integer paramater for machine-compiler specific double 
!       precision.
!
!       zero: floating point double precision 0.
!
!       eps: machine-compiler specific double precision unit roundoff.
!
!       main: main subroutine for computing roots of polynomial.
!
!       rcheck_lag: subroutine for checking backward error in approximation with
!       moduli greater than 1.
!
!       check_lag: subroutine for checking backward error in approximation with
!       moduli less than or equal to 1.
!
!       rmodify_lag: subroutine for modifying approximation with moduli greather
!       than 1.
!
!       modify_lag: subroutine for modifying approximation wtih moduli less than
!       or equal to 1.
!
!       estimates: subroutine that computes the initial root approximations.
!
!       conv_hull: subroutines that computes the indices of the upper envelope
!       of the convex hull of a 2D point set.
!
!       cross: returns 2D cross product of two vectors.
!********************************************************************************
module fpml
    implicit none
    integer, parameter          :: dp = KIND(1.0D0)
    real(kind=dp), parameter    :: eps = epsilon(1.0D0)
contains
    !************************************************
    !                       main                    *
    !************************************************
    ! Computes the roots of a polynomial of degree 
    ! deg whose coefficients are stored in p. The
    ! root approximations are stored in roots, the
    ! backward error in each approximation in berr,
    ! and the condition number of each root 
    ! approximation is stored in cond. 
    !************************************************
    subroutine main(p, deg, roots, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:), cond(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, c, z
        integer, parameter              :: itmax = 30
        
        ! main
        check = .true.
        do i=1,deg+1
            alpha(i) = abs(p(i))
        end do
        call estimates(alpha, deg, roots)
        do i=1,deg+1
            ralpha(i) = alpha(i)*(deg+2-i)
            alpha(i) = alpha(i)*i
        end do
        nz = 0
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    z = roots(j)
                    r = abs(z)
                    if(r > 1) then
                        call rcheck_lag(p, ralpha, deg, a, b, z, r, check(j), berr(j), cond(j))
                        if(check(j)) then
                            call rmodify_lag(p, deg, a, b, c, z, j, roots)
                            roots(j) = roots(j) - c
                        else
                            nz = nz + 1
                            if(nz==deg) return
                        end if
                    else
                        call check_lag(p, alpha, deg, a, b, z, r, check(j), berr(j), cond(j))
                        if(check(j)) then
                            call modify_lag(p, deg, a, b, c, z, j, roots)
                            roots(j) = roots(j) - c
                        else 
                            nz = nz + 1
                            if(nz==deg) return
                        end if
                    end if
                end if
            end do
        end do
    end subroutine main
    !************************************************
    !                       rcheck_lag              *
    !************************************************
    ! Computes backward error of root approximation
    ! with moduli greater than 1. 
    ! If the backward error is less than eps, then
    ! both backward error and condition number are
    ! computed. Otherwise, the correction term
    ! p'/p is computed.
    !************************************************
    subroutine rcheck_lag(p, ralpha, deg, a, b, z, r, check, berr, cond)
        ! argument variables
        integer, intent(in)             :: deg
        logical, intent(out)            :: check
        real(kind=dp), intent(in)       :: ralpha(:), r
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:), z
        complex(kind=dp), intent(out)   :: a, b
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: rr
        complex(kind=dp)                :: zz
        
        ! main
        zz = 1/z
        rr = 1/r
        a = p(1)
        b = deg*p(1)
        berr = ralpha(1)
        do k=2,deg
            a = zz*a + p(k)
            b = zz*b + (deg-k+1)*p(k)
            berr = rr*berr + ralpha(k)
        end do
        a = zz*a + p(deg+1)
        berr = rr*berr + ralpha(deg+1)
        if(abs(a)<eps*berr) then
            cond = berr/(rr*abs(b))
            berr = abs(a)/berr
            check = .false.
        else 
            b = b/a
        end if
    end subroutine rcheck_lag
    !************************************************
    !                       check_lag               *
    !************************************************
    ! Computes backward error of root approximation
    ! with moduli less than or equal to 1. 
    ! If the backward error is less than eps, then
    ! both backward error and condition number are
    ! computed. Otherwise, the correction term
    ! p'/p is computed.
    !************************************************
    subroutine check_lag(p, alpha, deg, a, b, z, r, check, berr, cond)
        ! argument variables
        integer, intent(in)             :: deg
        logical, intent(out)            :: check
        real(kind=dp), intent(in)       :: alpha(:), r
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:), z
        complex(kind=dp), intent(out)   :: a, b
        ! local variables
        integer                         :: k
        
        ! main
        a = p(deg+1)
        b = deg*p(deg+1)
        berr = alpha(deg+1)
        do k=deg,2,-1
            a = z*a + p(k)
            b = z*b + (k-1)*p(k)
            berr = r*berr + alpha(k)
        end do
        a = z*a + p(1)
        berr = r*berr + alpha(1)
        if(abs(a)<eps*berr) then
            cond = berr/(r*abs(b))
            berr = abs(a)/berr
            check = .false.
        else
            b = b/a
        end if
    end subroutine check_lag
    !************************************************
    !                       rmodify_lag             *
    !************************************************
    ! Computes modified Laguerre modification of
    ! the jth rooot approximation stored in z,
    ! which has moduli greater than 1.
    ! The coefficients of the polynomial of degree 
    ! deg are stored in p, all root approximations 
    ! are stored in roots. The values a and b come
    ! from rcheck_lag, c will store the second 
    ! derivative and be used to return the correction
    ! term. 
    !************************************************
    subroutine rmodify_lag(p, deg, a, b, c, z, j, roots)
        ! argument variables
        integer, intent(in)             :: deg, j
        complex(kind=dp), intent(in)    :: p(:), roots(:), a, b, z
        complex(kind=dp), intent(out)   :: c
        ! local variables
        integer                         :: k
        complex(kind=dp)                :: g, h, zz
        
        ! main
        zz = 1/z
        c = deg*(deg-1)*p(1)
        do k=2,deg-1
            c = zz*c + (deg-k+1)*(deg-k)*p(k)
        end do
        c = c/a
        g = zz*(deg-zz*b)
        h = zz**2*(deg-2*zz*b+zz**2*(b**2-c))
        do k=1,j-1
            c = 1/(z - roots(k))
            g = g - c
            h = h - c**2
        end do
        do k=j+1,deg
            c = 1/(z - roots(k))
            g = g - c
            h = h - c**2
        end do
        c = sqrt((deg-1)*(deg*h-g**2))
        h = g + c
        g = g - c
        if(abs(g)>abs(h)) then
            c = deg/g
        else
            c = deg/h
        end if
    end subroutine rmodify_lag
    !************************************************
    !                       modify_lag              *
    !************************************************
    ! Computes modified Laguerre modification of
    ! the jth rooot approximation stored in z,
    ! which has moduli less than or equal to 1. 
    ! The coefficients of the polynomial of degree 
    ! deg are stored in p, all root approximations 
    ! are stored in roots. The values a and b come
    ! from rcheck_lag, c will store the second 
    ! derivative and be used to return the correction
    ! term.
    !************************************************
    subroutine modify_lag(p, deg, a, b, c, z, j, roots)
        ! argument variables
        integer, intent(in)             :: deg, j
        complex(kind=dp), intent(in)    :: p(:), roots(:), a, b, z
        complex(kind=dp), intent(out)   :: c
        ! local variables
        integer                         :: k
        complex(kind=dp)                :: g, h
        
        ! main
        c = deg*(deg-1)*p(deg+1)
        do k=deg,3,-1
            c = z*c + (k-1)*(k-2)*p(k)
        end do
        c = c/a
        g = b
        h = b**2-c
        do k=1,j-1
            c = 1/(z - roots(k))
            g = g - c
            h = h - c**2
        end do
        do k=j+1,deg
            c = 1/(z - roots(k))
            g = g - c
            h = h - c**2
        end do
        c = sqrt((deg-1)*(deg*h-g**2))
        h = g + c
        g = g - c
        if(abs(g)>abs(h)) then
            c = deg/g
        else
            c = deg/h
        end if
    end subroutine modify_lag
    !************************************************
    !                       estimates               *
    !************************************************
    ! Computes initial estimates for the roots of an 
    ! univariate polynomial of degree deg, whose 
    ! coefficients moduli are stored in alpha. The 
    ! estimates are returned in the array roots.
    ! The computation is performed as follows: First
    ! the set (i,log(alpha(i))) is formed and the
    ! upper envelope of the convex hull of this set
    ! is computed, its indices are returned in the
    ! array h (in descending order). For i=c-1,1,-1
    ! there are h(i) - h(i+1) zeros placed on a 
    ! circle of radius alpha(h(i+1))/alpha(h(i))
    ! raised to the 1/(h(i)-h(i+1)) power. 
    !************************************************
    subroutine estimates(alpha, deg, roots)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(in)       :: alpha(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: c, i, j, k, nzeros
        real(kind=dp)                   :: ang, r, th
        integer, dimension(deg+1)       :: h
        real(kind=dp), dimension(deg+1) :: a
        real(kind=dp), parameter        :: pi2 = 6.283185307179586D0, sigma = 0.7D0
        ! intrinsic functions
        intrinsic                       :: log, cos, sin
        
        ! main
        do i=1,deg+1
            if(alpha(i)>eps) then
                a(i)=log(alpha(i))
            else
                a(i)=-1.0D30
            end if
        end do
        call conv_hull(deg+1, a, h, c)
        k=0
        th=pi2/deg
        do i=c-1,1,-1
            nzeros = h(i)-h(i+1)
            r = (alpha(h(i+1))/alpha(h(i)))**(1.0D0/nzeros)
            ang = pi2/nzeros
            DO j=1,nzeros
                roots(k+j) = r*(cos(ang*j+th*h(i)+sigma) + (0,1)*sin(ang*j+th*h(i)+sigma))
            ENDDO
            k = k+nzeros
        end do
    end subroutine estimates
    !************************************************
    !                       conv_hull               *
    !************************************************
    ! Computex upper envelope of the convex hull of 
    ! the points in the array a, which has size n.
    ! The number of vertices in the hull is equal to 
    ! c, and they are returned in the first c entries 
    ! of the array h. 
    ! The computation follows Andrew's monotone chain
    ! algorithm: Each consecutive three pairs are
    ! tested via cross to determine if they form
    ! a clockwise angle, if so that current point
    ! is rejected from the returned set. 
    !************************************************
    subroutine conv_hull(n, a, h, c)
        implicit none
        ! argument variables
        integer, intent(in)         :: n
        integer, intent(inout)      :: c
        integer, intent(inout)      :: h(:)
        real(kind=dp), intent(in)   :: a(:)
        ! local variables
        integer                     :: i
        
        ! main
        c=0
        do i=n,1,-1
            do while(c>=2 .and. cross(h, a, c, i)<eps)
                c = c - 1
            end do
            c = c + 1
            h(c) = i
        end do
    end subroutine conv_hull
    !************************************************
    !                       cross                   *
    !************************************************
    ! Returns 2D cross product of OA and OB vectors, 
    ! where
    ! O=(h(c-1),a(h(c-1))),
    ! A=(h(c),a(h(c))),
    ! B=(i,a(i)).
    ! If det>0, then OAB makes counter-clockwise turn.
    !************************************************
    function cross(h, a, c, i) result(det)
        implicit none
        ! argument variables
        integer, intent(in)         :: c, i
        integer, intent(in)         :: h(:)
        real(kind=dp), intent(in)   :: a(:)
        ! local variables
        real(kind=dp)               :: det
        
        ! main
        det = (a(i)-a(h(c-1)))*(h(c)-h(c-1)) - (a(h(c))-a(h(c-1)))*(i-h(c-1))
        return
    end function cross
end module fpml