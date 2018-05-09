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
    integer, parameter          :: dp = kind(1.0D0)
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
        complex(kind=dp)                :: b, c, z
        integer, parameter              :: itmax = 30
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! main
        check = .true.
        do i=1,deg+1
            alpha(i) = abs(p(i))
        end do
        call estimates(alpha, deg, roots)
        do i=1,deg+1
            ralpha(deg-i+2) = alpha(i)*(3.8*(deg-i+1) + 1)
            alpha(i) = alpha(i)*(3.8*(i-1) + 1)
        end do
        nz = 0
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    z = roots(j)
                    r = abs(z)
                    if(r > 1) then
                        call rcheck_lag(p, ralpha, deg, b, c, z, r, check(j), berr(j), cond(j))
                        if(check(j)) then
                            call modify_lag(deg, b, c, z, j, roots)
                            roots(j) = roots(j) - c
                        else
                            nz = nz + 1
                            if(nz==deg) return
                        end if
                    else
                        call check_lag(p, alpha, deg, b, c, z, r, check(j), berr(j), cond(j))
                        if(check(j)) then
                            call modify_lag(deg, b, c, z, j, roots)
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
    ! computed. Otherwise, the Laguerre correction terms
    ! Gj and Hj are computed and stored in variables
    ! b and c, respectively. 
    !************************************************
    subroutine rcheck_lag(p, ralpha, deg, b, c, z, r, check, berr, cond)
        ! argument variables
        integer, intent(in)             :: deg
        logical, intent(out)            :: check
        real(kind=dp), intent(in)       :: ralpha(:), r
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:), z
        complex(kind=dp), intent(out)   :: b, c
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: rr
        complex(kind=dp)                :: a, zz
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! main
        zz = 1/z
        rr = 1/r
        a = p(1)
        b = 0
        c = 0
        berr = ralpha(1)
        do k=deg,1,-1
            c = zz*c + b
            b = zz*b + a
            a = zz*a + p(deg-k+2)
            berr = rr*berr + ralpha(k)
        end do
        if(abs(a)<eps*berr) then
            cond = berr/(rr*abs(b))
            berr = abs(a)/berr
            check = .false.
        else 
            b = b/a
            c = 2*(c/a)
            c = zz**2*(deg-2*zz*b+zz**2*(b**2-c))
            b = zz*(deg-zz*b)
        end if
    end subroutine rcheck_lag
    !************************************************
    !                       check_lag               *
    !************************************************
    ! Computes backward error of root approximation
    ! with moduli less than or equal to 1. 
    ! If the backward error is less than eps, then
    ! both backward error and condition number are
    ! computed. Otherwise, the Laguerre correction terms
    ! Gj and Hj are computed and stored in variables
    ! b and c, respectively. 
    !************************************************
    subroutine check_lag(p, alpha, deg, b, c, z, r, check, berr, cond)
        ! argument variables
        integer, intent(in)             :: deg
        logical, intent(out)            :: check
        real(kind=dp), intent(in)       :: alpha(:), r
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:), z
        complex(kind=dp), intent(out)   :: b, c
        ! local variables
        integer                         :: k
        complex(kind=dp)                :: a
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! main
        a = p(deg+1)
        b = 0
        c = 0
        berr = alpha(deg+1)
        do k=deg,1,-1
            c = z*c + b
            b = z*b + a
            a = z*a + p(k)
            berr = r*berr + alpha(k)
        end do 
        if(abs(a)<eps*berr) then
            cond = berr/(r*abs(b))
            berr = abs(a)/berr
            check = .false.
        else
            b = b/a
            c = b**2 - 2*(c/a)
        end if
    end subroutine check_lag
    !************************************************
    !                       modify_lag              *
    !************************************************
    ! Computes modified Laguerre correction term of
    ! the jth rooot approximation.
    ! The coefficients of the polynomial of degree 
    ! deg are stored in p, all root approximations 
    ! are stored in roots. The values b, and c come
    ! from rcheck_lag or check_lag, c will be used 
    ! to return the correction term.
    !************************************************
    subroutine modify_lag(deg, b, c, z, j, roots)
        ! argument variables
        integer, intent(in)             :: deg, j
        complex(kind=dp), intent(in)    :: roots(:), z
        complex(kind=dp), intent(inout) :: b, c
        ! local variables
        integer                         :: k
        complex(kind=dp)                :: t
        ! intrinsic functions
        intrinsic                       :: abs, sqrt
        
        ! main
        do k=1,j-1
            t = 1/(z - roots(k))
            b = b - t
            c = c - t**2
        end do
        do k=j+1,deg
            t = 1/(z - roots(k))
            b = b - t
            c = c - t**2
        end do
        t = sqrt((deg-1)*(deg*c-b**2))
        c = b + t
        b = b - t
        if(abs(b)>abs(c)) then
            c = deg/b
        else
            c = deg/c
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
        real(kind=dp), parameter        :: pi2 = 6.2831853071795865, sigma = 0.7
        ! intrinsic functions
        intrinsic                       :: log, cos, sin, cmplx
        
        ! main
        do i=1,deg+1
            if(alpha(i)>0) then
                a(i) = log(alpha(i))
            else
                a(i) = -1D30
            end if
        end do
        call conv_hull(deg+1, a, h, c)
        k=0
        th=pi2/deg
        do i=c-1,1,-1
            nzeros = h(i)-h(i+1)
            r = (alpha(h(i+1))/alpha(h(i)))**(1D0/nzeros)
            ang = pi2/nzeros
            DO j=1,nzeros
                roots(k+j) = r*cmplx(cos(ang*j+th*h(i)+sigma),sin(ang*j+th*h(i)+sigma),kind=dp)
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