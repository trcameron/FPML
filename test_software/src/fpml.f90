!************************************************************************                       
!   FPML: Fourth order Parallelizable Modification of Laguerre's method
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 25 February 2018
!************************************************************************
!   This software is being provided as is, without any express or
!   implied warranty. In no event, may the author be liable for any error
!   in the software, or any misuse of it or any damage occuring out of 
!   its use. The entire risk of using the software lies with the party
!   doing so. 
!   Any use of the software constitutes acceptance of the terms of the
!   above statement. 
!************************************************************************
!   This module contains the following paramaters, functions, and 
!   subroutines:
!       dp: integer paramater for machine-compiler specific double 
!       precision.
!************************************************************************
module fpml
    implicit none
    integer, parameter          :: dp = KIND(0.0D0)
    real(kind=dp), parameter    :: zero = 0.0D0
    real(kind=dp), parameter    :: eps = epsilon(zero)
contains
    !************************************************
    !                       estimates               *
    !************************************************
    ! Returns initial estimates for the roots of an 
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
            if(alpha(i)>zero) then
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
    ! Returns upper envelope of the convex hull of 
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