!****************************************************
!   Thomas R. Cameron, Davidson College             
!   Last Modified: 25 February 2018
!****************************************************
!   Module: driver contains the paramters and 
!   subroutines associated with computing all root
!   approximations of a univariate polynomial.
!****************************************************
!   Contains the following paramaters, functions, and
!   subroutines.
!       itmax: maximum number of iterations allowed
!       for each method.
!
!       main_lseq: sequential laguerre's method.
!
!       main_lcon: concurrent laguerre's method.
!
!       main_aberth: aberth's method.
!*****************************************************
module driver
    use methods
    implicit none
    integer, parameter      :: itmax = 30
contains
    !************************************
    !               main_lseq           *
    !************************************
    subroutine main_lseq(p, deg, roots, berr)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        
        ! main
        do i=1,deg
            alpha(i) = abs(p(i))
            check(i) = .true.
        end do
        alpha(deg+1)=abs(p(deg+1))
        call estimates(alpha, deg, roots)
        do i=1,deg+1
            ralpha(i) = alpha(i)*(3.8*(deg+1-i)+1)
            alpha(i) = alpha(i)*(3.8*(i-1)+1)
        end do
        do j=1,deg
            do i=1,itmax
                if(check(j)) then
                    call laguerre_seq(p, alpha, ralpha, deg, j, check(j), roots, berr(j))
                    if(.not.check(j)) exit
                end if
            end do
        end do
    end subroutine main_lseq
    !************************************
    !               main_lcon           *
    !************************************
    subroutine main_lcon(p, deg, roots, berr)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        
        ! main
        do i=1,deg
            alpha(i) = abs(p(i))
            check(i) = .true.
        end do
        alpha(deg+1)=abs(p(deg+1))
        call estimates(alpha, deg, roots)
        do i=1,deg+1
            ralpha(i) = alpha(i)*(3.8*(deg+1-i)+1)
            alpha(i) = alpha(i)*(3.8*(i-1)+1)
        end do
        nz = 0
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    call laguerre_con(p, alpha, ralpha, deg, j, check(j), roots, berr(j))
                    if(.not.check(j)) then
                        nz = nz + 1
                        if(nz==deg) return
                    end if
                end if
            end do
        end do
    end subroutine main_lcon
    !************************************
    !               main_aberth         *
    !************************************
    subroutine main_aberth(p, deg, roots, berr)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        
        ! main
        do i=1,deg
            alpha(i) = abs(p(i))
            check(i) = .true.
        end do
        alpha(deg+1)=abs(p(deg+1))
        call estimates(alpha, deg, roots)
        do i=1,deg+1
            ralpha(i) = alpha(i)*(3.8*(deg+1-i)+1)
            alpha(i) = alpha(i)*(3.8*(i-1)+1)
        end do
        nz = 0
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    call aberth(p, alpha, ralpha, deg, j, check(j), roots, berr(j))
                    if(.not.check(j)) then
                        nz = nz + 1
                        if(nz==deg) return
                    end if
                end if
            end do
        end do
    end subroutine main_aberth
end module driver