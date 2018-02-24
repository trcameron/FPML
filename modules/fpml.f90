!****************************************************
!   Thomas R. Cameron, Davidson College             
!   Last Modified: 22 February 2018
!****************************************************
!   Module: fpml contains the paramters and 
!   subroutines associated with computing the roots
!   of an univariate polynomial.
!****************************************************
module fpml
    use modified_laguerre
    implicit none
    integer, parameter                  :: itmax = 30
contains
    
    subroutine main(p, deg, roots, berr)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha
        
        ! main
        do i=1,deg
            alpha(i) = abs(p(i))
            check(i) = .true.
        end do
        alpha(deg+1)=abs(p(deg+1))
        call estimates(alpha, deg, roots)
        do i=1,deg+1
            alpha(i) = alpha(i)*(3.8*(i-1) + 1)
        end do
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    call laguerre(p, alpha, deg, j, check(j), roots, berr(j))
                end if
            end do
        end do
    end subroutine main

end module fpml