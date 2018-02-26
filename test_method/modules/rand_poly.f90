!****************************************************
!   Thomas R. Cameron, Davidson College
!   Last Modified: 25 February 2018
!****************************************************
!   Module: rand_poly contains double precision
!   kind parameter and subroutines associated 
!   with creating a random polynomial.
!****************************************************
!   Contains the following paramaters, functions, and
!   subroutines.
!       dp: integer paramater for machine-compiler 
!       specific double precision. 
!
!       init_random_seed: seeds random number 
!       generator using system_clock count parameter.
!
!       dble_rand_poly: creates n real random coeffs
!       uniformly distributed between -1 and 1 using
!       random_number.
!
!       cmplx_rang_poly: creates n cmplx random
!       coeffs whose real and imaginary parts are
!       uniformly distributed between -1 and 1 using
!       random_number.
!****************************************************
module rand_poly
    implicit none
    integer, parameter  :: dp = KIND(0.0D0)
contains
    !************************************
    !               init_random_seed    *
    !************************************
    subroutine init_random_seed()
        implicit none
        ! local variables
        integer                             :: i, n , clock
        integer, dimension(:), allocatable  :: seed
        ! intrinsic subroutines
        intrinsic                           :: random_seed, system_clock
        
        ! main
        call random_seed(size = n)
        allocate(seed(n))
        
        call system_clock(count = clock)
        seed = clock + 37 * (/ (i - 1, i = 1,n) /)
        call random_seed(put = seed)
        
        deallocate(seed)
    end subroutine init_random_seed
    !************************************
    !               dble_rand_poly      *
    !************************************
    subroutine dble_rand_poly(n,x)
        implicit none
        ! argument variables
        integer, intent(in)         :: n
        real(kind=dp), intent(out)  :: x(:)
        ! local variables
        integer                     :: k
        real(kind=dp)               :: r
        ! intrinsic subroutines
        intrinsic                   :: random_number
        
        ! main
        do k=1,n
            call random_number(r)
            x(k) = -1 + 2*r
        end do
    end subroutine dble_rand_poly
    !************************************
    !               cmplx_rand_poly     *
    !************************************
    subroutine cmplx_rand_poly(n,x)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        complex(kind=dp), intent(out)   :: x(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r1, r2
        
        ! main
        do k=1,n
            call random_number(r1)
            call random_number(r2)
            x(k)=(-1+2*r1) + (0,1)*(-1+2*r2)
        end do
    end subroutine cmplx_rand_poly
end module rand_poly