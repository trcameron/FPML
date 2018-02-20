!************************************************
!	Thomas R. Cameron, Davidson College
!	Last Modified: 19 February 2018
!************************************************
!	Module: rand_poly contains double precision
!	kind parameter and subroutines associated 
!	with creating a random polynomial.
!************************************************
module rand_poly
	implicit none
	integer, parameter	:: dp = KIND(0.0D0)
contains
	subroutine init_random_seed()
		implicit none
		! local variables
		integer								:: i, n , clock
		integer, dimension(:), allocatable	:: seed
		! intrinsic subroutines
		intrinsic 							:: random_seed, system_clock
		
		! main
		call random_seed(size = n)
		allocate(seed(n))
		
		call system_clock(count = clock)
		seed = clock + 37 * (/ (i - 1, i = 1,n) /)
		call random_seed(put = seed)
		
		deallocate(seed)
	end subroutine init_random_seed
	
	subroutine real_rand_poly(n,x)
		implicit none
		! argument variables
		integer, intent(in)			:: n
		real(kind=dp), intent(out)	:: x(:)
		! local variables
		integer						:: k
		real(kind=dp)				:: r
		! intrinsic subroutines
		intrinsic					:: random_number
		
		! main
		do k=1,n
			call random_number(r)
			x(k) = -1 + 2 * r
		end do
	end subroutine real_rand_poly
end module rand_poly