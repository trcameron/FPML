!****************************************************
!	Thomas R. Cameron, Davidson College				
!	Last Modified: 21 February 2018
!****************************************************
!	Module: modified_laguerre contains the paramters 
!	and subroutines associated with computing the 
!	initial estimates for the roots of an univariate 
!	polynomial. 
!****************************************************
!	Contains the following paramaters, functions, and
!	subroutines.
!		modify: computes the modifications for both
!	 	Laguerre update terms.
!****************************************************
module modified_laguerre
	use initial_estimates
	implicit none
	
contains
	
	subroutine laguerre(p, alpha, deg, j, roots, error)
		implicit none
		! argument variables
		integer, intent(in)				:: deg, j
		real(kind=dp), intent(in)		:: alpha(:), error
		complex(kind=dp), intent(in)	:: p(:)
		complex(kind=dp), intent(inout)	:: roots(:)
		! local variables
		real(kind=dp)					:: r
		complex(kind=dp)				:: z
		
		! main
		z = roots(j)
		r = abs(z)
		if(r>1) then
			z=1/z
			r=1/z
		else
		end if
	end subroutine laguerre
	
	subroutine modify(n, j, roots, corr1, corr2)
		implicit none
		! argument variables
		integer, intent(in)				:: n, j
		complex(kind=dp), intent(in)	:: roots(:)
		complex(kind=dp), intent(out)	:: corr1, corr2
		! local variables
		integer							:: i
		complex(kind=dp)				:: z, zj
		
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