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
*****************************************************
module methods
    use initial_estimates
    implicit none
    real(kind=dp), parameter        :: eps = epsilon(0.0D0)
end module methods