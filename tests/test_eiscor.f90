program test_eiscor
    implicit none
    ! parameters
    integer, parameter          :: dp = KIND(1.0D0)
    real(kind=dp), parameter    :: eps = epsilon(1.0D0)
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, maxit
    ! testing variables
    integer                                     :: deg, it, j
    ! AMVW variables
    real(kind=dp), dimension(:), allocatable    :: radius
    complex(kind=dp), dimension(:), allocatable :: coeffs, p, roots
    
    ! read in optional arguments
    call get_command_argument(1,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') startDegree
    else
        startDegree=100
    end if
    call get_command_argument(2,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') endDegree
    else
        endDegree=1600
    end if
    call get_command_argument(3,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') maxit
    else
        maxit=10
    end if
    
    ! main test
    deg=startDegree
    do while(deg<=endDegree)
        allocate(coeffs(deg), p(deg+1), radius(deg), roots(deg))
        do it=1,maxit
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            coeffs = (/ (p(deg-j+1), j=0,deg)/)
            ! eiscor z_poly_roots
            call z_poly_roots(deg, coeffs, roots, radius, flag)
        end do
        deallocate(coeffs, p, radius, roots)
        deg=2*deg
    end do
contains
    !****************************************
    !               init_random_seed        *
    !****************************************
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
    !****************************************
    !               cmplx_rand_poly         *
    !****************************************
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
    
end program test_eiscor