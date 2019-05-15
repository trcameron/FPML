!********************************************************************************
!   CONVG: Test the convergence rate of FPML
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 1 November 2018
!********************************************************************************
! Tests the convergence rate of FPML on polynomials with random complex roots
! in the unit circle. 
!********************************************************************************
program convg
    use fpml
    use mpmodule
    implicit none
    ! testing variables
    integer                                     :: deg, j
    integer, parameter                          :: itmax = 10
    real(kind=dp), parameter                    :: pi = 3.141592653589793d0
    real(kind=dp), dimension(:), allocatable    :: coeffs, xr, xi
    real(kind=dp), dimension(:,:), allocatable  :: error
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    integer, dimension(:), allocatable          :: conv
    real(kind=dp), dimension(:), allocatable    :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    
    call mpinit
    
    ! Testing: convergence
    call init_random_seed()
    open(unit=1,file="data_files/convg.dat")
    write(1,'(A)') 'Iteration, Error-1, Error-2, Error-3, Error-4, Error-5, Error-6'
    allocate(error(itmax,6))
    deg = 4;
    j=1;
    do while(deg<=24)
        ! allocate
        allocate(exact_roots(deg), p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg), xr(deg), xi(deg), coeffs(deg))
        ! roots
        call rand_unitcmplx(exact_roots,deg)
        xr = dble(exact_roots)
        xi = aimag(exact_roots)
        ! polynomial
        call rootstocoeffs(deg,xr,xi,coeffs)
        p(1:deg) = (/ (cmplx(coeffs(j),0,kind=dp), j=1,deg)/)
        p(deg+1) = cmplx(1,0,kind=dp)
        ! error
        error(:,j)=0d0
        ! solve
        call conv_main(p, deg, roots, berr, cond, conv, itmax, error(:,j), exact_roots)
        ! deallocate
        deallocate(exact_roots, p, roots, berr, cond, conv, xr, xi, coeffs)
        ! update deg and j
        deg = deg+4
        j = j+1
    end do
    
    ! write results to file
    do j=1,itmax
        write(1,'(I10)', advance='no') j
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)',advance='no') error(j,1)
        write(1,'(A)',advance='no') ','
        write(1,'(ES15.2)',advance='no') error(j,2)
        write(1,'(A)',advance='no') ','
        write(1,'(ES15.2)',advance='no') error(j,3)
        write(1,'(A)',advance='no') ','
        write(1,'(ES15.2)',advance='no') error(j,4)
        write(1,'(A)',advance='no') ','
        write(1,'(ES15.2)',advance='no') error(j,5)
        write(1,'(A)',advance='no') ','
        write(1,'(ES15.2)') error(j,6)
    end do
    deallocate(error)
    ! close file
    close(1)
contains
    !************************************************
    !                       conv_main               *
    !************************************************
    ! Modified FPML for storing the convergence rate
    ! of the root approximations. Before each 
    ! iteration, the roots and check array are 
    ! sorted with respect to exact roots. Then the 
    ! max relative forward error is stored in err.
    !************************************************
    subroutine conv_main(p, deg, roots, berr, cond, conv, itmax, err, exact_roots)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg, itmax
        integer, intent(out)            :: conv(:)
        real(kind=dp), intent(out)      :: berr(:), cond(:), err(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:), exact_roots(:)
        ! local variables
        integer                         :: i, j, nz
        real(kind=dp), dimension(deg+1) :: alpha
        real(kind=dp)                   :: r
        complex(kind=dp)                :: b, c, z
        
        ! main
        conv = 0
        alpha = abs(p)
        call estimates(alpha, deg, roots)
        alpha = (/ (alpha(i)*(3.8*(i-1)+1),i=1,deg+1)/)
        nz = 0
        do i=1,itmax
            err(i) = maxrel_fwderr(roots, exact_roots, deg)
            do j=1,deg
                if(conv(j)==0) then
                    z = roots(j)
                    r = abs(z)
                    if(r > 1) then
                        call rcheck_lag(p, alpha, deg, b, c, z, r, conv(j), berr(j), cond(j))
                    else
                        call check_lag(p, alpha, deg, b, c, z, r, conv(j), berr(j), cond(j))
                    end if
                    if(conv(j)==0) then
                        call modify_lag(deg, b, c, z, j, roots)
                        roots(j) = roots(j) - c
                    else
                        nz = nz + 1
                        if(nz==deg) return
                    end if
                end if
            end do
        end do
    end subroutine conv_main
    !************************************************
    !                   init_random_seed            *
    !************************************************
    ! Initiate random seed using system_clock. This
    ! seed is then available for the random number
    ! generator in random_number for the life of
    ! the program.
    !************************************************
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
    !************************************************
    !                       maxrel_fwderr           *
    !************************************************
    ! Compute the maximum relative forward error
    ! in the approximation roots to exact_roots.
    !************************************************
    function maxrel_fwderr(roots, exact_roots, deg) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: i, j, k
        real(kind=dp)                   :: hd, res
        
        ! main
        res = 0d0
        do i=1,deg
            hd = abs(roots(i)-exact_roots(1))
            j = 1
            do k = 2,deg
                if (hd>abs(roots(i)-exact_roots(k))) then
                    hd = abs(roots(i)-exact_roots(k))
                    j = k
                end if
            end do
            if (res<hd/abs(exact_roots(j))) then
                res = hd/abs(exact_roots(j))
            end if
        end do
        return
    end function maxrel_fwderr
    !************************************************
    !                       rand_unitcmplx          *
    !************************************************
    ! creates a random array of unitary complex numbers.
    ! Complex numbers come in conjugate pairs, thus
    ! array should have size divisible by 2. 
    !************************************************
    subroutine rand_unitcmplx(array,size)
        implicit none
        ! argument variables
        integer, intent(in)             :: size
        complex(kind=dp), intent(out)   :: array(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: theta1, theta2, r
        
        ! main
        do k=1,size,2
            call random_number(theta1)
            call random_number(theta2)
            call random_number(r)
            theta1 = -1 + 2*theta1
            theta2 = -1 + 2*theta2
            array(k) = r*cmplx(cos(theta1*pi),sin(theta2*pi),kind=dp)
            array(k+1) = conjg(array(k))
        end do
    end subroutine rand_unitcmplx
end program convg