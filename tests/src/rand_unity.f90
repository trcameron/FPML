!********************************************************************************
!   RAND_UNITY: Compare FPML against Polzers and AMVW on roots of unity
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 21 March 2018
!********************************************************************************
! The speed and accuracy of FPML is compared against Polzeros and AMVW for
! computing the roots of a polynomial whose roots are random in the unit circle. 
!********************************************************************************
program rand_unity
    use fpml
    use poly_zeroes
    use mpmodule
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, itnum
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    real(kind=dp), parameter                    :: pi = 3.141592653589793D0
    real(kind=dp), dimension(:), allocatable    :: coeffs, xr, xi
    real(kind=dp), dimension(:,:), allocatable  :: results
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! Polzeros variables
    integer                                     :: iter
    integer, parameter                          :: nitmax=30
    real(kind=dp), parameter                    :: small=tiny(1.0D0), big=huge(1.0D0)
    logical, dimension(:), allocatable          :: h
    real(kind=dp), dimension(:), allocatable    :: radius
    complex(kind=dp), dimension(:), allocatable :: zeros
    ! AMVW variables
    real(kind=dp), dimension(:), allocatable    :: residuals
    complex(kind=dp), dimension(:), allocatable :: poly, eigs
    
    call mpinit
    
    ! read in optional arguments
    call get_command_argument(1,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') startDegree
    else
        startDegree=10
    end if
    call get_command_argument(2,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') endDegree
    else
        endDegree=100
    end if
    call get_command_argument(3,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') itnum
    else
        itnum=512
    end if
    
    ! Testing: polynomial with random unitary complex roots
    call init_random_seed()
    open(unit=1,file="data_files/rand_unity.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, Polzeros_time, Polzeros_err, AMVW_time, AMVW_err'
    allocate(results(itnum,6))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(exact_roots(deg), coeffs(deg), xr(deg), xi(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        allocate(zeros(deg), radius(deg), h(deg+1))
        allocate(poly(deg+1), eigs(deg), residuals(deg))
        do it=1,itnum
            ! roots
            call rand_unitcmplx(exact_roots,deg)
            xr = dble(exact_roots)
            xi = aimag(exact_roots)
            ! polynomial
            call rootstocoeffs(deg,xr,xi,coeffs)
            p(1:deg) = (/ (cmplx(coeffs(j),0,kind=dp), j=1,deg)/)
            p(deg+1) = cmplx(1,0,kind=dp)
            poly = (/ (p(deg-j+1), j=0,deg)/)
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            results(it,1) = dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,2) = maxrel_fwderr(roots, exact_roots, deg)
            ! Polzeros
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,4) = maxrel_fwderr(zeros, exact_roots, deg)
            ! AMVW
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call z_poly_roots(deg, poly, eigs, residuals, flag)
            call system_clock(count=clock_stop)
            results(it,5)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,6) = maxrel_fwderr(eigs, exact_roots, deg)
        end do
        deallocate(exact_roots, coeffs, xr, xi, p, roots, berr, cond, zeros, radius, h, poly, eigs, residuals)
        deg=deg+2
        ! write results to file
        write(1,'(ES15.2)', advance='no') sum(results(:,1))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,2))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,3))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,4))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,5))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(results(:,6))/itnum
    end do
    deallocate(results)
    ! close file
    close(1)
contains
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
            array(k) = r*cmplx(cos(theta1*pi),sin(theta2*pi))
            array(k+1) = conjg(array(k))
        end do
    end subroutine rand_unitcmplx
end program rand_unity