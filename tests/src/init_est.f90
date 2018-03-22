!********************************************************************************
!   INIT_EST: Initial Estimates of FPML
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 21 March 2018
!********************************************************************************
! Tests speed of initial estimates against the starting values of Polzeros. The
! values are the same in theory but computed differently. FPML uses Andrew's
! Monotone Chain algorithm to compute upper envelope of convex hull and 
! Polzers uses a Divide and Conquer Method. In addition, and example is done
! with a special polynomial.
!********************************************************************************
program init_est
    use fpml
    use poly_zeroes
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, itnum
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    real(kind=dp), dimension(:,:), allocatable  :: time
    ! FPML variables
    real(kind=dp), dimension(:), allocatable    :: alpha  
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! Polzeros variables
    integer                                     :: nz
    logical, dimension(:), allocatable          :: h
    real(kind=dp), parameter                    :: small = tiny(1.0D0), big = huge(1.0D0)
    real(kind=dp), dimension(:), allocatable    :: radius
    
    ! read in optional arguments
    call get_command_argument(1,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') startDegree
    else
        startDegree=1000
    end if
    call get_command_argument(2,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') endDegree
    else
        endDegree=1000000
    end if
    call get_command_argument(3,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') itnum
    else
        itnum=100
    end if
    
    ! Example: special polynomial
    deg = 10
    allocate(alpha(deg), p(deg+1), roots(deg), exact_roots(deg))
    p(1) = 1D0; p(2) = 3000D0; p(3) = 3000000D0; p(4) = 1000000000D0;
    p(5:10) = 0D0; p(11) = 1D0
    alpha = (/ (abs(p(j)), j=1,deg)/)
    exact_roots(1) = cmplx(-19.30654870154738D0,0,kind=dp)
    exact_roots(2) = cmplx(-0.001000000000100000D0,0,kind=dp)
    exact_roots(3) = cmplx(-12.03727486299114D0,-15.09480268810185D0,kind=dp)
    exact_roots(4) = cmplx(-12.03727486299114D0,+15.09480268810185D0,kind=dp)
    exact_roots(5) = cmplx(-0.0009999999999500000D0,-8.66025D-14,kind=dp)
    exact_roots(6) = cmplx(-0.0009999999999500000D0,+8.66025D-14,kind=dp)
    exact_roots(7) = cmplx(4.29663518608383D0,-18.82291107420092D0,kind=dp)
    exact_roots(8) = cmplx(4.29663518608383D0,+18.82291107420092D0,kind=dp)
    exact_roots(9) = cmplx(17.39541402768100D0,-8.37698350401508D0,kind=dp)
    exact_roots(10) = cmplx(17.39541402768100D0,+8.37698350401508D0,kind=dp)
    ! write estimates and exact roots to file
    open(unit=1,file="data_files/estimates_accuracy.dat")
    write(1,'(A)') 'est_real, est_imag, roots_real, roots_imag'
    call estimates(alpha, deg, roots)
    do j=1,deg
        write(1,'(ES15.2)', advance='no') dble(roots(j))
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') aimag(roots(j)) 
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') dble(exact_roots(j))
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') aimag(exact_roots(j)) 
    end do
    deallocate(alpha, p, roots, exact_roots)
    close(1)

    ! Test: random polynomials
    call init_random_seed()
    open(unit=1,file="data_files/estimates_time.dat")
    write(1,'(A)') 'Degree, estimates, start'
    allocate(time(itnum,2))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(p(deg+1), roots(deg))
        allocate(h(deg+1), radius(deg))
        do it=1,itnum
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            alpha = (/ (abs(p(j)), j=1,deg+1)/)
            ! FPML estimates
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call estimates(alpha, deg, roots)
            call system_clock(count=clock_stop)
            time(it,1) = dble(clock_stop-clock_start)/dble(clock_rate)
            ! Polzeros start
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call start(deg, alpha, roots, radius, nz, small, big, h)
            call system_clock(count=clock_stop)
            time(it,2) = dble(clock_stop-clock_start)/dble(clock_rate)
        end do
        deallocate(p, roots, h, radius)
        deg = 2*deg
        write(1,'(ES15.2)', advance='no') sum(time(:,1))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(time(:,2))/itnum
    end do
    deallocate(time)
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
    !                   cmplx_rand_poly             *
    !************************************************
    ! Creates array of random complex numbers of
    ! size n. Each complex number has real and
    ! imaginary parts uniformly distributed on (-1,1).
    !************************************************
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
end program init_est