!********************************************************************************
!   UNITY: Compare FPML against Polzers and AMVW on roots of unity
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 21 March 2018
!********************************************************************************
! The speed and accuracy of FPML is compared against Polzeros and AMVW for
! computing the roots of unity via the polynomial z^n-1.
!********************************************************************************
program unity
    use fpml
    use poly_zeroes
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, itnum
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    real(kind=dp), parameter                    :: pi = 3.141592653589793D0
    real(kind=dp), dimension(:), allocatable    :: err
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
    complex(kind=dp), dimension(:), allocatable :: coeffs, eigs
    
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
        read(arg, '(I10)') itnum
    else
        itnum=10
    end if
    
    ! Testing: roots of unity
    call init_random_seed()
    open(unit=1,file="data_files/unity.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, Polzeros_time, Polzeros_err, AMVW_time, AMVW_err'
    allocate(results(itnum,6))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(err(deg), exact_roots(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        allocate(zeros(deg), radius(deg), h(deg+1))
        allocate(coeffs(deg+1), eigs(deg), residuals(deg))
        ! polynomial and roots
        p(1) = -1D0
        p(2:deg) = 0D0
        p(deg+1) = 1D0
        coeffs = (/ (p(deg-j+1), j=0,deg)/)
        exact_roots = (/ (cmplx(cos(2*pi*j/deg),sin(2*pi*j/deg),kind=dp), j=1,deg)/)
        do it=1,itnum
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            results(it,1) = dble(clock_stop-clock_start)/dble(clock_rate)
            call sort(roots, exact_roots, deg)
            err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
            results(it,2) = sum(err)/deg
            ! Polzeros
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            call sort(zeros, exact_roots, deg)
            err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
            results(it,4) = sum(err)/deg
            ! AMVW
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call z_poly_roots(deg, coeffs, eigs, residuals, flag)
            call system_clock(count=clock_stop)
            results(it,5)=dble(clock_stop-clock_start)/dble(clock_rate)
            call sort(eigs, exact_roots, deg)
            err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
            results(it,6) = sum(err)/deg
        end do
        deallocate(err, exact_roots, p, roots, berr, cond, zeros, radius, h, coeffs, eigs, residuals)
        deg=2*deg
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
    !                       sort                    *
    !************************************************
    ! The roots are sorted with respect to exact_roots. 
    ! For each i, roots(i) is the root approximation 
    ! that is closest to exact_roots(i).
    !************************************************
    subroutine sort(roots, exact_roots, deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: i, j, k
        real(kind=dp)                   :: diff, x
        complex(kind=dp)                :: temp
        
        ! main
        do i=1,deg
            diff = abs(roots(i) - exact_roots(i))
            k = i
            do j=i+1, deg
                x = abs(roots(j) - exact_roots(i))
                if(x<diff) then
                    diff = x
                    k = j
                end if
            end do
            temp = roots(i)
            roots(i) = roots(k)
            roots(k) = temp
        end do
    end subroutine sort
end program unity