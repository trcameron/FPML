program test_software
    use fpml
    use poly_zeroes
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, maxit
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    real(kind=dp), dimension(:,:), allocatable  :: time
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
    
    ! testing
    call init_random_seed()
    open(unit=1,file="data_files/main_timing.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_berr, Polzeros_time, Polzeros_berr'
    allocate(time(maxit,4))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg), zeros(deg), radius(deg), h(deg+1))
        do it=1,maxit
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            time(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,2)=maxval(berr*cond)
            ! Polzeros
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
            call system_clock(count=clock_stop)
            time(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,4)=maxval(radius)
        end do
        deallocate(p, roots, berr, cond, zeros, radius, h)
        deg=2*deg
        write(1,'(ES15.2)', advance='no') sum(time(:,1))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,2))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,3))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(time(:,4))/maxit
    end do
    deallocate(time)
    close(1)
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
end program test_software