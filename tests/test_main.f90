program test_main
    use fpml
    use poly_zeroes
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, maxit
    ! testing variables
    integer                                     :: deg, it, clock_rate, clock_start, clock_stop
    real(kind=dp), dimension(:,:), allocatable  :: time
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: berr   
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
    open(unit=1,file="data_files/main_timing.dat")
    write(1,'(A)') 'Degree, FPML, Polzeros'
    allocate(time(maxit,2))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(p(deg+1), roots(deg), berr(deg), zeros(deg), radius(deg), h(deg+1))
        do it=1,maxit
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr)
            call system_clock(count=clock_stop)
            time(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            ! Polzeros
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
            call system_clock(count=clock_stop)
            time(it,2)=dble(clock_stop-clock_start)/dble(clock_rate)
        end do
        deallocate(p, roots, berr, zeros, radius, h)
        deg=2*deg
        write(1,'(ES15.2)', advance='no') sum(time(:,1))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(time(:,2))/maxit
    end do
    deallocate(time)
    close(1)
end program test_main