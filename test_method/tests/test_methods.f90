program test_methods
    use driver
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, maxit
    ! testing variables
    integer                                     :: deg, it, clock_rate, clock_start, clock_stop
    real(kind=dp), dimension(:,:), allocatable  :: time
    ! method variables
    real(kind=dp), dimension(:),    allocatable :: berr   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    
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
    open(unit=1,file="data/test_methods.dat")
    write(1,'(A)') 'Degree, lseq_time, lseq_berr, lcon_time, lcon_berr, aberth_time, aberth_berr'
    allocate(time(maxit,6))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)',advance='no') deg
        write(1,'(A)',advance='no') ','
        allocate(p(deg+1), roots(deg), berr(deg))
        do it=1,maxit
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            ! lseq
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main_lseq(p, deg, roots, berr)
            call system_clock(count=clock_stop)
            time(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,2)=maxval(berr)
            ! lcon
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main_lcon(p, deg, roots, berr)
            call system_clock(count=clock_stop)
            time(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,4)=maxval(berr)
            ! aberth
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main_aberth(p, deg, roots, berr)
            call system_clock(count=clock_stop)
            time(it,5)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,6)=maxval(berr)
        end do
        deallocate(p, roots, berr)
        deg=2*deg
        write(1,'(ES15.2)', advance='no') sum(time(:,1))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,2))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,3))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,4))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,5))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(time(:,6))/maxit
    end do
    deallocate(time)
    close(1)
end program test_methods