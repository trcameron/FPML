program test_modules
    use modified_laguerre
    use poly_zeroes
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, itmax
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    real(kind=dp), dimension(:,:), allocatable  :: time
    ! FPML variables
    real(kind=dp)                               :: berr
    real(kind=dp), dimension(:),    allocatable :: alpha, ralpha    
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! Polzeros variables
    integer                                     :: nz
    real(kind=dp), parameter                    :: small=tiny(1.0D0), big=huge(1.0D0)
    complex(kind=dp)                            :: abcorr, corr
    logical, dimension(:), allocatable          :: h
    real(kind=dp), dimension(:), allocatable    :: apoly, apolyr, radius
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
        read(arg, '(I10)') itmax
    else
        itmax=10
    end if
    
    ! accuracy testing
    open(unit=1,file="data_files/start_accuracy.dat")
    write(1,'(A)')  'CST_real, CST_imag, zeros_real, zeros_imag'
    deg=10
    allocate(p(deg+1), alpha(deg+1), roots(deg), zeros(deg))
    ! polynomial coefficients and zeros
    p(1)=1; p(2)=3000; p(3)=3000000; p(4)=1000000000
    p(5)=0; p(6)=0; p(7)=0; p(8)=0; p(9)=0; p(10)=0; p(11)=1;
    alpha=abs(p)
    zeros(1)=-19.30654870154738D0
    zeros(2)=-0.001000000000100000D0
    zeros(3)=-12.03727486299114D0 - (0,1)*15.09480268810185D0
    zeros(4)=-12.03727486299114D0 + (0,1)*15.09480268810185D0
    zeros(5)=-0.0009999999999500000D0 - (0,1)*8.66025D-14
    zeros(6)=-0.0009999999999500000D0 + (0,1)*8.66025D-14
    zeros(7)=4.29663518608383D0 - (0,1)*18.82291107420092D0
    zeros(8)=4.29663518608383D0 + (0,1)*18.82291107420092D0
    zeros(9)=17.39541402768100D0 - (0,1)*8.37698350401508D0
    zeros(10)=17.39541402768100D0 + (0,1)*8.37698350401508D0
    ! initial estimates
    call estimates(alpha, deg, roots)
    ! record results
    do j=1,deg
        write(1,'(ES15.2)', advance='no') dble(roots(j))
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') aimag(roots(j)) 
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') dble(zeros(j))
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') aimag(zeros(j)) 
    end do
    deallocate(p, alpha, roots, zeros)
    close(1)
    
    ! timing testing
    open(unit=1,file="data_files/start_timing.dat")
    open(unit=2,file="data_files/correction_timing.dat")
    write(1,'(A)') 'Degree, CST, BST'
    write(2,'(A)') 'Degree, CCT, BCT'
    allocate(time(itmax,4))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        write(2,'(I10)', advance='no') deg
        write(2,'(A)', advance='no') ','
        allocate(alpha(deg+1), ralpha(deg+1), roots(deg), zeros(deg), radius(deg), h(deg+1), p(deg+1), apoly(deg+1), apolyr(deg+1))
        ! iterations
        do it=1,itmax
            ! polynomial coefficients
            call cmplx_rand_poly(deg+1,p)
            do j=1,deg+1
                alpha(j)=abs(p(j))
                apoly(j)=abs(p(j))
            end do
            ! initialize variables
            do j=1,deg
                radius(j) = 0.0D0
                h(j) = .true.
            end do
            ! FPML start
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call estimates(alpha, deg, roots)
            call system_clock(count=clock_stop)
            time(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            ! Polzeros start
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call start(deg, alpha, zeros, radius, nz, small, big, h)
            call system_clock(count=clock_stop)
            time(it,2)=dble(clock_stop-clock_start)/dble(clock_rate)
            ! initialize variables
            do j=1,deg+1
                ralpha(j) = alpha(j)*(3.8*(deg+1-j)+1)
                alpha(j) = alpha(j)*(3.8*(j-1)+1)
                apolyr(deg-j+2) = eps*apoly(j)*(3.8*(deg-j+1) + 1)
                apoly(j) = eps*apoly(j)*(3.8*(j-1) + 1)
            end do
            do j=1,deg
                h(j) = .true.
                if(radius(j)==-1) h(j) = .false.
            end do
            ! FPML correction
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call laguerre(p, alpha, ralpha, deg, it, h(it), roots, berr)
            call system_clock(count=clock_stop)
            time(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            ! Polzeros correction
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call newton(deg, p, apoly, apolyr, zeros(it), small, radius(it), corr, h(it))
            if (h(it)) then
                call aberth(deg, it, zeros, abcorr)
                zeros(it) = zeros(it)-corr/(1-corr*abcorr)
            end if
            call system_clock(count=clock_stop)
            time(it,4)=dble(clock_stop-clock_start)/dble(clock_rate)
        end do
        deallocate(alpha, ralpha, roots, zeros, radius, h, p, apoly, apolyr)
        deg=2*deg
        write(1,'(ES15.2)', advance='no') sum(time(:,1))/itmax
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(time(:,2))/itmax
        write(2,'(ES15.2)', advance='no') sum(time(:,3))/itmax
        write(2,'(A)', advance='no') ','
        write(2,'(ES15.2)') sum(time(:,4))/itmax
    end do
    deallocate(time)
    close(1)
    close(2)
end program test_modules