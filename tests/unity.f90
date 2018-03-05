program unity
    use fpml
    use poly_zeroes
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, maxit
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
    integer, dimension(:), allocatable          :: its
    real(kind=dp), dimension(:), allocatable    :: reigs, ieigs, rcoeffs, icoeffs
    
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
    
    ! testing polynomials of the form z^n -1
    open(unit=1,file="data_files/unity.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, Polzeros_time, Polzeros_err, AMVW_time, AMVW_err'
    allocate(results(maxit,6))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(err(deg), exact_roots(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        allocate(zeros(deg), radius(deg), h(deg+1))
        allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg))
        ! polynomial and roots
        p(1) = -1D0
        p(2:deg) = 0D0
        p(deg+1) = 1D0
        rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
        icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
        exact_roots = (/ (cmplx(cos(2*pi*j/deg),sin(2*pi*j/deg),kind=dp), j=1,deg)/)
        do it=1,maxit
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
            call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
            call system_clock(count=clock_stop)
            results(it,5)=dble(clock_stop-clock_start)/dble(clock_rate)
            zeros = (/ (cmplx(reigs(j),ieigs(j),kind=dp), j=1,deg)/)
            call sort(zeros, exact_roots, deg)
            err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
            results(it,6) = sum(err)/deg
        end do
        deallocate(err, exact_roots, p, roots, berr, cond, zeros, radius, h, rcoeffs, icoeffs, reigs, ieigs, its)
        deg=2*deg
        write(1,'(ES15.2)', advance='no') sum(results(:,1))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,2))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,3))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,4))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,5))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(results(:,6))/maxit
    end do
contains
    !****************************************
    !               sort                    *
    !****************************************
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