program rand_poly
    use fpml
    use poly_zeroes
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, maxit
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    real(kind=dp)                               :: x, y
    real(kind=dp), dimension(:), allocatable    :: err
    real(kind=dp), dimension(:,:), allocatable  :: results
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
    
    ! testing random non-monic polynomials
    call irc()  ! initiate random seed
    open(unit=1,file="data_files/rand_poly.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, Polzeros_time, Polzeros_err, AMVW_time, AMVW_err'
    allocate(results(maxit,6))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(err(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        allocate(zeros(deg), radius(deg), h(deg+1))
        allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg))
        do it=1,maxit
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            x = dble(p(deg+1))
            y = aimag(p(deg+1))
            rcoeffs = (/ ((dble(p(deg-j+1))*x+aimag(p(deg-j+1))*y)/(x**2+y**2), j=1,deg)/)
            icoeffs = (/ ((aimag(p(deg-j+1))*x-dble(p(deg-j+1))*y)/(x**2+y**2), j=1,deg)/)
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            results(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,2)=maxval(berr*cond)
            ! Polzeros
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            call error(p, zeros, err, deg)
            results(it,4)=maxval(err)
            ! AMVW
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
            call system_clock(count=clock_stop)
            results(it,5)=dble(clock_stop-clock_start)/dble(clock_rate)
            call error(p, cmplx(reigs,ieigs, kind=dp), err, deg)
            results(it,6)=maxval(err)
        end do
        deallocate(err, p, roots, berr, cond, zeros, radius, h, rcoeffs, icoeffs, reigs, ieigs, its)
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
    deallocate(results)
    close(1)
contains
    !****************************************
    !               init_random_seed (irc)  *
    !****************************************
    subroutine irc()
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
    end subroutine irc
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
    !****************************************
    !               error                   *
    !****************************************
    subroutine error(p, roots, err, deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: err(:)
        complex(kind=dp), intent(in)    :: p(:), roots(:)
        ! local variables
        integer                         :: j, k
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, z
        
        ! main
        do j=1,deg
            z = roots(j)
            r = abs(z)
            if(r>1) then
                z = 1/z
                r = 1/r
                a = p(1)
                b = deg*p(1)
                do k=2,deg
                    a = z*a + p(k)
                    b = z*b + (deg-k+1)*p(k)
                end do
                a = z*a + p(deg+1)
            else
                a = p(deg+1)
                b = deg*p(deg+1)
                do k=deg,2,-1
                    a = z*a + p(k)
                    b = z*b + (k-1)*p(k)
                end do
                a = z*a + p(1)
            end if
            err(j) = abs(a)/(r*abs(b))
        end do
    end subroutine error
end program rand_poly