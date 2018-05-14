!********************************************************************************
!   RAND_POLY: Compare FPML against Polzers and AMVW on random polynomials
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 21 March 2018
!********************************************************************************
! The speed and accuracy of FPML is compared against NAG c02aff.
!********************************************************************************
program rand_poly
    use fpml
    use nag_library, only: c02aff
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, itnum
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    real(kind=dp), dimension(:), allocatable    :: err
    real(kind=dp), dimension(:,:), allocatable  :: results
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! NAG variables
    logical, parameter                          :: scal = .false.
    integer                                     :: ifail
    real(kind=dp), allocatable                  :: a(:,:), w(:), z(:,:)
    complex(kind=dp), allocatable               :: zeros(:)

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
    
    ! Test: random polynomials
    call init_random_seed()
    open(unit=1,file="../data_files/rand_poly_nag.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, NAG_time, NAG_err'
    allocate(results(itnum,4))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(err(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg))
        allocate(zeros(deg))
        do it=1,itnum
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            do j=0,deg
                a(1,j) = dble(p(deg+1-j))
                a(2,j) = aimag(p(deg+1-j))
            end do
            write(*,*) deg
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            results(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,2)=maxval(berr)
            write(*,*) 'FPML finished'
            ! NAG
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            ifail = 0
            call c02aff(a,deg,scal,z,w,ifail)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
            call error(p, zeros, err, deg)
            results(it,4)=maxval(err)
            write(*,*) 'NAG finished'
        end do
        deallocate(err, p, roots, berr, cond, zeros, a, w, z)
        deg=2*deg
        ! write results to file
        write(1,'(ES15.2)', advance='no') sum(results(:,1))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,2))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(:,3))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(results(:,4))/itnum
    end do
    deallocate(results)
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
    !************************************************
    !                   error                       *
    !************************************************
    ! Computes the relative forward error (upper bound)
    ! for each root of the polynomial p.
    !************************************************
    subroutine error(p, roots, err, deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: err(:)
        complex(kind=dp), intent(in)    :: p(:), roots(:)
        ! local variables
        integer                         :: j, k
        real(kind=dp)                   :: berr, r
        real(kind=dp), dimension(deg)   :: alpha
        complex(kind=dp)                :: a, b, z
        
        ! main
        do j=1,deg
            z = roots(j)
            r = abs(z)
            if(r>1) then
                alpha = (/ (abs(p(j))*(3.8*(deg-j+1)+1), j=deg+1,1,-1)/)
                z = 1/z
                r = 1/r
                a = p(1)
                b = 0
                berr = alpha(1)
                do k=deg,1,-1
                    b = z*b + a
                    a = z*a + p(deg-k+2)
                    berr = r*berr + alpha(k)
                end do
            else
                alpha = (/ (abs(p(j))*(3.8*(j-1)+1), j=1,deg+1)/)
                a = p(deg+1)
                b = 0
                berr = alpha(deg+1)
                do k=deg,1,-1
                    b = z*b + a
                    a = z*a + p(k)
                    berr = r*berr + alpha(k)
                end do
            end if
            err(j) = abs(a)/berr
        end do
    end subroutine error
end program rand_poly
