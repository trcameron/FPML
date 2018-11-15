!********************************************************************************
!   RAND_UNITY_NAG: Compare FPML against NAG on polynomials with
!   random complex roots in the unit circle. 
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 15 Novemeber 2018
!********************************************************************************
! The speed and accuracy of FPML is compared against NAG for
! computing the roots of a polynomial with random roots in the unit circle.  
!********************************************************************************
program rand_unity_nag
    use fpml
    use nag_library, only: c02aff
    use mpmodule
    implicit none
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, itnum
    ! testing variables
    integer                                     :: deg, it, j, clock_rate, clock_start, clock_stop
    real(kind=dp), dimension(:), allocatable    :: coeffs, err, xr, xi
    real(kind=dp), dimension(:,:), allocatable  :: results
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    integer, parameter                          :: nitmax=35
    logical, dimension(:), allocatable          :: conv
    real(kind=dp), dimension(:), allocatable    :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! NAG variables
    logical, parameter                          :: scal = .false.
    integer                                     :: ifail
    real(kind=dp), allocatable                  :: a(:,:), w(:), z(:,:)
    complex(kind=dp), allocatable               :: zeros(:)
    
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
    
    ! Testing: polynomial with random unitary complex roots in the unit circle
    call init_random_seed()
    open(unit=1,file="data_files/rand_unity_nag.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, NAG_time, NAG_err'
    allocate(results(itnum,4))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(exact_roots(deg), coeffs(deg), err(deg), xr(deg), xi(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg))
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg))
        allocate(zeros(deg))
        do it=1,itnum
            ! roots
            call rand_unitcmplx(exact_roots,deg)
            xr = dble(exact_roots)
            xi = aimag(exact_roots)
            ! polynomial
            call rootstocoeffs(deg,xr,xi,coeffs)
            p(1:deg) = (/ (cmplx(coeffs(j),0,kind=dp), j=1,deg)/)
            p(deg+1) = cmplx(1,0,kind=dp)
            do j=0,deg
                a(1,j) = dble(p(deg+1-j))
                a(2,j) = aimag(p(deg+1-j))
            end do
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond, conv, nitmax)
            call system_clock(count=clock_stop)
            results(it,1) = dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,2) = maxval(berr)
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
        end do
        deallocate(exact_roots, coeffs, err, xr, xi, p, roots, berr, cond, conv, a, w, z, zeros)
        ! write results to file
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,1))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,2))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,3))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(results(1:itnum,4))/itnum
        ! update deg
        deg=deg+2
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
    !                   error                       *
    !************************************************
    ! Computes the relative backward error for 
    ! each root approximation of p.
    !************************************************
    subroutine error(p, roots, err, deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: err(:)
        complex(kind=dp), intent(in)    :: p(:), roots(:)
        ! local variables
        integer                         :: j, k
        real(kind=dp)                   :: r, berr
        real(kind=dp), dimension(deg+1) :: alpha
        complex(kind=dp)                :: a, z
        
        ! main
        alpha = abs(p)
        alpha = (/ (alpha(j)*(3.8*(j-1)+1),j=1,deg+1)/)
        do j=1,deg
            z = roots(j)
            r = abs(z)
            if(r>1) then
                z = 1/z
                r = 1/r
                a = p(1)
                berr = alpha(1)
                do k=2,deg+1
                    a = z*a + p(k)
                    berr = r*berr + alpha(k)
                end do
            else
                a = p(deg+1)
                berr = alpha(deg+1)
                do k=deg,1,-1
                    a = z*a + p(k)
                    berr = r*berr + alpha(k)
                end do 
            end if
            err(j) = abs(a)/berr
        end do
    end subroutine error
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
        real(kind=dp), parameter        :: pi = 3.141592653589793d0
        
        ! main
        do k=1,size,2
            call random_number(theta1)
            call random_number(theta2)
            call random_number(r)
            theta1 = -1 + 2*theta1
            theta2 = -1 + 2*theta2
            array(k) = r*cmplx(cos(theta1*pi),sin(theta2*pi),kind=dp)
            array(k+1) = conjg(array(k))
        end do
    end subroutine rand_unitcmplx
end program rand_unity_nag
