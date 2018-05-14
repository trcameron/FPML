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
    use nag_library, only: c02aff
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
    
    ! Testing: roots of unity
    open(unit=1,file="../data_files/unity_nag.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, NAG_time, NAG_err'
    allocate(results(itnum,4))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(err(deg), exact_roots(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg))
        allocate(zeros(deg))
        ! polynomial and roots
        p(1) = -1D0
        p(2:deg) = 0D0
        p(deg+1) = 1D0
        do j=0,deg
            a(1,j) = dble(p(deg+1-j))
            a(2,j) = aimag(p(deg+1-j))
        end do
        exact_roots = (/ (cmplx(cos(2*pi*j/deg),sin(2*pi*j/deg),kind=dp), j=1,deg)/)
        do it=1,itnum
            write(*,*) deg
            ! FPML
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            results(it,1) = dble(clock_stop-clock_start)/dble(clock_rate)
            call sort(roots, exact_roots, deg)
            err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
            results(it,2) = sum(err)/deg
            write(*,*) 'FPML finished'
            ! NAG
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            ifail = 0
            call c02aff(a,deg,scal,z,w,ifail)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
            call sort(zeros, exact_roots, deg)
            err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
            results(it,4) = sum(err)/deg
            write(*,*) 'NAG finished'
        end do
        deallocate(err, exact_roots, p, roots, berr, cond, zeros, a, w, z)
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
contains
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
