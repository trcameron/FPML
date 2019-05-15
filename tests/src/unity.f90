!********************************************************************************
!   UNITY: Compare FPML against Polzers and AMVW on roots of unity.
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 11 Novemeber 2018
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
    real(kind=dp), dimension(:,:), allocatable  :: results
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    integer, parameter                          :: nitmax=30
    integer, dimension(:), allocatable          :: conv
    real(kind=dp), dimension(:), allocatable    :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! Polzeros variables
    integer                                     :: iter
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
        startDegree=80
    end if
    call get_command_argument(2,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') endDegree
    else
        endDegree=20480
    end if
    call get_command_argument(3,arg,status=flag)
    if(flag==0) then
        read(arg, '(I10)') itnum
    else
        itnum=512
    end if
    
    ! Testing: roots of unity
    open(unit=1,file="data_files/unity.dat")
    write(1,'(A)') 'Degree, FPML_time, FPML_err, Polzeros_time, Polzeros_err, AMVW_time, AMVW_err'
    allocate(results(itnum,6))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)', advance='no') deg
        write(1,'(A)', advance='no') ','
        allocate(exact_roots(deg))
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg))
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
            call main(p, deg, roots, berr, cond, conv, nitmax)
            call system_clock(count=clock_stop)
            results(it,1) = dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,2) = maxrel_fwderr(roots, exact_roots, deg)
            ! Polzeros
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,4) = maxrel_fwderr(zeros, exact_roots, deg)
            ! AMVW
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call z_poly_roots(deg, coeffs, eigs, residuals, flag)
            call system_clock(count=clock_stop)
            results(it,5)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,6) = maxrel_fwderr(eigs, exact_roots, deg)
        end do
        deallocate(exact_roots, p, roots, berr, cond, conv, zeros, radius, h, coeffs, eigs, residuals)
        ! write results to file
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,1))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,2))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,3))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,4))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(results(1:itnum,5))/itnum
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(results(1:itnum,6))/itnum
        ! update deg and itnum
        deg=2*deg
        itnum=itnum/2
    end do
    deallocate(results)
    ! close file
    close(1)
contains
    !************************************************
    !                       maxrel_fwderr           *
    !************************************************
    ! Compute the maximum relative forward error
    ! in the approximation roots to exact_roots.
    !************************************************
    function maxrel_fwderr(roots, exact_roots, deg) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: i, j, k
        real(kind=dp)                   :: hd, res
        
        ! main
        res = 0d0
        do i=1,deg
            hd = abs(roots(i)-exact_roots(1))
            j = 1
            do k = 2,deg
                if (hd>abs(roots(i)-exact_roots(k))) then
                    hd = abs(roots(i)-exact_roots(k))
                    j = k
                end if
            end do
            if (res<hd/abs(exact_roots(j))) then
                res = hd/abs(exact_roots(j))
            end if
        end do
        return
    end function maxrel_fwderr
end program unity
