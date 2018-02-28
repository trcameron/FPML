program fpml_driver
    use fpml
    implicit none
    ! read in variables
    character(len=32)                           :: arg, poly_file
    integer                                     :: flag
    real(kind=dp)                               :: aux
    ! testing variables
    integer                                     :: deg, i, clock_rate, clock_start, clock_stop
    ! method variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    
    ! read in optional arguments
    call get_command_argument(1,arg,status=flag)
    if(flag==0) then
        ! load polynomial
        read(arg, '(A)') poly_file
        open(unit=1,file="data_files/"//poly_file, err=10, status='OLD')
        read(1,*) deg
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        do i=1, deg+1
            read(1,*) aux
            p(i) = aux
        end do
        ! run main in fpml
        open(unit=2,file="data_files/results.dat")
        write(2,'(A)') 'root_real, root_imag, berr, cond'
        call system_clock(count_rate=clock_rate)
        call system_clock(count=clock_start)
        call main(p, deg, roots, berr, cond)
        call system_clock(count=clock_stop)
        ! record time and other results
        write(*,'(A)', advance='no') 'Elapsed Time = '
        write(*,'(ES15.2)') dble(clock_stop-clock_start)/dble(clock_rate)
        do i=1,deg
            write(2,'(ES15.2)', advance='no') dble(roots(i))
            write(2,'(A)', advance='no') ','
            write(2,'(ES15.2)', advance='no') aimag(roots(i))
            write(2,'(A)', advance='no') ',' 
            write(2,'(ES15.2)', advance='no') berr(i)
            write(2,'(A)', advance='no') ','
            write(2,'(ES15.2)') cond(i)
        end do
        deallocate(p, roots, berr, cond)
        close(1)
        close(2)
        return
        ! statement only for poly_file open error
        10  write(*,*) 'Warning: File does not exist.'
        return
    else
        ! random polynomial
        write(*,*) 'Warning: No file was given. Now running random problem of degree 4096.'
        call init_random_seed()
        deg = 3200
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        call cmplx_rand_poly(deg+1,p)
        ! run main in fpml
        open(unit=2,file="data_files/results.dat")
        write(2,'(A)') 'root_real, root_imag, berr, cond'
        call system_clock(count_rate=clock_rate)
        call system_clock(count=clock_start)
        call main(p, deg, roots, berr, cond)
        call system_clock(count=clock_stop)
        ! record time and other results
        write(*,'(A)', advance='no') 'Elapsed Time = '
        write(*,'(ES15.2)') dble(clock_stop-clock_start)/dble(clock_rate)
        do i=1,deg
            write(2,'(ES15.2)', advance='no') dble(roots(i))
            write(2,'(A)', advance='no') ','
            write(2,'(ES15.2)', advance='no') aimag(roots(i))
            write(2,'(A)', advance='no') ',' 
            write(2,'(ES15.2)', advance='no') berr(i)
            write(2,'(A)', advance='no') ','
            write(2,'(ES15.2)') cond(i)
        end do
        deallocate(p, roots, berr, cond)
        close(2)
        return
    end if
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
end program fpml_driver