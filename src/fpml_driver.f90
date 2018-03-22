!********************************************************************************
!   FPML_DRIVER: Driver Program for FPML
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 20 March 2018
!********************************************************************************
! Loads polynomial stored in data_files and computes the roots of said polynomial
! and returns the results along with backward error and condtion numbers into
! data_files/results.dat. The polynomial should be stored with degree on first
! line, then subsequent lines should contain the coefficients from leading to
! constant.
!********************************************************************************
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
        write(*,*) 'Starting: running on '//poly_file
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
        write(*,*) 'Finished: results are located in data_files/results.dat.'
        return
        ! statement only for poly_file open error
        10  write(*,*) 'Warning: File does not exist.'
        return
    else
        ! poly1.dat
        write(*,*) 'Warning: no file was given, running on poly1.dat.'
        open(unit=1,file="data_files/poly1.dat", status='OLD')
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
        write(*,*) 'Finished: results are located in data_files/results.dat.'
        return
    end if
end program fpml_driver