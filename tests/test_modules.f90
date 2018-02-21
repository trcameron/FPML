program test_modules
	use initial_estimates
	use poly_zeroes
	implicit none
	! test variables
	integer										:: clock, clock_rate, clock_start, clock_stop,&
												&it, j, deg, startDegree, endDegree, itmax, flag, nz
	logical, dimension(:), allocatable			:: h
	real(kind=dp), dimension(:), allocatable	:: alpha, p, radius
	real(kind=dp), dimension(:,:), allocatable	:: time
	complex(kind=dp), dimension(:), allocatable	:: root, roots
	character(len=32)							:: arg
	real(kind=dp), parameter					:: small=tiny(1.0D0), large=huge(1.0D0)
	
	call init_random_seed()
	
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
	
	! test start accuracy
	open(unit=1,file="data_files/start_accuracy.dat")
	write(1,'(A)') 'CST values, BST values'
	deg=10
	allocate(p(deg+1), root(deg), roots(deg), radius(deg), h(deg+1))
	p(1)=1
	p(2)=3000
	p(3)=3000000
	p(4)=1000000000
	p(5)=0
	p(6)=0
	p(7)=0
	p(8)=0
	p(9)=0
	p(10)=0
	p(11)=1
	call estimates(p, deg, roots)
	call start(deg, p, root, radius, nz, small, large, h)
	do j=1,deg
		write(1,'(ES15.2)', advance='no') abs(roots(j))
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)') abs(root(j))
	end do
	deallocate(p, root, roots, radius, h)
	! test start timing
	open(unit=1,file="data_files/start_timing.dat")
	write(1,'(A)') 'Degree, CST, BST'
	allocate(time(itmax,2))
	deg=startDegree
	do while(deg<=endDegree)
		write(1,'(I10)', advance='no') deg
		write(1,'(A)', advance='no') ','
		allocate(alpha(deg+1), p(deg+1), roots(deg), radius(deg), h(deg+1))
		do it=1,itmax
			call real_rand_poly(deg+1,p)
			do j=1,deg+1
				alpha(j)=abs(p(j))
			end do
			! cam start
			call system_clock(count_rate=clock_rate)
			call system_clock(count=clock_start)
			call estimates(alpha, deg, roots)
			call system_clock(count=clock_stop)
			time(it,1)=(dble(clock_stop-clock_start)/dble(clock_rate))
			! bin start
			call system_clock(count_rate=clock_rate)
			call system_clock(count=clock_start)
			call start(deg, alpha, roots, radius, nz, small, large, h)
			call system_clock(count=clock_stop)
			time(it,2)=(dble(clock_stop-clock_start)/dble(clock_rate))
		end do
		write(1,'(ES15.2)', advance='no') sum(time(:,1))/itmax
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)') sum(time(:,2))/itmax
		deg=2*deg
		deallocate(alpha, p, roots, radius, h)
	end do
	deallocate(time)
	close(1)
end program test_modules