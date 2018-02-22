program test_modules
	use modified_laguerre
	use poly_zeroes
	implicit none
	! test variables
	integer										:: clock, clock_rate, clock_start, clock_stop,&
												&it, j, deg, startDegree, endDegree, itmax, flag, nz
	logical, dimension(:), allocatable			:: h
	real(kind=dp)								:: berr
	complex(kind=dp)							:: abcorr, corr
	real(kind=dp), dimension(:), allocatable	:: alpha, apoly, apolyr, radius
	real(kind=dp), dimension(:,:), allocatable	:: time
	complex(kind=dp), dimension(:), allocatable	:: p, root, roots
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
	write(1,'(A)') 'CST_real, CST_imag, BST_real, BST_imag'
	deg=10
	allocate(alpha(deg+1), root(deg), roots(deg), radius(deg), h(deg+1))
	alpha(1)=1
	alpha(2)=3000
	alpha(3)=3000000
	alpha(4)=1000000000
	alpha(5)=0
	alpha(6)=0
	alpha(7)=0
	alpha(8)=0
	alpha(9)=0
	alpha(10)=0
	alpha(11)=1
	call estimates(alpha, deg, roots)
	call start(deg, alpha, root, radius, nz, small, large, h)
	do j=1,deg
		write(1,'(ES15.2)', advance='no') dble(roots(j))
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)', advance='no') aimag(roots(j)) 
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)', advance='no') dble(root(j))
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)') aimag(root(j)) 
	end do
	deallocate(alpha, root, roots, radius, h)
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
			call cmplx_rand_poly(deg+1,p)
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
	
	! test correction term
	open(unit=1,file="data_files/correction_timing.dat")
	write(1,'(A)') 'Degree, CCT, BCT'
	allocate(time(itmax,2))
	deg=startDegree
	do while(deg<=endDegree)
		write(1,'(I10)', advance='no') deg
		write(1,'(A)', advance='no') ','
		allocate(alpha(deg+1), apoly(deg+1), apolyr(deg+1), p(deg+1), roots(deg), radius(deg), h(deg))
		call cmplx_rand_poly(deg+1,p)
		do j=1,deg+1
			alpha(j)=abs(p(j))
		end do
		do j=1,deg+1
		    apolyr(deg-j+2) = eps*alpha(j)*(3.8*(deg-j+1) + 1)
		    apoly(j) = eps*alpha(j)*(3.8*(j-1) + 1)
		end do
		do it=1,itmax
			! reset roots
			call estimates(alpha, deg, roots)
			! cam correction
			call system_clock(count_rate=clock_rate)
			call system_clock(count=clock_start)
			call laguerre(p, alpha, deg, it, h(it), roots, berr)
			call system_clock(count=clock_stop)
			time(it,1)=(dble(clock_stop-clock_start)/dble(clock_rate))
			! bin correction
			call system_clock(count_rate=clock_rate)
			call system_clock(count=clock_start)
			call newton(deg, p, apoly, apolyr, roots(it), small, radius(it), corr, h(it))
			if (h(it)) then
				call aberth(deg, it, roots, abcorr)
				roots(it) = roots(it)-corr/(1-corr*abcorr)
			end if
			call system_clock(count=clock_stop)
			time(it,2)=(dble(clock_stop-clock_start)/dble(clock_rate))
		end do
		write(1,'(ES15.2)', advance='no') sum(time(:,1))/itmax
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)') sum(time(:,2))/itmax
		deg=2*deg
		deallocate(alpha, apoly, apolyr, p, roots, radius, h)
	end do
	deallocate(time)
	close(1)
end program test_modules