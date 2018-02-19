program test_modules
	use rand_poly
	implicit none
	! test variables
	integer										:: it, j, deg, startDegree, endDegree, itmax, flag
	real(kind=dp), dimension(:), allocatable	:: p
	character(len=32)							:: arg
	
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
	
	! rand_poly test
	open(unit=1,file="data_files/rand_poly.dat")
	deg=startDegree
	do while(deg<=endDegree)
		allocate(p(deg+1))
		do it=1,itmax
			call double_rand_poly(deg+1,p)
			do j=1,deg+1
				write(1,'(ES15.2)') p(j)
			end do
		end do
		deg=2*deg
		deallocate(p)
	end do
	close(1)
end program test_modules