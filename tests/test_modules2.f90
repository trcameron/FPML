program test_modules
	use modified_laguerre
	use poly_zeroes
	implicit none
	! read in variables
	character(len=32)				:: arg
	integer							:: flag, startDegree, endDegree, itmax
	! testing variables
	integer							:: deg, it, j, clock, clock_rate, clock_start, clock_stop
	complex(kind=dp)				:: a, b, c
	! FPML variables
	real(kind=dp), dimension(:),	allocatable	:: alpha	
	complex(kind=dp), dimension(:), allocatable	:: p, roots, zeros
	! Polzeros variables
	
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
	
	! accuracy testing
	open(unit=1,file="data_files/start_accuracy.dat")
	open(unit=2, file="data_files/eval_accuracy.dat")
	write(1,'(A)')	'CST_real, CST_imag, zeros_real, zeros_imag'
	deg=10
	allocate(p(deg+1), alpha(deg+1), roots(deg), zeros(deg))
	! polynomial coefficients and zeros
	p(1)=1; p(2)=3000; p(3)=3000000; p(4)=1000000000
	p(5)=0; p(6)=0; p(7)=0; p(8)=0; p(9)=0; p(10)=0; p(11)=1;
	alpha=abs(p)
	zeros(1)=-19.30654870154738D0
	zeros(2)=-0.001000000000100000D0
	zeros(3)=-12.03727486299114D0 - (0,1)*15.09480268810185D0
	zeros(4)=-12.03727486299114D0 + (0,1)*15.09480268810185D0
	zeros(5)=-0.0009999999999500000D0 - (0,1)*8.66025D-14
	zeros(6)=-0.0009999999999500000D0 + (0,1)*8.66025D-14
	zeros(7)=4.29663518608383D0 - (0,1)*18.82291107420092D0
	zeros(8)=4.29663518608383D0 + (0,1)*18.82291107420092D0
	zeros(9)=17.39541402768100D0 - (0,1)*8.37698350401508D0
	zeros(10)=17.39541402768100D0 + (0,1)*8.37698350401508D0
	! initial estimates
	call estimates(alpha, deg, roots)
	! record results
	do j=1,deg
		write(1,'(ES15.2)', advance='no') dble(roots(j))
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)', advance='no') aimag(roots(j)) 
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)', advance='no') dble(zeros(j))
		write(1,'(A)', advance='no') ','
		write(1,'(ES15.2)') aimag(zeros(j)) 
	end do
	deallocate(p, alpha, roots, zeros)
	close(1)
end program test_modules