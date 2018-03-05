program module
    use fpml
    implicit none
    ! testing variables
    integer                                     :: deg, j
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: alpha  
    complex(kind=dp), dimension(:), allocatable :: p, roots
    
    ! test polynomial
    deg = 10
    allocate(alpha(deg), p(deg+1), roots(deg), exact_roots(deg))
    p(1) = 1D0; p(2) = 3000D0; p(3) = 3000000D0; p(4) = 1000000000D0;
    p(5:10) = 0D0; p(11) = 1D0
    alpha = (/ (abs(p(j)), j=1,deg)/)
    exact_roots(1) = cmplx(-19.30654870154738D0,0,kind=dp)
    exact_roots(2) = cmplx(-0.001000000000100000D0,0,kind=dp)
    exact_roots(3) = cmplx(-12.03727486299114D0,-15.09480268810185D0,kind=dp)
    exact_roots(4) = cmplx(-12.03727486299114D0,+15.09480268810185D0,kind=dp)
    exact_roots(5) = cmplx(-0.0009999999999500000D0,-8.66025D-14,kind=dp)
    exact_roots(6) = cmplx(-0.0009999999999500000D0,+8.66025D-14,kind=dp)
    exact_roots(7) = cmplx(4.29663518608383D0,-18.82291107420092D0,kind=dp)
    exact_roots(8) = cmplx(4.29663518608383D0,+18.82291107420092D0,kind=dp)
    exact_roots(9) = cmplx(17.39541402768100D0,-8.37698350401508D0,kind=dp)
    exact_roots(10) = cmplx(17.39541402768100D0,+8.37698350401508D0,kind=dp)
    ! estimates accuracy
    open(unit=1,file="data_files/estimates_accuracy.dat")
    write(1,'(A)') 'est_real, est_imag, roots_real, roots_imag'
    call estimates(alpha, deg, roots)
    do j=1,deg
        write(1,'(ES15.2)', advance='no') dble(roots(j))
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') aimag(roots(j)) 
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') dble(exact_roots(j))
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') aimag(exact_roots(j)) 
    end do
    deallocate(alpha, p, roots, exact_roots)

end program module