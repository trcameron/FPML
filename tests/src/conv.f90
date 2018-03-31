!********************************************************************************
!   CONV: Test the convergence rate of FPML
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 21 March 2018
!********************************************************************************
! Tests the convergence rate of FPML on three special polynomials. The first
! polynomial is Z^5-1, the second polynomial is the 10th degree Chebyshev, and 
! the third polynomial is 1+z+z^2+...+z^20.
!********************************************************************************
program conv
    use fpml
    implicit none
    ! testing variables
    integer                                     :: deg, j
    integer, parameter                          :: itmax = 8
    real(kind=dp), parameter                    :: pi = 3.141592653589793D0
    real(kind=dp), dimension(:,:), allocatable  :: err
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    
    ! Test: 3 special polynomials
    open(unit=1,file="data_files/conv.dat")
    write(1,'(A)') 'Iteration, Error-1, Error-2, Error-3'
    allocate(err(itmax,3))
    ! poly1: roots of unity deg 5
    deg = 5
    allocate(exact_roots(deg), p(deg+1), roots(deg), berr(deg), cond(deg))
    p(1) = -1D0
    p(2:deg) = 0D0
    p(deg+1) = 1D0
    exact_roots = (/ (cmplx(cos(2*pi*j/deg),sin(2*pi*j/deg),kind=dp), j=1,deg)/)
    err(:,1) = 0D0
    call conv_main(p, deg, roots, berr, cond, err(:,1), exact_roots)
    deallocate(exact_roots, p, roots, berr, cond)
    
    ! Poly 2: Chebyshev polynomial deg 10
    deg = 10
    allocate(exact_roots(deg), p(deg+1), roots(deg), berr(deg), cond(deg))
    p(1) = -1D0
    p(2) = 0D0
    p(3) = 50D0
    p(4) = 0D0
    p(5) = -400D0
    p(6) = 0D0
    p(7) = 1120D0
    p(8) = 0D0
    p(9) = -1280D0
    p(10) = 0D0
    p(11) = 512D0
    exact_roots = (/ (cmplx(cos((2D0*j-1D0)*pi/(2D0*deg)),0,kind=dp), j=1,deg)/)
    err(:,2) = 0D0
    call conv_main(p, deg, roots, berr, cond, err(:,2), exact_roots)
    deallocate(exact_roots, p, roots, berr, cond)
    
    ! Poly 3: 1+z+z^2+...+z^20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), roots(deg), berr(deg), cond(deg))
    p = (/ (cmplx(1,0,kind=dp), j=1,deg+1)/)
    exact_roots = (/ (cmplx(cos(2.0D0*j*pi/21.0D0),sin(2.0D0*j*pi/21.0D0),kind=dp), j=1,deg)/)
    err(:,3) = 0D0
    call conv_main(p, deg, roots, berr, cond, err(:,3), exact_roots)
    deallocate(exact_roots, p, roots, berr, cond)
    
    ! write results to file
    do j=1,itmax
        write(1,'(I10)', advance='no') j
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)',advance='no') err(j,1)
        write(1,'(A)',advance='no') ','
        write(1,'(ES15.2)',advance='no') err(j,2)
        write(1,'(A)',advance='no') ','
        write(1,'(ES15.2)') err(j,3)
    end do
    deallocate(err)
    close(1)
contains
    !************************************************
    !                       conv_main               *
    !************************************************
    ! Modified FPML for storing the convergence rate
    ! of the root approximations. Before each 
    ! iteration, the roots and check array are 
    ! sorted with respect to exact roots. Then the 
    ! max relative forward error is stored in err.
    !************************************************
    subroutine conv_main(p, deg, roots, berr, cond, err, exact_roots)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:), cond(:), err(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:), exact_roots(:)
        ! local variables
        integer                         :: i, j, nz
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, c, z
        
        ! main
        alpha = (/ (abs(p(i)), i = 1,deg+1) /)
        check = (/ (.true., i = 1,deg) /)
        call estimates(alpha, deg, roots)
        ralpha = (/ (alpha(i)*(1.0D0*(deg+1-i)+1), i=1,deg+1)/)
        alpha = (/ (alpha(i)*(1.0D0*(i-1)+1), i=1,deg+1)/)
        nz = 0
        do i=1,itmax
            call sort(roots, exact_roots, deg, check)
            err(i) = maxval(abs(roots-exact_roots)/abs(exact_roots))
            do j=1,deg
                if(check(j)) then
                    z = roots(j)
                    r = abs(z)
                    if(r > 1) then
                        call rcheck_lag(p, ralpha, deg, a, b, z, r, check(j), berr(j), cond(j))
                        if(check(j)) then
                            call rmodify_lag(p, deg, a, b, c, z, j, roots)
                            roots(j) = roots(j) - c
                        else
                            nz = nz + 1
                            if(nz==deg) return
                        end if
                    else
                        call check_lag(p, alpha, deg, a, b, z, r, check(j), berr(j), cond(j))
                        if(check(j)) then
                            call modify_lag(p, deg, a, b, c, z, j, roots)
                            roots(j) = roots(j) - c
                        else 
                            nz = nz + 1
                            if(nz==deg) return
                        end if
                    end if
                end if
            end do
        end do
    end subroutine conv_main
    !************************************************
    !                       sort                    *
    !************************************************
    ! The roots and check array are sorted with 
    ! respect to exact_roots. For each i, roots(i) 
    ! is the root approximation that is closest to
    ! exact_roots(i). Then check(i) is the logical
    ! element that keeps track of the convergence of
    ! that root.
    !************************************************
    subroutine sort(roots, exact_roots, deg, check)
        implicit none
        ! argument variables
        logical, intent(inout)          :: check(:)
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        logical                         :: t
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
            if(k>i) then
                temp = roots(i)
                roots(i) = roots(k)
                roots(k) = temp
                t = check(i)
                check(i) = check(k)
                check(k) = t
            end if
        end do
    end subroutine sort
end program conv