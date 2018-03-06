program conv
    use fpml
    implicit none
    ! testing variables
    integer                                     :: deg, j
    integer, parameter                          :: itmax = 10
    real(kind=dp), parameter                    :: pi = 3.141592653589793D0
    real(kind=dp), dimension(:), allocatable    :: err
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    
    ! poly1: roots of unity deg 20
    !open(unit=1,file="data_files/conv1.dat")
    !write(1,'(A)') 'iteration, error'
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(itmax))
    p(1) = -1D0
    p(2:deg) = 0D0
    p(deg+1) = 1D0
    exact_roots = (/ (cmplx(cos(2*pi*j/deg),sin(2*pi*j/deg),kind=dp), j=1,deg)/)
    err = 0D0
    allocate(roots(deg), berr(deg), cond(deg))
    call conv_main(p, deg, roots, berr, cond, err, exact_roots)
    do j=1,itmax
        write(*,*) err(j)
    end do
    deallocate(exact_roots, p, err)
    deallocate(roots, berr, cond)
    
    write(*,*) ' '
    ! Poly 2: prescribed roots of varying scale deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(itmax))
    p(1) = 1D0/1024D0
    p(2) = -1048575D0/524288D0
    p(3) = 183251413675D0/134217728D0
    p(4) = -6862582190715075D0/17179869184D0
    p(5) = 59965700687947706355D0/1099511627776D0
    p(6) = -126769425631762997934675D0/35184372088832D0
    p(7) = 65934186820253621481357075D0/562949953421312D0
    p(8) = -8505510099812717171095062675D0/4503599627370496D0
    p(9) = 273210326382611632738979052435D0/18014398509481984D0
    p(10) = -2189425218271613769209626653075D0/36028797018963968D0
    p(11) = 4380990637147598617372537398675D0/36028797018963968D0
    p(12) = -2189425218271613769209626653075D0/18014398509481984D0
    p(13) = 273210326382611632738979052435D0/4503599627370496D0
    p(14) = -8505510099812717171095062675D0/562949953421312D0
    p(15) = 65934186820253621481357075D0/35184372088832D0
    p(16) = -126769425631762997934675D0/1099511627776D0
    p(17) = 59965700687947706355D0/17179869184D0
    p(18) = -6862582190715075D0/134217728D0
    p(19) = 183251413675D0/524288D0
    p(20) = -1048575D0/1024D0
    p(21) = 1D0
    exact_roots = (/ (cmplx(2.0D0**(-10+(j-1)),0,kind=dp), j=1,deg)/)
    err = 0D0
    allocate(roots(deg), berr(deg), cond(deg))
    call conv_main(p, deg, roots, berr, cond, err, exact_roots)
    do j=1,itmax
        write(*,*) err(j)
    end do
    deallocate(exact_roots, p, err)
    deallocate(roots, berr, cond)
    
    write(*,*) ' '
    ! Poly 3: z^i for i=0,20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(itmax))
    p = (/ (cmplx(1,0,kind=dp), j=1,deg+1)/)
    exact_roots = (/ (cmplx(cos(2.0D0*j*pi/21.0D0),sin(2.0D0*j*pi/21.0D0),kind=dp), j=1,deg)/)
    err = 0D0
    allocate(roots(deg), berr(deg), cond(deg))
    call conv_main(p, deg, roots, berr, cond, err, exact_roots)
    do j=1,itmax
        write(*,*) err(j)
    end do
    deallocate(exact_roots, p, err)
    deallocate(roots, berr, cond)
    
contains
    
    !************************************************
    !                       conv_main               *
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
        ralpha = (/ (alpha(i)*(3.8*(deg+1-i)+1), i=1,deg+1)/)
        alpha = (/ (alpha(i)*(3.8*(i-1)+1), i=1,deg+1)/)
        nz = 0
        do i=1,itmax
            call sort_err(roots, exact_roots, deg, err(i))
            if(nz==deg) return
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
                        end if
                    end if
                end if
            end do
        end do
    end subroutine conv_main
    !****************************************
    !               sort_err                *
    !****************************************
    subroutine sort_err(roots, exact_roots, deg, err)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: err
        complex(kind=dp), intent(in)    :: roots(:), exact_roots(:)
        ! local variables
        integer                         :: i, j, k
        real(kind=dp)                   :: diff, x
        complex(kind=dp)                :: temp
        complex(kind=dp)                :: temp_roots(deg)
        
        ! main
        temp_roots = roots
        do i=1,deg
            diff = abs(temp_roots(i) - exact_roots(i))
            k = i
            do j=i+1, deg
                x = abs(temp_roots(j) - exact_roots(i))
                if(x<diff) then
                    diff = x
                    k = j
                end if
            end do
            temp = temp_roots(i)
            temp_roots(i) = temp_roots(k)
            temp_roots(k) = temp
        end do
        err = maxval(abs(temp_roots - exact_roots))
    end subroutine sort_err
end program conv