!****************************************************
!   Thomas R. Cameron, Davidson College             
!   Last Modified: 4 March 2018
!****************************************************
!****************************************************
!   Contains the following paramaters, functions, and
!   subroutines.
!       dp: integer paramater for machine-compiler 
!       specific double precision.
!
!       zero: 0 to 15 digits.
!
!       eps: machine double precision.
!
!       itmax: maximum number of iterations allowed
!       for each method.
!       
!       cross: computes 2D cross product.
!
!       conv_hull: computes integer array that 
!       represents the vertices of the upper envelope 
!       of the convex hull of a set of points.
!
!       estimates:  computes initial estimates for 
!       the roots of an univariate polynomial.
!
!       laguerre_seq: a sequential modification of
!       laguerre's method that allows all of the
!       root approximations to be computed one at a 
!       time.
!
!       laguerre_con: a concurrent modification of
!       laguerre's method that allows progress to be
!       made on all root approximations on a single
!       iteration.
!
!       aberth: aberth method which is a concurrent
!       modification of newton's method. 
!
!       modify_lseq: modification term for sequential
!       modification of laguerres method.
!
!       modify_lcon: modification term for concurrent
!       modification of laguerres method.
!
!       modify_aberth: modification term for aberth's
!       method.
!
!       main_lseq: sequential laguerre's method.
!
!       main_lcon: concurrent laguerre's method.
!
!       main_aberth: aberth's method.
!****************************************************
program methods
    implicit none
    integer, parameter                          :: dp = KIND(1.0D0), itmax = 50
    real(kind=dp), parameter                    :: zero = 0.0D0, eps = epsilon(1.0D0)
    ! read in variables
    character(len=32)                           :: arg
    integer                                     :: flag, startDegree, endDegree, maxit
    ! testing variables
    integer                                     :: deg, it, clock_rate, clock_start, clock_stop
    real(kind=dp), dimension(:,:), allocatable  :: time
    ! method variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond
    complex(kind=dp), dimension(:), allocatable :: p, roots
    
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
        read(arg, '(I10)') maxit
    else
        maxit=10
    end if
    
    ! testing
    call irc()  ! initiate random seed
    open(unit=1,file="data_files/methods.dat")
    write(1,'(A)') 'Degree, lseq_time, lseq_err, lcon_time, lcon_err, aberth_time, aberth_err, lhan_time, lhan_err'
    allocate(time(maxit,8))
    deg=startDegree
    do while(deg<=endDegree)
        write(1,'(I10)',advance='no') deg
        write(1,'(A)',advance='no') ','
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg))
        do it=1,maxit
            ! polynomial
            call cmplx_rand_poly(deg+1,p)
            ! lseq
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main_lseq(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            time(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,2)=maxval(berr*cond)
            ! lcon
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main_lcon(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            time(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,4)=maxval(berr*cond)
            ! aberth
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main_aberth(p, deg, roots, berr, cond)
            call system_clock(count=clock_stop)
            time(it,5)=dble(clock_stop-clock_start)/dble(clock_rate)
            time(it,6)=maxval(berr*cond)
            ! lhan
            !call system_clock(count_rate=clock_rate)
            !call system_clock(count=clock_start)
            !call main_lhan(p, deg, roots, berr, cond)
            !call system_clock(count=clock_stop)
            !time(it,7)=dble(clock_stop-clock_start)/dble(clock_rate)
            !time(it,8)=maxval(berr*cond)
        end do
        deallocate(p, roots, berr, cond)
        deg=2*deg
        write(1,'(ES15.2)', advance='no') sum(time(:,1))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,2))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,3))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,4))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,5))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,6))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') sum(time(:,7))/maxit
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') sum(time(:,8))/maxit
    end do
    deallocate(time)
    close(1)
contains
    !****************************************
    !               init_random_seed (irc)  *
    !****************************************
    subroutine irc()
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
    end subroutine irc
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
    !************************************
    !               main_lseq           *
    !************************************
    subroutine main_lseq(p, deg, roots, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:), cond(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        
        ! main
        alpha = (/ (abs(p(i)), i = 1,deg+1) /)
        check = (/ (.true., i = 1,deg) /)
        call estimates(alpha, deg, roots)
        ralpha = (/ (alpha(i)*(3.8*(deg+1-i)+1), i=1,deg+1)/)
        alpha = (/ (alpha(i)*(3.8*(i-1)+1), i=1,deg+1)/)
        do j=1,deg
            do i=1,itmax
                if(check(j)) then
                    call laguerre_seq(p, alpha, ralpha, deg, j, check(j), roots, berr(j), cond(j))
                    if(.not.check(j)) exit
                end if
            end do
        end do
    end subroutine main_lseq
    !************************************
    !               main_lcon           *
    !************************************
    subroutine main_lcon(p, deg, roots, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:), cond(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        
        ! main
        alpha = (/ (abs(p(i)), i = 1,deg+1) /)
        check = (/ (.true., i = 1,deg) /)
        call estimates(alpha, deg, roots)
        ralpha = (/ (alpha(i)*(3.8*(deg+1-i)+1), i=1,deg+1)/)
        alpha = (/ (alpha(i)*(3.8*(i-1)+1), i=1,deg+1)/)
        nz = 0
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    call laguerre_con(p, alpha, ralpha, deg, j, check(j), roots, berr(j), cond(j))
                    if(.not.check(j)) then
                        nz = nz + 1
                        if(nz==deg) return
                    end if
                end if
            end do
        end do
    end subroutine main_lcon
    !************************************
    !               main_aberth         *
    !************************************
    subroutine main_aberth(p, deg, roots, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:), cond(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        
        ! main
        alpha = (/ (abs(p(i)), i = 1,deg+1) /)
        check = (/ (.true., i = 1,deg) /)
        call estimates(alpha, deg, roots)
        ralpha = (/ (alpha(i)*(3.8*(deg+1-i)+1), i=1,deg+1)/)
        alpha = (/ (alpha(i)*(3.8*(i-1)+1), i=1,deg+1)/)
        nz = 0
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    call aberth(p, alpha, ralpha, deg, j, check(j), roots, berr(j), cond(j))
                    if(.not.check(j)) then
                        nz = nz + 1
                        if(nz==deg) return
                    end if
                end if
            end do
        end do
    end subroutine main_aberth
    !************************************
    !               main_lhan           *
    !************************************
    subroutine main_lhan(p, deg, roots, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(out)      :: berr(:), cond(:)
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        logical, dimension(deg)         :: check
        real(kind=dp), dimension(deg+1) :: alpha, ralpha
        
        ! main
        alpha = (/ (abs(p(i)), i = 1,deg+1) /)
        check = (/ (.true., i = 1,deg) /)
        call estimates(alpha, deg, roots)
        ralpha = (/ (alpha(i)*(3.8*(deg+1-i)+1), i=1,deg+1)/)
        alpha = (/ (alpha(i)*(3.8*(i-1)+1), i=1,deg+1)/)
        nz = 0
        do i=1,itmax
            do j=1,deg
                if(check(j)) then
                    call laguerre_han(p, alpha, ralpha, deg, j, check(j), roots, berr(j), cond(j))
                    if(.not.check(j)) then
                        nz = nz + 1
                        if(nz==deg) return
                    end if
                end if
            end do
        end do
    end subroutine main_lhan
    !************************************
    !               laguerre_seq        *
    !************************************
    subroutine laguerre_seq(p, alpha, ralpha, deg, j, check, roots, berr, cond)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:), ralpha(:)
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, c, g, h, z, corr1, corr2
        
        ! main
        z = roots(j)
        r = abs(z)
        if(r>1) then
            z = 1/z
            r = 1/r
            a = p(1)
            b = deg*p(1)
            berr = ralpha(1)
            do k=2,deg
                a = z*a + p(k)
                b = z*b + (deg-k+1)*p(k)
                berr = r*berr + ralpha(k)
            end do
            a = z*a + p(deg+1)
            berr = r*berr + ralpha(deg+1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            c = deg*(deg-1)*p(1)
            do k=2,deg-1
                c = z*c+(deg-k+1)*(deg-k)*p(k)
            end do
            b = b/a
            c = c/a
            g = z*(deg-z*b)
            h = z**2*(deg-2*z*b+z**2*(b**2-c))
        else
            a = p(deg+1)
            b = deg*p(deg+1)
            berr = alpha(deg+1)
            do k=deg,2,-1
                a = z*a + p(k)
                b = z*b + (k-1)*p(k)
                berr = r*berr + alpha(k)
            end do
            a = z*a + p(1)
            berr = r*berr + alpha(1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            c = deg*(deg-1)*p(deg+1)
            do k=deg,3,-1
                c = z*c+(k-1)*(k-2)*p(k)
            end do
            b = b/a
            c = c/a
            g = b
            h = b**2-c
        end if
        call modify_lseq(j, roots, corr1, corr2)
        g = g-corr1
        h = h-corr2
        z = sqrt((deg-1)*(deg*h-g**2))
        corr1 = g-z
        corr2 = g+z
        if(abs(corr1)>abs(corr2)) then
            roots(j) = roots(j)-deg/corr1
        else
            roots(j) = roots(j)-deg/corr2
        end if
    end subroutine laguerre_seq
    !************************************
    !               laguerre_con        *
    !************************************
    subroutine laguerre_con(p, alpha, ralpha, deg, j, check, roots, berr, cond)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:), ralpha(:)
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, c, g, h, z, corr1, corr2
        
        ! main
        z = roots(j)
        r = abs(z)
        if(r>1) then
            z = 1/z
            r = 1/r
            a = p(1)
            b = deg*p(1)
            berr = ralpha(1)
            do k=2,deg
                a = z*a + p(k)
                b = z*b + (deg-k+1)*p(k)
                berr = r*berr + ralpha(k)
            end do
            a = z*a + p(deg+1)
            berr = r*berr + ralpha(deg+1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            c = deg*(deg-1)*p(1)
            do k=2,deg-1
                c = z*c+(deg-k+1)*(deg-k)*p(k)
            end do
            b = b/a
            c = c/a
            g = z*(deg-z*b)
            h = z**2*(deg-2*z*b+z**2*(b**2-c))
        else
            a = p(deg+1)
            b = deg*p(deg+1)
            berr = alpha(deg+1)
            do k=deg,2,-1
                a = z*a + p(k)
                b = z*b + (k-1)*p(k)
                berr = r*berr + alpha(k)
            end do
            a = z*a + p(1)
            berr = r*berr + alpha(1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            c = deg*(deg-1)*p(deg+1)
            do k=deg,3,-1
                c = z*c+(k-1)*(k-2)*p(k)
            end do
            b = b/a
            c = c/a
            g = b
            h = b**2-c
        end if
        call modify_lcon(deg, j, roots, corr1, corr2)
        g = g-corr1
        h = h-corr2
        z = sqrt((deg-1)*(deg*h-g**2))
        corr1 = g-z
        corr2 = g+z
        if(abs(corr1)>abs(corr2)) then
            roots(j) = roots(j)-deg/corr1
        else
            roots(j) = roots(j)-deg/corr2
        end if
    end subroutine laguerre_con
    !************************************
    !               aberth              *
    !************************************
    subroutine aberth(p, alpha, ralpha, deg, j, check, roots, berr, cond)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:), ralpha(:)
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, g, z, corr
        
        ! main
        z = roots(j)
        r = abs(z)
        if(r>1) then
            z = 1/z
            r = 1/r
            a = p(1)
            b = deg*p(1)
            berr = ralpha(1)
            do k=2,deg
                a = z*a + p(k)
                b = z*b + (deg-k+1)*p(k)
                berr = r*berr + ralpha(k)
            end do
            a = z*a + p(deg+1)
            berr = r*berr + ralpha(deg+1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            g = b/a
            g = 1/(z*(deg-z*g))
        else
            a = p(deg+1)
            b = deg*p(deg+1)
            berr = alpha(deg+1)
            do k=deg,2,-1
                a = z*a + p(k)
                b = z*b + (k-1)*p(k)
                berr = r*berr + alpha(k)
            end do
            a = z*a + p(1)
            berr = r*berr + alpha(1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            g = a/b
        end if
        call modify_aberth(deg, j, roots, corr)
        roots(j) = roots(j) - g/(1-g*corr)
    end subroutine aberth
    !************************************
    !               laguerre_han        *
    !************************************
    subroutine laguerre_han(p, alpha, ralpha, deg, j, check, roots, berr, cond)
        implicit none
        ! argument variables
        logical, intent(out)            :: check
        integer, intent(in)             :: deg, j
        real(kind=dp), intent(in)       :: alpha(:), ralpha(:)
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: p(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r
        complex(kind=dp)                :: a, b, c, g, h, z, corr
        
        ! main
        z = roots(j)
        r = abs(z)
        if(r>1) then
            z = 1/z
            r = 1/r
            a = p(1)
            b = deg*p(1)
            berr = ralpha(1)
            do k=2,deg
                a = z*a + p(k)
                b = z*b + (deg-k+1)*p(k)
                berr = r*berr + ralpha(k)
            end do
            a = z*a + p(deg+1)
            berr = r*berr + ralpha(deg+1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            c = deg*(deg-1)*p(1)
            do k=2,deg-1
                c = z*c+(deg-k+1)*(deg-k)*p(k)
            end do
            b = b/a
            c = c/a
            g = z*(deg-z*b)
            h = z**2*(deg-2*z*b+z**2*(b**2-c))
        else
            a = p(deg+1)
            b = deg*p(deg+1)
            berr = alpha(deg+1)
            do k=deg,2,-1
                a = z*a + p(k)
                b = z*b + (k-1)*p(k)
                berr = r*berr + alpha(k)
            end do
            a = z*a + p(1)
            berr = r*berr + alpha(1)
            if(abs(a)<berr*eps) then
                cond = berr/(r*abs(b))
                berr = abs(a)/berr
                check = .false.
                return
            end if
            c = deg*(deg-1)*p(deg+1)
            do k=deg,3,-1
                c = z*c+(k-1)*(k-2)*p(k)
            end do
            b = b/a
            c = c/a
            g = b
            h = b**2-c
        end if
        call modify_lhan(deg, j, roots, corr)
        z = sqrt((deg-1)*(deg*h-g**2-deg*corr))
        h = g - z
        g = g + z
        if(abs(g)>abs(h)) then
            roots(j) = roots(j) - deg/g
        else 
            roots(j) = roots(j) - deg/h
        end if
    end subroutine laguerre_han
    !************************************
    !               modify_lseq         *
    !************************************
    subroutine modify_lseq(j, roots, corr1, corr2)
        implicit none
        ! argument variables
        integer, intent(in)             :: j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr1, corr2
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: z, zj
        
        ! main
        corr1 = cmplx(zero, zero, kind=dp)
        corr2 = cmplx(zero, zero, kind=dp)
        zj = roots(j)
        do i=1,j-1
            z = 1/(zj - roots(i))
            corr1 = corr1 + z
            corr2 = corr2 + z**2
        end do
    end subroutine modify_lseq
    !************************************
    !               modify_lcon         *
    !************************************
    subroutine modify_lcon(n, j, roots, corr1, corr2)
        implicit none
        ! argument variables
        integer, intent(in)             :: n, j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr1, corr2
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: z, zj
        
        ! main
        corr1 = cmplx(zero, zero, kind=dp)
        corr2 = cmplx(zero, zero, kind=dp)
        zj = roots(j)
        do i=1,j-1
            z = 1/(zj - roots(i))
            corr1 = corr1 + z
            corr2 = corr2 + z**2
        end do
        do i=j+1,n
            z = 1/(zj - roots(i))
            corr1 = corr1 + z
            corr2 = corr2 + z**2
        end do
    end subroutine modify_lcon
    !************************************
    !               modify_aberth       *
    !************************************
    subroutine modify_aberth(n, j, roots, corr)
        implicit none
        ! argument variables
        integer, intent(in)             :: n, j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: z, zj
        
        ! main
        corr = cmplx(zero, zero, kind=dp)
        zj = roots(j)
        do i=1,j-1
            z = 1/(zj - roots(i))
            corr = corr + z
        end do
        do i=j+1,n
            z = 1/(zj - roots(i))
            corr = corr + z
        end do
    end subroutine modify_aberth
    !************************************
    !               modify_lhan         *
    !************************************
    subroutine modify_lhan(n, j, roots, corr)
        implicit none
        ! argument variables
        integer, intent(in)             :: n, j
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: corr
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: beta, z, zj
        
        ! main
        beta = cmplx(zero, zero, kind=dp)
        zj = roots(j)
        do i=1,j-1
            z = 1/(zj-roots(i))
            beta = beta + z
        end do
        do i=j+1,n
            z = 1/(zj-roots(i))
            beta = beta + z
        end do
        beta = beta/(n-1)
        corr = cmplx(zero, zero, kind=dp)
        do i=1,j-1
            z = 1/(zj-roots(i))
            corr = corr + (z-beta)**2
        end do
        do i=j+1,n
            z = 1/(zj-roots(i))
            corr = corr + (z-beta)**2
        end do
    end subroutine modify_lhan
    !************************************************
    !                       estimates               *
    !************************************************
    subroutine estimates(alpha, deg, roots)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(in)       :: alpha(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: c, i, j, k, nzeros
        real(kind=dp)                   :: ang, r, th
        integer, dimension(deg+1)       :: h
        real(kind=dp), dimension(deg+1) :: a
        real(kind=dp), parameter        :: pi2 = 6.283185307179586D0, sigma = 0.7D0
        ! intrinsic functions
        intrinsic                       :: log, cos, sin
        
        ! main
        do i=1,deg+1
            if(alpha(i)>eps) then
                a(i)=log(alpha(i))
            else
                a(i)=-1.0D30
            end if
        end do
        call conv_hull(deg+1, a, h, c)
        k=0
        th=pi2/deg
        do i=c-1,1,-1
            nzeros = h(i)-h(i+1)
            r = (alpha(h(i+1))/alpha(h(i)))**(1.0D0/nzeros)
            ang = pi2/nzeros
            DO j=1,nzeros
                roots(k+j) = r*(cos(ang*j+th*h(i)+sigma) + (0,1)*sin(ang*j+th*h(i)+sigma))
            ENDDO
            k = k+nzeros
        end do
    end subroutine estimates
    !************************************************
    !                       conv_hull               *
    !************************************************
    subroutine conv_hull(n, a, h, c)
        implicit none
        ! argument variables
        integer, intent(in)         :: n
        integer, intent(inout)      :: c
        integer, intent(inout)      :: h(:)
        real(kind=dp), intent(in)   :: a(:)
        ! local variables
        integer                     :: i
        
        ! main
        c=0
        do i=n,1,-1
            do while(c>=2 .and. cross(h, a, c, i)<eps)
                c = c - 1
            end do
            c = c + 1
            h(c) = i
        end do
    end subroutine conv_hull
    !************************************************
    !                       cross                   *
    !************************************************
    function cross(h, a, c, i) result(det)
        implicit none
        ! argument variables
        integer, intent(in)         :: c, i
        integer, intent(in)         :: h(:)
        real(kind=dp), intent(in)   :: a(:)
        ! local variables
        real(kind=dp)               :: det
        
        ! main
        det = (a(i)-a(h(c-1)))*(h(c)-h(c-1)) - (a(h(c))-a(h(c-1)))*(i-h(c-1))
        return
    end function cross
end program methods