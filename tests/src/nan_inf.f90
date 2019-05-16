!********************************************************************************
!   NAN_INF: Test FPML ability to avoid overflow/underflow and produce meaningful results.
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 16 May 2019
!********************************************************************************
! The speed and accuracy of FPML is compared against Polzeros and AMVW for
! computing the roots of the polynomial sum(i+1)x^i. 
!********************************************************************************
program nan_inf
    use fpml
    implicit none
    ! testing variables
    integer                                     :: deg, j, k
    real(kind=dp)                               :: a, r
    real(kind=dp), parameter                    :: pi2 = 6.2831853071795864769_dp
    complex(kind=dp), allocatable               :: exact_roots(:)
    ! FPML variables
    integer, parameter                          :: nitmax = 30
    integer, allocatable                        :: conv(:)
    real(kind=dp), allocatable                  :: berr(:), cond(:)
    complex(kind=dp), allocatable               :: p(:), roots(:)
    
    ! Test 1: p(z) = -a^-1 + az^deg, a = 2.0_dp**(deg)
    open(unit=1,file="data_files/nan_inf1.dat")
    write(1,'(A)') 'deg, relative forward error, backward error, condition number'
    do deg=2,538
        write(1,'(I10,A)',advance='no') deg, ', '
        ! allocate
        allocate(exact_roots(deg), conv(deg), berr(deg), cond(deg), p(deg+1), roots(deg))
        ! polynomial and roots
        a = 2.0_dp**(deg)
        r = 2.0_dp**(-2)
        exact_roots = (/ (r*cmplx(cos(pi2*k/real(deg,kind=dp)),sin(pi2*k/real(deg,kind=dp)),kind=dp), k=1,deg)/)
        p = cmplx(0,0,kind=dp)
        p(1) = cmplx(-a**(-1),0,kind=dp)
        p(deg+1) = cmplx(a,0,kind=dp)
        ! fpml
        call main(p, deg, roots, berr, cond, conv, nitmax)
        ! relative forward error, backward error, and condition number
        write(1,'(ES15.2,ES15.2,ES15.2)') max_rel_err(roots, exact_roots, deg), maxval(berr), maxval(cond)
        ! deallocate
        deallocate(exact_roots, conv, berr, cond, p, roots)
    end do
    ! close file
    close(1)
    
    ! Test 2: p(z) = -a^-1 + az^10, a = 2.0_dp**(j), j=2,538
    open(unit=1,file="data_files/nan_inf2.dat")
    write(1,'(A)') 'power, relative forward error, backward error, condition number'
    deg = 10
    do j=2,538
        write(1,'(I10,A)',advance='no') j, ', '
        ! allocate
        allocate(exact_roots(deg), conv(deg), berr(deg), cond(deg), p(deg+1), roots(deg))
        ! polynomial and roots
        a = 2.0_dp**(j)
        r = 2.0_dp**(-2.0_dp*j/real(deg,kind=dp))
        exact_roots = (/ (r*cmplx(cos(pi2*k/real(deg,kind=dp)),sin(pi2*k/real(deg,kind=dp)),kind=dp), k=1,deg)/)
        p = cmplx(0,0,kind=dp)
        p(1) = cmplx(-a**(-1),0,kind=dp)
        p(deg+1) = cmplx(a,0,kind=dp)
        ! fpml
        call main(p, deg, roots, berr, cond, conv, nitmax)
        ! relative forward error, backward error, and condition number
        write(1,'(ES15.2,ES15.2,ES15.2)') max_rel_err(roots, exact_roots, deg), maxval(berr), maxval(cond)
        ! deallocate
        deallocate(exact_roots, conv, berr, cond, p, roots)
    end do
    ! close file
    close(1)
contains
    
    !************************************************
    !                       max_rel_err             *
    !************************************************
    ! Compute maximum relative error between roots 
    ! and exact roots.
    !************************************************
    function max_rel_err(roots, exact_roots, deg) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:), roots(:)
        ! local variables
        logical                         :: m(deg)
        integer                         :: i, j, k
        real(kind=dp)                   :: delta, err, res, relerr, rmax
        ! intrinsic procedures
        intrinsic                       :: abs, max, huge
        
        ! initialize
        m = .True.
        rmax = huge(1.0_dp)
        res = eps
        ! main loop
        do i=1,deg
            ! find corresponding exact root (yet to be matched)
            err = rmax
            do j=1,deg
                if(m(j)) then
                    delta = abs(exact_roots(j) - roots(i))
                    if(delta < err) then
                        err = delta
                        k = j
                    end if
                end if
            end do
            ! mark corresponding root as matched
            m(k) = .False.
            ! calculate relative error on this root and update max
            ! zero roots give unhelpful relative errors, so revert to absolute error
            if(abs(roots(i))<eps .or. abs(exact_roots(k))<eps) then
                relerr = err
            else
                relerr = err/abs(exact_roots(k))
            end if
            res = max(relerr, res)
        end do
        res = max(res,eps)
    end function max_rel_err
end program nan_inf