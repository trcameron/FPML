!********************************************************************************
!   SPEC_POLY: Compare FPML against Polzeros and AMVW for special polynomials
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 21 March 2018
!********************************************************************************
! The accuracy of FPML is compared against NAG.
!********************************************************************************
program spec_poly
    use fpml
    use nag_library, only: c02aff
    implicit none
    ! testing variables
    integer                                     :: deg, j
    real(kind=dp), parameter                    :: pi = 3.141592653589793D0
    real(kind=dp), dimension(:), allocatable    :: err
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! NAG variables
    logical, parameter                          :: scal = .false.
    integer                                     :: ifail
    real(kind=dp), allocatable                  :: a(:,:), w(:), z(:,:)
    complex(kind=dp), allocatable               :: zeros(:)
    
    ! Testing: special polynomials
    open(unit=1,file="../data_files/spec_poly_nag.dat")
    write(1,'(A)') 'Poly No., FPML, NAG'
    ! Poly 1: Wilkinson deg 10
    deg = 10
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 3628800d0
    p(2) = -10628640d0
    p(3) = 12753576d0
    p(4) = -8409500d0
    p(5) = 3416930d0
    p(6) = -902055d0
    p(7) = 157773d0
    p(8) = -18150d0
    p(9) = 1320d0
    p(10) = -55d0
    p(11) = 1d0
    exact_roots = (/ (cmplx(j,0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '1, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 2: Wilkinson deg 15
    deg = 15
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1307674368000d0
    p(2) = 4339163001600d0
    p(3) = -6165817614720d0
    p(4) = 5056995703824d0
    p(5) = -2706813345600d0
    p(6) = 1009672107080d0
    p(7) = -272803210680d0
    p(8) = 54631129553d0
    p(9) = -8207628000d0
    p(10) = 928095740d0
    p(11) = -78558480d0
    p(12) = 4899622d0
    p(13) = -218400d0
    p(14) = 6580d0
    p(15) = -120d0
    p(16) = 1d0
    exact_roots = (/ (cmplx(j,0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '2, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 3: Wilkinson deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 2432902008176640000d0
    p(2) = -8752948036761600000d0
    p(3) = 13803759753640704000d0
    p(4) = -12870931245150988800d0
    p(5) = 8037811822645051776d0
    p(6) = -3599979517947607200d0
    p(7) = 1206647803780373360d0
    p(8) = -311333643161390640d0
    p(9) = 63030812099294896d0
    p(10) = -10142299865511450d0
    p(11) = 1307535010540395d0
    p(12) = -135585182899530d0
    p(13) = 11310276995381d0
    p(14) = -756111184500d0
    p(15) = 40171771630d0
    p(16) = -1672280820d0
    p(17) = 53327946d0
    p(18) = -1256850d0
    p(19) = 20615d0
    p(20) = -210d0
    p(21) = 1d0
    exact_roots = (/ (cmplx(j,0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '3, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 4: scaled and shifted Wilkinson deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -4.73793336560919375d-3
    p(2) = -4.74980788532250000d-3
    p(3) = 5.70184204081479750d-1
    p(4) = 5.72803665215850000d-1
    p(5) = -1.02724049585427219d+1
    p(6) = -1.04417101502222400d+1
    p(7) = 6.52029582074772000d+1
    p(8) = 6.79833441052960000d+1
    p(9) = -1.89914877310094000d+2
    p(10) = -2.07429286552800000d+2
    p(11) = 2.82271585674000000d+2
    p(12) = 3.34966323120000000d+2
    p(13) = -2.16691326540000000d+2
    p(14) = -3.01185872000000000d+2
    p(15) = 7.47163600000000000d+1
    p(16) = 1.50388800000000000d+2
    p(17) = -1.06590000000000000d0
    p(18) = -3.87600000000000000d+1
    p(19) = -5.70000000000000000d0
    p(20) = 4.00000000000000000d0
    p(21) = 1.00000000000000000d0
    exact_roots = (/ (cmplx(dble(-21+2*(j-1))/10d0,0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '4, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 5: reverse Wilkinson deg 10
    deg = 10
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 2.7557319223985890653d-7
    p(2) = -1.5156525573192239859d-5
    p(3) = 3.6375661375661375661d-4
    p(4) = -5.0016534391534391534d-3
    p(5) = 4.3478009259259259259d-2
    p(6) = -2.4858217592592592593d-1
    p(7) = 9.4161430776014109347d-1
    p(8) = -2.3174327601410934744d0
    p(9) = 3.5145436507936507937d0
    p(10) = -2.9289682539682539683d0
    p(11) = 1.0000000000000000000d0
    exact_roots = (/ (cmplx(1d0/dble(j),0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '5, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 6: reverse Wilkinson deg 15
    deg = 15
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -7.6471637318198164759d-13
    p(2) = 9.1765964781837797711d-11
    p(3) = -5.0318337355374392411d-9
    p(4) = 1.6701405590294479183d-7
    p(5) = -3.7468211658026472841d-6
    p(6) = 6.0074955908289241623d-5
    p(7) = -7.0973000825844741189d-4
    p(8) = 6.2765075165868816662d-3
    p(9) = -4.1777319254605134235d-2
    p(10) = 2.0861708186360964139d-1
    p(11) = -7.7211279182922701441d-1
    p(12) = 2.0699444845278178612d0
    p(13) = -3.8671674138251519204d0
    p(14) = 4.7151016840302554588d0
    p(15) = -3.3182289932289932290d0
    p(16) = 1.0000000000000000000d0
    exact_roots = (/ (cmplx(1d0/dble(j),0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '6, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 7: reverse Wilkinson deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 4.1103176233121648585d-19
    p(2) = -8.6316670089555462028d-17
    p(3) = 8.4734197804580278558d-15
    p(4) = -5.1660527048598944024d-13
    p(5) = 2.1919479625883946872d-11
    p(6) = -6.8736053255729181655d-10
    p(7) = 1.6511874089046065090d-8
    p(8) = -3.1078571268337857844d-7
    p(9) = 4.6488830858656684917d-6
    p(10) = -5.5729816673194132685d-5
    p(11) = 5.3743841969218428034d-4
    p(12) = -4.1688073878128312445d-3
    p(13) = 2.5907665778340944144d-2
    p(14) = -1.2796801602162448034d-1
    p(15) = 4.9597057330093876840d-1
    p(16) = -1.4797059256182881762d0
    p(17) = 3.3037959587484829179d0
    p(18) = -5.2903615525383294991d0
    p(19) = 5.6737836983356572771d0
    p(20) = -3.5977396571436819115d0
    p(21) = 1.0000000000000000000d0
    exact_roots = (/ (cmplx(1d0/dble(j),0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '7, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 8: prescribed roots of varying scale deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 9.7656250000000000000d-4
    p(2) = -1.9999980926513671875d0
    p(3) = 1.3653294270858168602d+3
    p(4) = -3.9945485714794340311d+5
    p(5) = 5.4538487063789676863d+7
    p(6) = -3.6030037799651778964d+9
    p(7) = 1.1712264370845137315d+11
    p(8) = -1.8886026297987783920d+12
    p(9) = 1.5166219745766464518d+13
    p(10) = -6.0768757200502614680d+13
    p(11) = 1.2159691690071246554d+14
    p(12) = -1.2153751440100522936d+14
    p(13) = 6.0664878983065858073d+13
    p(14) = -1.5108821038390227136d+13
    p(15) = 1.8739622993352219703d+12
    p(16) = -1.1529612095888569268d+11
    p(17) = 3.4904631720825393192d+9
    p(18) = -5.1130221714936755598d+7
    p(19) = 3.4952433333396911621d+5
    p(20) = -1.0239990234375000000d+3
    p(21) = 1.0000000000000000000d0
    exact_roots = (/ (cmplx(2d0**(-10+(j-1)),0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '8, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 9: prescribed roots of varying scale -3 deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1.5349613000519736851d+17
    p(2) = 5.5520411097852325851d+17
    p(3) = 7.3565012552135469363d+17
    p(4) = 2.6698063737425267052d+17
    p(5) = -4.1248658652425396234d+17
    p(6) = -6.4942238157287993484d+17
    p(7) = -4.4266798932313593679d+17
    p(8) = -1.7238614011013362591d+17
    p(9) = -3.5607181634056450042d+16
    p(10) = -7.6746591012597006809d+14
    p(11) = 1.5838185919451432452d+15
    p(12) = 3.9297291452486582463d+14
    p(13) = 3.2081670784934426381d+13
    p(14) = -1.7306860761064211470d+12
    p(15) = -3.9964584708819407707d+11
    p(16) = -2.9542831313947303943d+9
    p(17) = 1.3377184301714449443d+9
    p(18) = -3.3801062211972735822d+7
    p(19) = 2.9286638899803161621d+5
    p(20) = -9.6399902343750000000d+2
    p(21) = 1.0000000000000000000d0
    exact_roots = (/ (cmplx(2d0**(-10+(j-1))-3d0,0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '9, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 10: Chebyshev polynomial deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d0
    p(2) = 0d0
    p(3) = -200d0
    p(4) = 0d0
    p(5) = 6600d0
    p(6) = 0d0
    p(7) = -84480d0
    p(8) = 0d0
    p(9) = 549120d0
    p(10) = 0d0
    p(11) = -2050048d0
    p(12) = 0d0
    p(13) = 4659200d0
    p(14) = 0d0
    p(15) = -6553600d0
    p(16) = 0d0
    p(17) = 5570560d0
    p(18) = 0d0
    p(19) = -2621440d0
    p(20) = 0d0
    p(21) = 524288d0
    exact_roots = (/ (cmplx(cos((2*j-1)*pi/40d0),0,kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '10, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 11: z^i for i=0,20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p = (/ (cmplx(1,0,kind=dp), j=1,deg+1)/)
    exact_roots = (/ (cmplx(cos(2*j*pi/21d0),sin(2*j*pi/21d0),kind=dp), j=1,deg)/)
    write(1, '(A)', advance='no') '11, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 12: C. Traverso 24 MPSolve
    deg = 24
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -54765291428198020791747503747742749163073958404455022926495744d0
    p(2) = -4052135566767965847649766745769409681058667331648450681896960d0
    p(3) = -31969984081155943263834965670035075493639295858977076674560d0
    p(4) = 575060225471570237690073740639182419333523437771848417280d0
    p(5) = 7337981286595499156409929740830030318565357725459415040d0
    p(6) = 6611223380089859336490797585290455483968982077145088d0
    p(7) = -195514288747757987122118583800597358656801082441728d0
    p(8) = -726907419403715013562762609680450059293446635520d0
    p(9) = 197178719520196724204974332265013056299335680d0
    p(10) = 5968852409133617129605588058090797893943296d0
    p(11) = 16576506891508825500182005531742679597056d0
    p(12) = 23375026506968330494765978581548924928d0
    p(13) = 2206941937668751746514177591607296d0
    p(14) = -75617855277818001758431020580864d0
    p(15) = -204797687173976372829472423936d0
    p(16) = -143150263927579584306872320d0
    p(17) =  20214880144364480233472d0
    p(18) =  453786251090072698880d0
    p(19) =  1265052493274939392d0
    p(20) = -968887355572224d0
    p(21) =  1015406084096d0
    p(22) = -3949133824d0
    p(23) =  3284992d0
    p(24) = -1728d0
    p(25) = 1d0
    exact_roots(1) = cmplx(-3.52d2, 0d0, kind=dp)
    exact_roots(2) = cmplx(-3.52d2, 0d0, kind=dp)
    exact_roots(3) = cmplx(-2.8371450777d2, -2.9920517772d2, kind=dp)
    exact_roots(4) = cmplx(-2.8371450777d2,  2.9920517772d2, kind=dp)
    exact_roots(5) = cmplx(-2.7867414048d2,  6.1005469197d2, kind=dp)
    exact_roots(6) = cmplx(-2.7867414048d2, -6.1005469197d2, kind=dp)
    exact_roots(7) = cmplx(-2.74892372d2, 0d0, kind=dp)
    exact_roots(8) = cmplx(-2.014171531d2, 0d0, kind=dp)
    exact_roots(9) = cmplx(-1.255366582d2, 0d0, kind=dp)
    exact_roots(10) = cmplx(-9.599999999d1, 0d0, kind=dp)
    exact_roots(11) = cmplx(-8.8692435121d1,  5.5009607430d2, kind=dp)
    exact_roots(12) = cmplx(-8.869243512d1, -5.5009607430d2, kind=dp)
    exact_roots(13) = cmplx(-1.6000000000d1, 0d0, kind=dp)
    exact_roots(14) = cmplx(8.23178509855d1, 0d0, kind=dp)
    exact_roots(15) = cmplx(8.8692435121d1, -5.50096074303d2, kind=dp)
    exact_roots(16) = cmplx(8.8692435121d1,  5.5009607430d2, kind=dp)
    exact_roots(17) = cmplx(1.9293739373d2,  1.60865921259d3, kind=dp)
    exact_roots(18) = cmplx(1.929373937d2, -1.6086592125d3, kind=dp)
    exact_roots(19) = cmplx(2.0141715312d2, 0d0, kind=dp)
    exact_roots(20) = cmplx(2.7489237213d2, 0d0, kind=dp)
    exact_roots(21) = cmplx(7.52d2, 0d0, kind=dp)
    exact_roots(22) = cmplx(7.52d2, 0d0, kind=dp)
    exact_roots(23) = cmplx(9.1106065d2,  1.5722d0, kind=dp)
    exact_roots(24) = cmplx(9.1106065d2, -1.5722d0, kind=dp)
    write(1, '(A)', advance='no') '12, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 13: Mandelbrot 31 MPSolve
    deg = 31
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d0
    p(2) = 1d0
    p(3) = 2d0
    p(4) = 5d0
    p(5) = 14d0
    p(6) = 42d0
    p(7) = 100d0
    p(8) = 221d0
    p(9) = 470d0
    p(10) = 958d0
    p(11) = 1860d0
    p(12)= 3434d0
    p(13) = 6036d0
    p(14) = 10068d0
    p(15) = 15864d0
    p(16) = 23461d0
    p(17) = 32398d0
    p(18) = 41658d0
    p(19) = 49700d0
    p(20) = 54746d0
    p(21) = 55308d0
    p(22) = 50788d0
    p(23) = 41944d0
    p(24) = 30782d0
    p(25) = 19788d0
    p(26) = 10948d0
    p(27) = 5096d0
    p(28) = 1932d0
    p(29) = 568d0
    p(30) = 120d0
    p(31) = 16d0        
    p(32) = 1d0
    exact_roots(1) = cmplx(-1.996376137,0d0,kind=dp)
    exact_roots(2) = cmplx(-1.966773216,0d0,kind=dp)
    exact_roots(3) = cmplx(-1.907280091,0d0,kind=dp)
    exact_roots(4) = cmplx(-1.772892903,0d0,kind=dp)
    exact_roots(5) = cmplx(-1.754877666,0d0,kind=dp)
    exact_roots(6) = cmplx(-1.47601464272,0d0,kind=dp)
    exact_roots(7) = cmplx(-1.284084925525, 4.272688960406d-1,kind=dp)
    exact_roots(8) = cmplx(-1.284084925525,-4.272688960406d-1,kind=dp)
    exact_roots(9) = cmplx(-1.138000666650,-2.403324012620d-1,kind=dp)
    exact_roots(10) = cmplx(-1.138000666650, 2.403324012620d-1,kind=dp)
    exact_roots(11) = cmplx(-1d0,0d0,kind=dp)
    exact_roots(12) = cmplx(-5.968916446451269d-1, 6.629807445770295d-1,kind=dp)
    exact_roots(13) = cmplx(-5.968916446451269d-1,-6.629807445770295d-1,kind=dp)
    exact_roots(14) = cmplx(-2.17526747030511d-1,-1.11445426587329,kind=dp)
    exact_roots(15) = cmplx(-2.17526747030511d-1, 1.11445426587329,kind=dp)
    exact_roots(16) = cmplx(-1.6359826155202d-1, 1.09778064288827,kind=dp)
    exact_roots(17) = cmplx(-1.6359826155202d-1,-1.09778064288827,kind=dp)
    exact_roots(18) = cmplx(-1.225611668766536d-1,-7.4486176661974423d-1,kind=dp)
    exact_roots(19) = cmplx(-1.225611668766536d-1, 7.4486176661974423d-1,kind=dp)
    exact_roots(20) = cmplx(-1.13418655949436d-1,-8.605694725015730d-1,kind=dp)
    exact_roots(21) = cmplx(-1.13418655949436d-1,8.605694725015730d-1,kind=dp)
    exact_roots(22) = cmplx(-1.5570386020902d-2, 1.020497366498289d0,kind=dp)
    exact_roots(23) = cmplx(-1.5570386020902d-2,-1.020497366498289d0,kind=dp)
    exact_roots(24) = cmplx(3.59892739012579001d-1, 6.84762020211812856d-1,kind=dp)
    exact_roots(25) = cmplx(3.59892739012579001d-1,-6.84762020211812856d-1,kind=dp)
    exact_roots(26) = cmplx(3.8900684056977123543d-1,-2.1585065087081910777d-1,kind=dp)
    exact_roots(27) = cmplx(3.8900684056977123543d-1, 2.1585065087081910777d-1,kind=dp)
    exact_roots(28) = cmplx(3.96534570032415023d-1, 6.04181810488988837d-1,kind=dp)
    exact_roots(29) = cmplx(3.96534570032415023d-1,-6.04181810488988837d-1,kind=dp)
    exact_roots(30) = cmplx(4.433256333996235387d-1, 3.729624166628465083d-1,kind=dp)
    exact_roots(31) = cmplx(4.433256333996235387d-1,-3.729624166628465083d-1,kind=dp)
    write(1, '(A)', advance='no') '13, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 14: Mandelbrot 63 MPSolve
    deg = 63
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d0
    p(2) = 1d0
    p(3) = 2d0
    p(4) = 5d0
    p(5) = 14d0
    p(6) = 42d0
    p(7) = 132d0
    p(8) = 365d0
    p(9) = 950d0
    p(10) = 2398d0
    p(11) = 5916d0
    p(12) = 14290d0
    p(13) = 33708d0
    p(14) = 77684d0
    p(15) = 175048d0
    p(16) = 385741d0
    p(17) = 831014d0
    p(18) = 1749654d0
    p(19) = 3598964d0
    p(20) = 7228014d0
    p(21) = 14162220d0
    p(22) = 27049196d0
    p(23) = 50323496d0
    p(24) = 91143114d0
    p(25) = 160617860d0
    p(26) = 275276716d0
    p(27) = 458591432d0
    p(28) = 742179284d0
    p(29) = 1166067016d0
    p(30) = 1777171560d0
    p(31) = 2625062128d0
    p(32) = 3754272037d0
    p(33) = 5193067630d0
    p(34) = 6939692682d0
    p(35) = 8948546308d0
    p(36) = 11120136162d0
    p(37) = 13299362332d0
    p(38) = 15286065700d0
    p(39) = 16859410792d0
    p(40) = 17813777994d0
    p(41) = 17999433372d0
    p(42) = 17357937708d0
    p(43) = 15941684776d0
    p(44) = 13910043524d0
    p(45) = 11500901864d0
    p(46) = 8984070856d0
    p(47) = 6609143792d0
    p(48) = 4562339774d0
    p(49) = 2943492972d0
    p(50) = 1766948340d0
    p(51) = 981900168d0
    p(52) = 502196500d0
    p(53) = 234813592d0
    p(54) = 99582920d0
    p(55) = 37945904d0
    p(56) = 12843980d0
    p(57) = 3807704d0
    p(58) = 971272d0
    p(59) = 208336d0
    p(60) = 36440d0
    p(61) = 4976d0
    p(62) = 496d0
    p(63) = 32d0
    p(64) = 1d0
    exact_roots(1) = cmplx(-1.999095682327018473210d0,0d0,kind=dp)
    exact_roots(2) = cmplx(-1.9918141725491222157325609498622881d0,0d0,kind=dp)
    exact_roots(3) = cmplx(-1.977179587006257387346088520662828616836d0,0d0,kind=dp)
    exact_roots(4) = cmplx(-1.953705894284396245427622199013653238901d0,0d0,kind=dp)
    exact_roots(5) = cmplx(-1.927147709363950262460068188946594278007d0,0d0,kind=dp)
    exact_roots(6) = cmplx(-1.8848035715866817923294780929158396496359d0,0d0,kind=dp)
    exact_roots(7) = cmplx(-1.8323152027512291920848975260425181432293d0,0d0,kind=dp)
    exact_roots(8) = cmplx(-1.76926167027683114607548022863625740038777d0, 5.6919500395600315304900187298015859319654d-2,kind=dp)
    exact_roots(9) = cmplx(-1.76926167027683114607548022863625740038777d0,-5.6919500395600315304900187298015859319654d-2,kind=dp)
    exact_roots(10) = cmplx(-1.674066091474787971565296029172325596206403d0,0d0,kind=dp)
    exact_roots(11) = cmplx(-1.5748891397523009698199655524959742837719482d0,0d0,kind=dp)
    exact_roots(12) = cmplx(-1.408446485740072654917577008805998851928020904d0, &
         -1.36171997304659915684707793608163610038822995d-1,kind=dp)
    exact_roots(13) = cmplx(-1.408446485740072654917577008805998851928020904d0, &
         1.36171997304659915684707793608163610038822995d-1,kind=dp)
    exact_roots(14) = cmplx(-1.29255806103352208716418470636149411998013630326d0, &
         4.3819881608663183712973712432734844004535476504d-1,kind=dp)
    exact_roots(15) = cmplx(-1.29255806103352208716418470636149411998013630326d0, &
         -4.3819881608663183712973712432734844004535476504d-1,kind=dp)
    exact_roots(16) = cmplx(-1.26228728143847254301011194120806575232050489502d0, &
         4.0810432411269038329016065742601506306041169168d-1,kind=dp)
    exact_roots(17) = cmplx(-1.26228728143847254301011194120806575232050489502d0, &
         -4.0810432411269038329016065742601506306041169168d-1,kind=dp)
    exact_roots(18) = cmplx(-1.25273588401203794629581100256433997387062287256d0, &
         -3.4247064788975089386187578687092843396383393805d-1,kind=dp)
    exact_roots(19) = cmplx(-1.25273588401203794629581100256433997387062287256d0, &
         3.4247064788975089386187578687092843396383393805d-1,kind=dp)
    exact_roots(20) = cmplx(-1.02819385245481759930249745596731843328070508279421d0, &
         -3.61376517118561592479460832997830315786692639704085d-1,kind=dp)
    exact_roots(21) = cmplx(-1.02819385245481759930249745596731843328070508279421d0, &
         3.61376517118561592479460832997830315786692639704085d-1,kind=dp)
    exact_roots(22) = cmplx(-6.23532485956252757990016587001026776428072703359878868d-1, &
         6.81064414225239608090835812686561539088332735217609127d-1,kind=dp)
    exact_roots(23) = cmplx(-6.23532485956252757990016587001026776428072703359878868d-1, &
         -6.81064414225239608090835812686561539088332735217609127d-1,kind=dp)
    exact_roots(24) = cmplx(-6.2243629504129358796016350694723840189750985673649588591d-1, &
         4.2487843647562918431157443880525338683545992964599689876d-1,kind=dp)
    exact_roots(25) = cmplx(-6.2243629504129358796016350694723840189750985673649588591d-1, &
         -4.2487843647562918431157443880525338683545992964599689876d-1,kind=dp)
    exact_roots(26) = cmplx(-5.308278048599427289214772971196026578135170949646890946d-1, &
         6.682887255592057714440924655647011851367651843270734380d-1,kind=dp)
    exact_roots(27) = cmplx(-5.308278048599427289214772971196026578135170949646890946d-1, &
         -6.682887255592057714440924655647011851367651843270734380d-1,kind=dp)
    exact_roots(28) = cmplx(-2.72102461488938894219383324518026874585894699621947085d-1, &
         -8.42364690294128145503155708242929569550778268698265965d-1,kind=dp)
    exact_roots(29) = cmplx(-2.72102461488938894219383324518026874585894699621947085d-1, &
         8.42364690294128145503155708242929569550778268698265965d-1,kind=dp)
    exact_roots(30) = cmplx(-2.24915951286740054685326255204118310792682454680693d-1, &
         1.11626015745499183500126825424467009109873946082435d0,kind=dp)
    exact_roots(31) = cmplx(-2.24915951286740054685326255204118310792682454680693d-1, &
         -1.11626015745499183500126825424467009109873946082435d0,kind=dp)
    exact_roots(32) = cmplx(-2.0728383545566641282413385018667121332401155604017d-1, &
         1.11748077249496291137377567312207879579746389236127d0,kind=dp)
    exact_roots(33) = cmplx(-2.0728383545566641282413385018667121332401155604017d-1, &
         -1.11748077249496291137377567312207879579746389236127d0,kind=dp)
    exact_roots(34) = cmplx(-1.7457822113571696945156643266162905020167505710204d-1, &
         1.07142767145403118922964631021955987671322451961088d0,kind=dp)
    exact_roots(35) = cmplx(-1.7457822113571696945156643266162905020167505710204d-1, &
         -1.07142767145403118922964631021955987671322451961088d0,kind=dp)
    exact_roots(36) = cmplx(-1.57516053475965356164335109644674141293297577896685d-1, &
         -1.10900651411360717797175198615475582901468585712356d0,kind=dp) 
    exact_roots(37) = cmplx(-1.57516053475965356164335109644674141293297577896685d-1, &
         1.10900651411360717797175198615475582901468585712356d0,kind=dp)
    exact_roots(38) = cmplx(-1.274999735463630001995395653459879637298616757217284d-1, &
         9.874609094894567922074076807929788675642068522522938d-1,kind=dp)
    exact_roots(39) = cmplx(-1.274999735463630001995395653459879637298616757217284d-1, &
         -9.874609094894567922074076807929788675642068522522938d-1,kind=dp)
    exact_roots(40) = cmplx(-1.42334819203540667677618453136202688358025283954839d-2, &
         -1.0329147752136441093950134026551104360994260360822540d0,kind=dp)
    exact_roots(41) = cmplx(-1.42334819203540667677618453136202688358025283954839d-2, &
         1.0329147752136441093950134026551104360994260360822540d0,kind=dp)
    exact_roots(42) = cmplx(-6.98356849626139181796649107548406610452886379651341d-3, &
         -1.0036038622882895485307049669513531297649273745391915d0,kind=dp)
    exact_roots(43) = cmplx(-6.98356849626139181796649107548406610452886379651341d-3, &
         1.0036038622882895485307049669513531297649273745391915d0,kind=dp)
    exact_roots(44) = cmplx( 1.4895466603687646529815779208794106185666477731693128d-2, &
         -8.481487619084165277193311117832376290806619901265058603d-1,kind=dp)
    exact_roots(45) = cmplx( 1.4895466603687646529815779208794106185666477731693128d-2, &
         8.481487619084165277193311117832376290806619901265058603d-1,kind=dp)
    exact_roots(46) = cmplx( 1.211927861059064863147044434105037593859287800520963579338d-1, &
         6.1061169221075421167538724415035774824319702690063863369691d-1,kind=dp)
    exact_roots(47) = cmplx( 1.211927861059064863147044434105037593859287800520963579338d-1, &
         -6.1061169221075421167538724415035774824319702690063863369691d-1,kind=dp)
    exact_roots(48) = cmplx( 3.52482539722363278193253964052161589243593334212239870706d-1, &
         -6.98337239583330331258141954760484537633150485928512286760d-1,kind=dp)
    exact_roots(49) = cmplx( 3.52482539722363278193253964052161589243593334212239870706d-1, &
         6.98337239583330331258141954760484537633150485928512286760d-1,kind=dp)
    exact_roots(50) = cmplx( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
         -1.4474937132163286474711018201298830556966056842762643026975894d-1,kind=dp)
    exact_roots(51) = cmplx( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
         1.4474937132163286474711018201298830556966056842762643026975894d-1,kind=dp)
    exact_roots(52) = cmplx( 3.76893240379311323690004017968512473363482317941533875341d-1, &
         6.78568693190448141957540792996773280196881194582788907016d-1,kind=dp)
    exact_roots(53) = cmplx( 3.76893240379311323690004017968512473363482317941533875341d-1, &
         -6.78568693190448141957540792996773280196881194582788907016d-1,kind=dp)
    exact_roots(54) = cmplx( 3.865391765961580265082930869043677799284877313516569138807d-1, &
         5.693247113031029032137923571351905081619323911951388853856d-1,kind=dp)
    exact_roots(55) = cmplx( 3.865391765961580265082930869043677799284877313516569138807d-1, &
         -5.693247113031029032137923571351905081619323911951388853856d-1,kind=dp)
    exact_roots(56) = cmplx( 4.12916024722700479197334566382612257174765142865547121703d-1, &
         6.148067601433856949545497204007997358291659758563137777616d-1,kind=dp)
    exact_roots(57) = cmplx( 4.12916024722700479197334566382612257174765142865547121703d-1, &
         -6.148067601433856949545497204007997358291659758563137777616d-1,kind=dp)
    exact_roots(58) = cmplx( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
         2.267599044353486186978765599716989721202321914603899690444951d-1,kind=dp)
    exact_roots(59) = cmplx( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
         -2.267599044353486186978765599716989721202321914603899690444951d-1,kind=dp)
    exact_roots(60) = cmplx( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
         -3.96170128033165002412596877271155937712569079351815707744770d-1,kind=dp)
    exact_roots(61) = cmplx( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
         3.96170128033165002412596877271155937712569079351815707744770d-1,kind=dp)
    exact_roots(62) = cmplx( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
         3.47758700883481983632188723200264206004781117755664551397643d-1,kind=dp)
    exact_roots(63) = cmplx( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
         -3.47758700883481983632188723200264206004781117755664551397643d-1,kind=dp)   
    write(1, '(A)', advance='no') '14, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 15: Jenkins Traub p1(z) with a=1e-8
    deg = 3
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d-16
    p(2) = -1d-16
    p(3) = -1d0
    p(4) = 1d0
    exact_roots(1) = cmplx(1d-8,0d0,kind=dp)
    exact_roots(2) = cmplx(-1d-8,0d0,kind=dp)
    exact_roots(3) = cmplx(1d0,0d0,kind=dp)
    write(1, '(A)', advance='no') '15, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 16: Jenkins Traub p1(z) with a=1e-15
    deg = 3
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d-30
    p(2) = -1d-30
    p(3) = -1d0
    p(4) = 1d0
    exact_roots(1) = cmplx(1d-15,0d0,kind=dp)
    exact_roots(2) = cmplx(-1d-15,0d0,kind=dp)
    exact_roots(3) = cmplx(1d0,0d0,kind=dp)
    write(1, '(A)', advance='no') '16, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 17: Jenkins Traub p1(z) with a=1e+8
    deg = 3
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d+16
    p(2) = -1d+16
    p(3) = -1d0
    p(4) = 1d0
    exact_roots(1) = cmplx(1d+8,0d0,kind=dp)
    exact_roots(2) = cmplx(-1d+8,0d0,kind=dp)
    exact_roots(3) = cmplx(1d0,0d0,kind=dp)
    write(1, '(A)', advance='no') '17, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 18: Jenkins Traub p1(z) with a=1e+15
    deg = 3
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d+30
    p(2) = -1d+30
    p(3) = -1d0
    p(4) = 1d0
    exact_roots(1) = cmplx(1d+15,0d0,kind=dp)
    exact_roots(2) = cmplx(-1d+15,0d0,kind=dp)
    exact_roots(3) = cmplx(1d0,0d0,kind=dp)
    write(1, '(A)', advance='no') '18, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 19: Jenkins Traub p3(z) deg 10
    deg = 10
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1.0000000000000000000d-55
    p(2) = -1.1111111110000000000d-45
    p(3) = 1.1223344544332211000d-36
    p(4) = -1.1234579011109875432d-28
    p(5) = 1.1235701457797754097d-21
    p(6) = -1.1235802580122097520d-15
    p(7) = 1.1235701457797754097d-10
    p(8) = -1.1234579011109875432d-6
    p(9) = 1.1223344544332211000d-3
    p(10) = -1.1111111110000000000d-1
    p(11) = 1.0000000000000000000d0
    exact_roots(1) = cmplx(1d-1,0d0,kind=dp)
    exact_roots(2) = cmplx(1d-2,0d0,kind=dp)
    exact_roots(3) = cmplx(1d-3,0d0,kind=dp)
    exact_roots(4) = cmplx(1d-4,0d0,kind=dp)
    exact_roots(5) = cmplx(1d-5,0d0,kind=dp)
    exact_roots(6) = cmplx(1d-6,0d0,kind=dp)
    exact_roots(7) = cmplx(1d-7,0d0,kind=dp)
    exact_roots(8) = cmplx(1d-8,0d0,kind=dp)
    exact_roots(9) = cmplx(1d-9,0d0,kind=dp)
    exact_roots(10) = cmplx(1d-10,0d0,kind=dp)
    write(1, '(A)', advance='no') '19, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 20: Jenkins Traub p3(z) deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1.0000000000000000000d-210
    p(2) = -1.1111111111111111111d-190
    p(3) = 1.1223344556677890010d-171
    p(4) = -1.1234579135813703702d-153
    p(5) = 1.1235702706084312021d-136
    p(6) = -1.1235815064234953247d-120
    p(7) = 1.1235826300061242073d-105
    p(8) = -1.1235827423643872079d-91
    p(9) = 1.1235827536001023856d-78
    p(10) = -1.1235827547225615576d-66
    p(11) = 1.1235827548236840055d-55
    p(12) = -1.1235827547225615576d-45
    p(13) = 1.1235827536001023856d-36
    p(14) = -1.1235827423643872079d-28
    p(15) = 1.1235826300061242073d-21
    p(16) = -1.1235815064234953247d-15
    p(17) = 1.1235702706084312021d-10
    p(18) = -1.1234579135813703702d-6
    p(19) = 1.1223344556677890010d-3
    p(20) = -1.1111111111111111111d-1
    p(21) = 1.0000000000000000000d0
    exact_roots(1) = cmplx(1d-1,0d0,kind=dp)
    exact_roots(2) = cmplx(1d-2,0d0,kind=dp)
    exact_roots(3) = cmplx(1d-3,0d0,kind=dp)
    exact_roots(4) = cmplx(1d-4,0d0,kind=dp)
    exact_roots(5) = cmplx(1d-5,0d0,kind=dp)
    exact_roots(6) = cmplx(1d-6,0d0,kind=dp)
    exact_roots(7) = cmplx(1d-7,0d0,kind=dp)
    exact_roots(8) = cmplx(1d-8,0d0,kind=dp)
    exact_roots(9) = cmplx(1d-9,0d0,kind=dp)
    exact_roots(10) = cmplx(1d-10,0d0,kind=dp)
    exact_roots(11) = cmplx(1d-11,0d0,kind=dp)
    exact_roots(12) = cmplx(1d-12,0d0,kind=dp)
    exact_roots(13) = cmplx(1d-13,0d0,kind=dp)
    exact_roots(14) = cmplx(1d-14,0d0,kind=dp)
    exact_roots(15) = cmplx(1d-15,0d0,kind=dp)
    exact_roots(16) = cmplx(1d-16,0d0,kind=dp)
    exact_roots(17) = cmplx(1d-17,0d0,kind=dp)
    exact_roots(18) = cmplx(1d-18,0d0,kind=dp)
    exact_roots(19) = cmplx(1d-19,0d0,kind=dp)
    exact_roots(20) = cmplx(1d-20,0d0,kind=dp)
    write(1, '(A)', advance='no') '20, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 21: Jenkins Traub p4(z)
    deg = 6
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 2.1d-4
    p(2) = -7.37d-3
    p(3) = 9.69d-2
    p(4) = -5.86d-1
    p(5) = 1.64d0
    p(6) = -2.1d0
    p(7) = 1d0
    exact_roots(1) = 1d-1
    exact_roots(2) = 1d-1
    exact_roots(3) = 1d-1
    exact_roots(4) = 5d-1
    exact_roots(5) = 6d-1
    exact_roots(6) = 7d-1
    write(1, '(A)', advance='no') '21, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 22: Jenkins Traub p5(z)
    deg = 10
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 2.88d-8
    p(2) = -1.848d-6
    p(3) = 5.204d-5
    p(4) = -8.458d-4
    p(5) = 8.777d-3
    p(6) = -6.072d-2
    p(7) = 2.835d-1
    p(8) = -8.820d-1
    p(9) = 1.750d0
    p(10) = -2d0
    p(11) = 1d0
    exact_roots(1) = 1d-1
    exact_roots(2) = 1d-1
    exact_roots(3) = 1d-1
    exact_roots(4) = 1d-1
    exact_roots(5) = 2d-1
    exact_roots(6) = 2d-1
    exact_roots(7) = 2d-1
    exact_roots(8) = 3d-1
    exact_roots(9) = 3d-1
    exact_roots(10) = 4d-1
    write(1, '(A)', advance='no') '22, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 23: Jenkins Traub p6(z)
    deg = 5
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -9.9900798978020040d-2
    p(2) = 1.3987105877382204d0
    p(3) = -4.5967287785602000d0
    p(4) = 6.3969289898000000d0
    p(5) = -4.0990100000000000d0
    p(6) = 1.0000000000000000d0
    exact_roots(1) = 1d-1
    exact_roots(2) = 1001d-3
    exact_roots(3) = 998d-3
    exact_roots(4) = 100002d-5
    exact_roots(5) = 99999d-5
    write(1, '(A)', advance='no') '23, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 24: Jenkins Traub p7(z) with a=0
    deg = 7
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1.000000d-7
    p(2) = 1.131100d-4
    p(3) = -1.345431d-2
    p(4) = 3.477743d-1
    p(5) = -3.477743d0
    p(6) = 1.345431d1
    p(7) = -1.131100d1
    p(8) = 1.000000d0
    exact_roots(1) = 1d-3
    exact_roots(2) = 1d-2
    exact_roots(3) = 1d-1
    exact_roots(4) = 1d-1
    exact_roots(5) = 1d-1
    exact_roots(6) = 1d0
    exact_roots(7) = 1d+1
    write(1, '(A)', advance='no') '24, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 25: Jenkins Traub p7(z) with a=10^(-10)
    deg = 7
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1.0000000000000000010d-7
    p(2) = 1.1311000000000000011d-4
    p(3) = -1.3454310000000000011d-2
    p(4) = 3.4777430000000000011d-1
    p(5) = -3.4777430000000000001d0
    p(6) = 1.3454310000000000000d+1
    p(7) = -1.1311000000000000000d+1
    p(8) = 1.0000000000000000000d0
    exact_roots(1) = 1d-3
    exact_roots(2) = 1d-2
    exact_roots(3) = 1d-1
    exact_roots(4) = cmplx(1d-1,1d-10,kind=dp)
    exact_roots(5) = cmplx(1d-1,-1d-10,kind=dp)
    exact_roots(6) = 1d0
    exact_roots(7) = 1d+1
    write(1, '(A)', advance='no') '25, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 26: Jenkins Traub p7(z) with a=10^(-6)
    deg = 7
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1.0000000001000000d-7
    p(2) = 1.1311000001111100d-4
    p(3) = -1.3454310001122211d-2
    p(4) = 3.4777430001122211d-1
    p(5) = -3.4777430000111110d0
    p(6) = 1.3454310000001000d+1
    p(7) = -1.1311000000000000d+1
    p(8) = 1.0000000000000000d0
    exact_roots(1) = 1d-3
    exact_roots(2) = 1d-2
    exact_roots(3) = 1d-1
    exact_roots(4) = cmplx(1d-1,1d-6,kind=dp)
    exact_roots(5) = cmplx(1d-1,-1d-6,kind=dp)
    exact_roots(6) = 1d0
    exact_roots(7) = 1d+1
    write(1, '(A)', advance='no') '26, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 27: Jenkins Traub p8(z)
    deg = 5
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1d0
    p(2) = 5d0
    p(3) = 10d0
    p(4) = 10d0
    p(5) = 5d0
    p(6) = 1d0
    exact_roots(1) = -1d0
    exact_roots(2) = -1d0
    exact_roots(3) = -1d0
    exact_roots(4) = -1d0
    exact_roots(5) = -1d0
    write(1, '(A)', advance='no') '27, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 28: Jenkins Traub p9(z)
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1d0
    p(2) = 0d0
    p(3) = 0d0
    p(4) = 0d0
    p(5) = 0d0
    p(6) = 0d0
    p(7) = 0d0
    p(8) = 0d0
    p(9) = 0d0
    p(10) = 0d0
    p(11) = 1d+20
    p(12) = 0d0
    p(13) = 0d0
    p(14) = 0d0
    p(15) = 0d0
    p(16) = 0d0
    p(17) = 0d0
    p(18) = 0d0
    p(19) = 0d0
    p(20) = 0d0
    p(21) = 1d0
    exact_roots(1) = -1d-2
    exact_roots(2) = 1d-2
    exact_roots(3) = cmplx(0d0,-1d+2,kind=dp)
    exact_roots(4) = cmplx(0d0,1d+2,kind=dp)
    exact_roots(5) = cmplx(-95.105651629515357212d0,+30.901699437494742410d0,kind=dp)
    exact_roots(6) = cmplx(-95.105651629515357212d0,-30.901699437494742410d0,kind=dp)
    exact_roots(7) = cmplx(-0.0080901699437494742410d0,+0.0058778525229247312917d0,kind=dp)
    exact_roots(8) = cmplx(-0.0080901699437494742410d0,-0.0058778525229247312917d0,kind=dp)
    exact_roots(9) = cmplx(-58.778525229247312917d0,+80.901699437494742410d0,kind=dp)
    exact_roots(10) = cmplx(-58.778525229247312917d0,-80.901699437494742410d0,kind=dp)
    exact_roots(11) = cmplx(-0.0030901699437494742410d0,+0.0095105651629515357212d0,kind=dp)
    exact_roots(12) = cmplx(-0.0030901699437494742410d0,-0.0095105651629515357212d0,kind=dp)    
    exact_roots(13) = cmplx(+95.105651629515357212d0,+30.901699437494742410d0,kind=dp)
    exact_roots(14) = cmplx(+95.105651629515357212d0,-30.901699437494742410d0,kind=dp)
    exact_roots(15) = cmplx(+0.0080901699437494742410d0,+0.0058778525229247312917d0,kind=dp)
    exact_roots(16) = cmplx(+0.0080901699437494742410d0,-0.0058778525229247312917d0,kind=dp)
    exact_roots(17) = cmplx(+58.778525229247312917d0,+80.901699437494742410d0,kind=dp)
    exact_roots(18) = cmplx(+58.778525229247312917d0,-80.901699437494742410d0,kind=dp)
    exact_roots(19) = cmplx(+0.0030901699437494742410d0,+0.0095105651629515357212d0,kind=dp)
    exact_roots(20) = cmplx(+0.0030901699437494742410d0,-0.0095105651629515357212d0,kind=dp)
    write(1, '(A)', advance='no') '28, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 29: Jenkins Traub p10(z) with a=1e+3
    deg = 3
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1d0
    p(2) = 1001001d-3
    p(3) = -1001001d-3
    p(4) = 1d0
    exact_roots(1) = 1d+3
    exact_roots(2) = 1d0
    exact_roots(3) = 1d-3
    write(1, '(A)', advance='no') '29, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 30: Jenkins Traub p10(z) with a=1e+6
    deg = 3
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1.000000000000d0
    p(2) = 1.000001000001d+6
    p(3) = -1.000001000001d+6
    p(4) = 1.000000000000d0
    exact_roots(1) = 1d+6
    exact_roots(2) = 1d0
    exact_roots(3) = 1d-6
    write(1, '(A)', advance='no') '30, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 31: Jenkins Traub p10(z) with a=1e+9
    deg = 3
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1.0000000000000000000d0
    p(2) = 1.0000000010000000010d+9
    p(3) = -1.0000000010000000010d+9
    p(4) = 1.0000000000000000000d0
    exact_roots(1) = 1d+9
    exact_roots(2) = 1d0
    exact_roots(3) = 1d-9
    write(1, '(A)', advance='no') '31, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 32: Jenkins Traub p11(z) with m=15
    deg = 60
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 4.2391158275216251d-2
    p(2) = -9.4584738186193285d-2
    p(3) = 1.0794581538155193d-2
    p(4) = 0.16078797907622777d0    
    p(5) = -0.14845902865058211d0    
    p(6) = -0.10102576217529596d0   
    p(7) = 0.21209884906475510d0  
    p(8) = 3.8416291071329871d-2
    p(9) = -0.24858045633123330d0     
    p(10) = 2.0307451027496867d-2
    p(11) = 0.27253121225404131d0    
    p(12) = -7.5673575543483473d-2
    p(13) = -0.29022672239698227d0    
    p(14) = 0.12893546540368780d0  
    p(15) = 0.30496711848391356d0   
    p(16) = -0.18123908396656324d0    
    p(17) = -0.31881142664054229d0   
    p(18) = 0.23353202157910169d0  
    p(19) = 0.33327316630240250d0   
    p(20) = -0.28658432122905075d0    
    p(21) = -0.34964976585288843d0   
    p(22) = 0.34101274057937758d0  
    p(23) = 0.36920471966525736d0   
    p(24) = -0.39728789823384797d0    
    p(25) = -0.39328601327814472d0   
    p(26) = 0.45571962435254554d0  
    p(27) = 0.42341663675939667d0   
    p(28) = -0.51641706563343759d0    
    p(29) = -0.46137376166906691d0   
    p(30) = 0.57921767479125386d0  
    p(31) = 0.50926408481343588d0   
    p(32) = -0.64357519421250398d0    
    p(33) = -0.56959723662847717d0   
    p(34) = 0.70839103653420699d0  
    p(35) = 0.64535381307635498d0   
    p(36) = -0.77176518544352035d0    
    p(37) = -0.74003701874365146d0   
    p(38) = 0.83063030145887728d0  
    p(39) = 0.85768372384334957d0   
    p(40) = -0.88021348963651835d0    
    p(41) = -1.0027857350532181d0  
    p(42) = 0.91323978476130963d0   
    p(43) = 1.1800223533802476d0  
    p(44) = -0.91874251967352283d0     
    p(45) = -1.3936019538936961d0  
    p(46) = 0.88026658614537412d0   
    p(47) = 1.6457841113586349d0  
    p(48) = -0.77312502631224089d0     
    p(49) = -1.9336248882093454d0  
    p(50) = 0.56019157805003328d0   
    p(51) = 2.2416409891507394d0  
    p(52) = -0.18559348051633595d0     
    p(53) = -2.5242465611040741d0  
    p(54) = -0.43344871139455649d0    
    p(55) = 2.6589984570372667d0 
    p(56) = 1.4072439804450414d0   
    p(57) = -2.2977425637976268d0    
    p(58) = -2.7650680357818582d0   
    p(59) = 0.20626025335611423d0   
    p(60) = 2.0081136687728178d0
    p(61) = 1d0
    do j=-14,0
       exact_roots(15+j) = 0.9d0*exp(cmplx(0d0,j*pi/30d0,kind=dp))
    end do
    do j=1,15
       exact_roots(15+j) = conjg(exact_roots(16-j))
    end do
    do j=16,30
       exact_roots(15+j) = exp(cmplx(0d0,j*pi/30d0,kind=dp))
    end do
    do j=31,45
       exact_roots(15+j) = conjg(exact_roots(76-j))
    end do
    write(1, '(A)', advance='no') '32, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! NAG
    allocate(a(2,0:deg),w(4*(deg+1)),z(2,deg))
    allocate(zeros(deg))
    do j=0,deg
        a(1,j) = dble(p(deg+1-j))
        a(2,j) = aimag(p(deg+1-j))
    end do
    ifail = 0
    call c02aff(a,deg,scal,z,w,ifail)
    zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(a, w, z, zeros)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    close(1)
contains
    !************************************************
    !                       sort                    *
    !************************************************
    ! The roots are sorted with respect to exact_roots. 
    ! For each i, roots(i) is the root approximation 
    ! that is closest to exact_roots(i).
    !************************************************
    subroutine sort(roots, exact_roots, deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
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
            temp = roots(i)
            roots(i) = roots(k)
            roots(k) = temp
        end do
    end subroutine sort
end program spec_poly
