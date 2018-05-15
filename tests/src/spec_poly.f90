!********************************************************************************
!   SPEC_POLY: Compare FPML against Polzeros and AMVW for special polynomials
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 21 March 2018
!********************************************************************************
! The accuracy of FPML is compared against Polzeros and AMVW on 14 notoriously
! difficult special polynomials.  
!********************************************************************************
program spec_poly
    use fpml
    use poly_zeroes
    implicit none
    ! testing variables
    integer                                     :: deg, j
    real(kind=dp), parameter                    :: pi = 3.141592653589793D0
    real(kind=dp), dimension(:), allocatable    :: err
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    real(kind=dp), dimension(:),    allocatable :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! Polzeros variables
    integer                                     :: iter
    integer, parameter                          :: nitmax=30
    real(kind=dp), parameter                    :: small=tiny(1.0D0), big=huge(1.0D0)
    logical, dimension(:), allocatable          :: h
    real(kind=dp), dimension(:), allocatable    :: radius
    complex(kind=dp), dimension(:), allocatable :: zeros
    ! AMVW variables
    integer                                     :: flag
    real(kind=dp), dimension(:), allocatable    :: residuals
    complex(kind=dp), dimension(:), allocatable :: coeffs, eigs
    
    ! Testing: special polynomials
    call init_random_seed()
    open(unit=1,file="data_files/spec_poly.dat")
    write(1,'(A)') 'Poly No., FPML, Polzeros, AMVW'
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1)/p(deg+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (dble(p(deg-j+1)), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 32: Jenkins Traub p11(z) with m=15
    deg = 60
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 4.2391158275216203514d-2
    p(2) = -9.4584738186193050855d-2
    p(3) = 1.07945815381548962735d-2
    p(4) = 1.6078797907622771809d-1
    p(5) = -1.4845902865058159236d-1
    p(6) = -1.01025762175296532794d-1
    p(7) = 2.1209884906475519856d-1
    p(8) = 3.8416291071330318748d-2
    p(9) = -2.4858045633123371111d-1
    p(10) = 2.0307451027496707512d-2
    p(11) = 2.7253121225404153306d-1
    p(12) = -7.5673575543483108178d-2
    p(13) = -2.9022672239698259232d-1
    p(14) = 1.28935465403687304842d-1
    p(15) = 3.0496711848391418142d-1
    p(16) = -1.8123908396656292076d-1
    p(17) = -3.1881142664054288599d-1
    p(18) = 2.3353202157910117428d-1
    p(19) = 3.3327316630240310058d-1
    p(20) = -2.8658432122905003452d-1
    p(21) = -3.4964976585288927202d-1
    p(22) = 3.4101274057937680974d-1
    p(23) = 3.6920471966525847300d-1
    p(24) = -3.9728789823384720065d-1
    p(25) = -3.9328601327814619351d-1
    p(26) = 4.5571962435254494354d-1
    p(27) = 4.2341663675939806950d-1
    p(28) = -5.1641706563343710483d-1
    p(29) = -4.6137376166906815362d-1
    p(30) = 5.7921767479125346222d-1
    p(31) = 5.0926408481343738485d-1
    p(32) = -6.4357519421250384691d-1
    p(33) = -5.6959723662847920200d-1
    p(34) = 7.0839103653420727685d-1
    p(35) = 6.4535381307635736854d-1
    p(36) = -7.7176518544352138654d-1
    p(37) = -7.4003701874365393997d-1
    p(38) = 8.3063030145887878565d-1
    p(39) = 8.5768372384335260519d-1
    p(40) = -8.8021348963651948138d-1
    p(41) = -1.00278573505322181238d0
    p(42) = 9.1323978476131014358d-1
    p(43) = 1.18002235338025109601d0
    p(44) = -9.1874251967352285198d-1
    p(45) = -1.39360195389369934621d0
    p(46) = 8.8026658614537399624d-1
    p(47) = 1.6457841113586375191d0
    p(48) = -7.7312502631224041499d-1
    p(49) = -1.9336248882093467243d0
    p(50) = 5.6019157805003299369d-1
    p(51) = 2.2416409891507395007d0
    p(52) = -1.8559348051633771754d-1
    p(53) = -2.5242465611040740194d0
    p(54) = -4.3344871139455217442d-1
    p(55) = 2.6589984570372695878d0
    p(56) = 1.40724398044503771124d0
    p(57) = -2.2977425637976343484d0
    p(58) = -2.7650680357818599628d0
    p(59) = 2.0626025335611973886d-1
    p(60) = 2.0081136687728211063d0
    p(61) = 1.0000000000000000000d0
    do j=-14,0
       exact_roots(15+j) = 0.9d0*cmplx(cos(j*pi/30d0),sin(j*pi/30d0),kind=dp)
    end do
    do j=1,15
       exact_roots(15+j) = conjg(exact_roots(16-j))
    end do
    do j=16,30
       exact_roots(15+j) = cmplx(cos(j*pi/30d0),sin(j*pi/30d0),kind=dp)
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
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 33: Jenkins Traub p11(z) with m=20
    deg = 80
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1.4780882941434592332d-2
    p(2) = -4.3442163898274625504d-2
    p(3) = 3.0811061677962020381d-2
    p(4) = 4.8174145573905845167d-2
    p(5) = -9.9895422090259650436d-2
    p(6) = 2.0367593752690566382d-2
    p(7) = 1.02730704300754634317d-1
    p(8) = -7.2421959976329501386d-2
    p(9) = -8.6179510586397895483d-2
    p(10) = 1.11130623035162433422d-1
    p(11) = 6.2288430530155258569d-2
    p(12) = -1.40861180023385154071d-1
    p(13) = -3.5124501490762462669d-2
    p(14) = 1.6457808423330592019d-1
    p(15) = 6.2141637784362712705d-3
    p(16) = -1.8424846472346489127d-1
    p(17) = 2.3887246006411263876d-2
    p(18) = 2.0123660776207915149d-1
    p(19) = -5.5042151404477955885d-2
    p(20) = -2.1655010742652691907d-1
    p(21) = 8.7308247479862632512d-2
    p(22) = 2.3098952256062750832d-1
    p(23) = -1.20839896048866662802d-1
    p(24) = -2.4524316316862245476d-1
    p(25) = 1.5584030932109201157d-1
    p(26) = 2.5995126187880752178d-1
    p(27) = -1.9253542069511053711d-1
    p(28) = -2.7575320187036790325d-1
    p(29) = 2.3115653041457324242d-1
    p(30) = 2.9332567468690245068d-1
    p(31) = -2.7192522443725521915d-1
    p(32) = -3.1341647125406680241d-1
    p(33) = 3.1503664809187655702d-1
    p(34) = 3.3687681253896948062d-1
    p(35) = -3.6063813116693822381d-1
    p(36) = -3.6469402480999342566d-1
    p(37) = 4.0880023982723262898d-1
    p(38) = 3.9802557620913000863d-1
    p(39) = -4.5947689308160295174d-1
    p(40) = -4.3823476903404861240d-1
    p(41) = 5.1245032822053625773d-1
    p(42) = 4.8692752114894290267d-1
    p(43) = -5.6725542355753450833d-1
    p(44) = -5.4598844473131688426d-1
    p(45) = 6.2307611618233901689d-1
    p(46) = 6.1761253333670921719d-1
    p(47) = -6.7860426870892201356d-1
    p(48) = -7.0432572851500706072d-1
    p(49) = 7.3184818906851594346d-1
    p(50) = 8.0898269490870113069d-1
    p(51) = -7.7987392727599626297d-1
    p(52) = -9.3472202131370181256d-1
    p(53) = 8.1845735150340386619d-1
    p(54) = 1.08484562323120574344d0
    p(55) = -8.4161895105779053806d-1
    p(56) = -1.26256657697771489173d0
    p(57) = 8.4100707730363691321d-1
    p(58) = 1.4705312179546839242d0
    p(59) = -8.0509137325098417157d-1
    p(60) = -1.7099546866515839107d0
    p(61) = 7.1813332727313867947d-1
    p(62) = 1.9790907332021869408d0
    p(63) = -5.5893356802509694516d-1
    p(64) = -2.2705405932583590595d0
    p(65) = 2.9946485119551650496d-1
    p(66) = 2.5664992503440600430d0
    p(67) = 9.6178378249591328894d-2
    p(68) = -2.8302464072142201120d0
    p(69) = -6.7115047960723573009d-1
    p(70) = 2.9906015117111130672d0
    p(71) = 1.4693731680026281359d0
    p(72) = -2.9128354841688884049d0
    p(73) = -2.5098266205260663539d0
    p(74) = 2.3435135157927505630d0
    p(75) = 3.6936432309637431823d0
    p(76) = -8.1367672571926600698d-1
    p(77) = -4.4341996816502822389d0
    p(78) = -2.3759711962084459258d0
    p(79) = 1.6884620531827973326d0
    p(80) = 2.6451699579357078932d0
    p(81) = 1.0000000000000000000d0
    do j=-19,0
       exact_roots(20+j) = 0.9d0*cmplx(cos(j*pi/40d0),sin(j*pi/40d0),kind=dp)
    end do
    do j=1,20
       exact_roots(20+j) = conjg(exact_roots(21-j))
    end do
    do j=21,40
       exact_roots(20+j) = cmplx(cos(j*pi/40d0),sin(j*pi/40d0),kind=dp)
    end do
    do j=41,60
       exact_roots(20+j) = conjg(exact_roots(101-j))
    end do
    write(1, '(A)', advance='no') '33, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    ! Poly 34: Jenkins Traub p11(z) with m=25
    deg = 100
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 5.1537752073201133104d-3
    p(2) = -1.8794395712668274115d-2
    p(3) = 2.2752526258720053689d-2
    p(4) = 6.4581638228381859949d-3
    p(5) = -4.8039958324421547026d-2
    p(6) = 4.0240863451162136221d-2
    p(7) = 2.6095249860867401354d-2
    p(8) = -6.1946632433651940035d-2
    p(9) = 2.1840603116015807165d-3
    p(10) = 6.8695720449026285311d-2
    p(11) = -2.9149094869278078081d-2
    p(12) = -6.6692624627556374514d-2
    p(13) = 5.3329286518131136703d-2
    p(14) = 5.9263725735899408512d-2
    p(15) = -7.4740169736381119245d-2
    p(16) = -4.8233685965723387626d-2
    p(17) = 9.3762959225871958187d-2
    p(18) = 3.4630923115122010239d-2
    p(19) = -1.10837738134972347479d-1
    p(20) = -1.9040148604704520169d-2
    p(21) = 1.26380197864395981459d-1
    p(22) = 1.7857435712930584793d-3
    p(23) = -1.40766067120767116344d-1
    p(24) = 1.6968338564966387623d-2
    p(25) = 1.5433651367071499382d-1
    p(26) = -3.7161740664692577621d-2
    p(27) = -1.6740956506441383716d-1
    p(28) = 5.8804027003499456932d-2
    p(29) = 1.8029277046097346512d-1
    p(30) = -8.1953843480833198185d-2
    p(31) = -1.9329577904675050835d-1
    p(32) = 1.06705273424908972793d-1
    p(33) = 2.0674270594407362717d-1
    p(34) = -1.33178116826723790942d-1
    p(35) = -2.2098455685728492470d-1
    p(36) = 1.6151004182819932570d-1
    p(37) = 2.3641209704390528239d-1
    p(38) = -1.9184917895817060259d-1
    p(39) = -2.5346956105423180649d-1
    p(40) = 2.2434601565956576337d-1
    p(41) = 2.7266957071909902168d-1
    p(42) = -2.5914353146859951616d-1
    p(43) = -2.9460957752704096511d-1
    p(44) = 2.9636445601061668842d-1
    p(45) = 3.1999006898082906000d-1
    p(46) = -3.3609436040764361776d-1
    p(47) = -3.4963466223155718933d-1
    p(48) = 3.7835901515858605506d-1
    p(49) = 3.8451202746001286928d-1
    p(50) = -4.2309405896374659512d-1
    p(51) = -4.2575930076595750118d-1
    p(52) = 4.7010450995971843902d-1
    p(53) = 4.7470620674075662874d-1
    p(54) = -5.1901099473057072024d-1
    p(55) = -5.3289843351860568408d-1
    p(56) = 5.6917875054216602781d-1
    p(57) = 6.0211776844622274157d-1
    p(58) = -6.1962445504166280069d-1
    p(59) = -6.8439493341906568241d-1
    p(60) = 6.6889475086228471555d-1
    p(61) = 7.8200869156377478495d-1
    p(62) = -7.1490898795279670687d-1
    p(63) = -8.9746123657035980720d-1
    p(64) = 7.5475730001174195500d-1
    p(65) = 1.03341452919741634319d0
    p(66) = -7.8444389607781890954d-1
    p(67) = -1.19256421590425444021d0
    p(68) = 7.9856488479337281444d-1
    p(69) = 1.37741569200646230553d0
    p(70) = -7.8991107631503704308d-1
    p(71) = -1.5899086851641703197d0
    p(72) = 7.4899104928059258771d-1
    p(73) = 1.8308092781904491147d0
    p(74) = -6.6348231488955776107d-1
    p(75) = -2.0987467737915542391d0
    p(76) = 5.1764653616279144143d-1
    p(77) = 2.3887100692541683076d0
    p(78) = -2.9180421854838802418d-1
    p(79) = -2.6897239661998793654d0
    p(80) = -3.7912840308102096606d-2
    p(81) = 2.9812867354059439601d0
    p(82) = 4.9905974577397915717d-1
    p(83) = -3.2279541138860647541d0
    p(84) = -1.12062745072342458261d0
    p(85) = 3.3712113823912271193d0
    p(86) = 1.9269152823109811647d0
    p(87) = -3.3175978429932854817d0
    p(88) = -2.9229144316109087713d0
    p(89) = 2.9224723753541657012d0
    p(90) = 4.0608779869665237091d0
    p(91) = -1.9720807603154555627d0
    p(92) = -5.1640066820857745174d0
    p(93) = 1.8242284752185294546d-1
    p(94) = 5.7489667411904366016d0
    p(95) = 2.6908596365654123217d0
    p(96) = -4.6105673032705917658d0
    p(97) = -6.1157142849158914545d0
    p(98) = -9.1350538924594045766d-1
    p(99) = 3.5759313373596545128d0
    p(100) = 3.2820515953773958039d0
    p(101) = 1.0000000000000000000d0
    do j=-24,0
       exact_roots(25+j) = 0.9d0*cmplx(cos(j*pi/50d0),sin(j*pi/50d0),kind=dp)
    end do
    do j=1,25
       exact_roots(25+j) = conjg(exact_roots(26-j))
    end do
    do j=26,50
       exact_roots(25+j) = cmplx(cos(j*pi/50d0),sin(j*pi/50d0),kind=dp)
    end do
    do j=51,75
       exact_roots(25+j) = conjg(exact_roots(126-j))
    end do
    write(1, '(A)', advance='no') '34, '
    ! FPML
    allocate(roots(deg), berr(deg), cond(deg))
    call main(p, deg, roots, berr, cond)
    call sort(roots, exact_roots, deg)
    err = (/ (abs(roots(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(roots, berr, cond)
    ! Polzeros
    allocate(zeros(deg), radius(deg), h(deg+1))
    call polzeros(deg, p, eps, big, small, nitmax, zeros, radius, h, iter)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)', advance='no') sum(err)/deg
    write(1, '(A)', advance='no') ','
    deallocate(zeros, radius, h)
    ! AMVW
    allocate(coeffs(deg+1), eigs(deg), residuals(deg))
    coeffs = (/ (p(deg-j+1), j=0,deg)/)
    call z_poly_roots(deg, coeffs, eigs, residuals, flag)
    call sort(eigs, exact_roots, deg)
    err = (/ (abs(eigs(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(coeffs, eigs, residuals)
    ! deallocate polynomial arrays
    deallocate(exact_roots, p, err)
    close(1)
contains
    !************************************************
    !                   init_random_seed            *
    !************************************************
    ! Initiate random seed using system_clock. This
    ! seed is then available for the random number
    ! generator in random_number for the life of
    ! the program.
    !************************************************
    subroutine init_random_seed()
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
    end subroutine init_random_seed
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