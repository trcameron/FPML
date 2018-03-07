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
    integer, dimension(:), allocatable          :: its
    real(kind=dp), dimension(:), allocatable    :: reigs, ieigs, rcoeffs, icoeffs
    
    ! testing specific polynomials
    call init_random_seed()
    open(unit=1,file="data_files/spec_poly.dat")
    write(1,'(A)') 'Poly No., FPML, Polzeros, AMVW'
    ! Poly 1: Wilkinson deg 10
    deg = 10
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 3628800D0
    p(2) = -10628640D0
    p(3) = 12753576D0
    p(4) = -8409500D0
    p(5) = 3416930D0
    p(6) = -902055D0
    p(7) = 157773D0
    p(8) = -18150D0
    p(9) = 1320D0
    p(10) = -55D0
    p(11) = 1D0
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 2: Wilkinson deg 15
    deg = 15
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1307674368000D0
    p(2) = 4339163001600D0
    p(3) = -6165817614720D0
    p(4) = 5056995703824D0
    p(5) = -2706813345600D0
    p(6) = 1009672107080D0
    p(7) = -272803210680D0
    p(8) = 54631129553D0
    p(9) = -8207628000D0
    p(10) = 928095740D0
    p(11) = -78558480D0
    p(12) = 4899622D0
    p(13) = -218400D0
    p(14) = 6580D0
    p(15) = -120D0
    p(16) = 1D0
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 3: Wilkinson deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 2432902008176640000D0
    p(2) = -8752948036761600000D0
    p(3) = 13803759753640704000D0
    p(4) = -12870931245150988800D0
    p(5) = 8037811822645051776D0
    p(6) = -3599979517947607200D0
    p(7) = 1206647803780373360D0
    p(8) = -311333643161390640D0
    p(9) = 63030812099294896D0
    p(10) = -10142299865511450D0
    p(11) = 1307535010540395D0
    p(12) = -135585182899530D0
    p(13) = 11310276995381D0
    p(14) = -756111184500D0
    p(15) = 40171771630D0
    p(16) = -1672280820D0
    p(17) = 53327946D0
    p(18) = -1256850D0
    p(19) = 20615D0
    p(20) = -210D0
    p(21) = 1D0
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 4: scaled and shifted Wilkinson deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -758069338497471D0/160000000000000000D0
    p(2) = -1899923154129D0/400000000000000D0
    p(3) = 2280736816325919D0/4000000000000000D0
    p(4) = 11456073304317D0/20000000000000D0
    p(5) = -102724049585427219D0/10000000000000000D0
    p(6) = -65260688438889D0/6250000000000D0
    p(7) = 163007395518693D0/2500000000000D0
    p(8) = 4248959006581D0/62500000000D0
    p(9) = -94957438655047D0/500000000000D0
    p(10) = -259286608191D0/1250000000D0
    p(11) = 141135792837D0/500000000D0
    p(12) = 4187079039D0/12500000D0
    p(13) = -10834566327D0/50000000D0
    p(14) = -18824117D0/62500D0
    p(15) = 1867909D0/25000D0
    p(16) = 93993D0/625D0
    p(17) = -10659D0/10000D0
    p(18) = -969/25D0
    p(19) = -57D0/10D0
    p(20) = 4D0
    p(21) = 1D0
    exact_roots = (/ (cmplx(-2.1D0+0.2D0*(j-1),0,kind=dp), j=1,deg)/)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 5: reverse Wilkinson deg 10
    deg = 10
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1D0/3628800D0
    p(2) = -11D0/725760D0
    p(3) = 11D0/30240D0
    p(4) = -121D0/24192D0
    p(5) = 7513D0/172800D0
    p(6) = -8591D0/34560D0
    p(7) = 341693D0/362880D0
    p(8) = -84095D0/36288D0
    p(9) = 177133D0/50400D0
    p(10) = -7381D0/2520D0
    p(11) = 1D0
    exact_roots = (/ (cmplx(1.0D0/dble(deg-j+1),0,kind=dp), j=1,deg)/)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 6: reverse Wilkinson deg 15
    deg = 15
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = -1D0/1307674368000D0
    p(2) = 1D0/10897286400D0
    p(3) = -47D0/9340531200D0
    p(4) = 1D0/5987520D0
    p(5) = -26921D0/7185024000D0
    p(6) = 109D0/1814400D0
    p(7) = -324509D0/457228800D0
    p(8) = 4783D0/762048D0
    p(9) = -54576553D0/1306368000D0
    p(10) = 2271089D0/10886400D0
    p(11) = -277382447D0/359251200D0
    p(12) = 2065639D0/997920D0
    p(13) = -35118025721D0/9081072000D0
    p(14) = 13215487D0/2802800D0
    p(15) = -1195757D0/360360D0
    p(16) = 1D0
    exact_roots = (/ (cmplx(1.0D0/dble(deg-j+1),0,kind=dp), j=1,deg)/)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 7: reverse Wilkinson deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1D0/2432902008176640000D0
    p(2) = -1D0/11585247657984000D0
    p(3) = 31D0/3658499260416000D0
    p(4) = -1D0/1935713894400D0
    p(5) = 3931D0/179338199040000D0
    p(6) = -587D0/853991424000D0
    p(7) = 12437081D0/753220435968000D0
    p(8) = -31849D0/102478970880D0
    p(9) = 384794917D0/82771476480000D0
    p(10) = -7321967D0/131383296000D0
    p(11) = 2965638101D0/5518098432000D0
    p(12) = -109542331D0/26276659200D0
    p(13) = 12196364570297D0/470762772480000D0
    p(14) = -573738838201D0/4483454976000D0
    p(15) = 6670985204447D0/13450364928000D0
    p(16) = -6670985204447D0/13450364928000D0
    p(17) = 52460655692911D0/15878903040000D0
    p(18) = -13334148911D0/2520460800D0
    p(19) = 665690574539D0/117327450240D0
    p(20) = -55835135D0/15519504D0
    p(21) = 1D0
    exact_roots = (/ (cmplx(1.0D0/dble(deg-j+1),0,kind=dp), j=1,deg)/)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 8: prescribed roots of varying scale deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 9: prescribed roots of varying scale -3 deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 2765140455576880316286330097421875D0/18014398509481984D0
    p(2) = 20003336218539558834627071739613125D0/36028797018963968D0
    p(3) = 26504589049384252861409184537893125D0/36028797018963968D0
    p(4) = 4809495595975287378276611244229875D0/18014398509481984D0
    p(5) = -1857674437365958001629359052983525D0/4503599627370496D0
    p(6) = -11698953582630728570229643313343213D0/18014398509481984D0
    p(7) = -7974397567058086827152963496557445D0/18014398509481984D0
    p(8) = -776358156363835911942964026680595D0/4503599627370496D0
    p(9) = -641441959755400789084497279447735D0/18014398509481984D0
    p(10) = -27650873494903018971933761124915D0/36028797018963968D0
    p(11) = 57063078564052886214162370894557D0/36028797018963968D0
    p(12) = 7079170685683534011791963440605D0/18014398509481984D0
    p(13) = 144483000592453610567900626155D0/4503599627370496D0
    p(14) = -974289645931023019776572595D0/562949953421312D0
    p(15) = -14061288187707477104464845D0/35184372088832D0
    p(16) = -3248268654710998505043D0/1099511627776D0
    p(17) = 22981827635371262835D0/17179869184D0
    p(18) = -4536701774077635D0/134217728D0
    p(19) = 153546333355D0/524288D0
    p(20) = -987135D0/1024D0
    p(21) = 1D0
    exact_roots = (/ (cmplx(2.0D0**(-10+(j-1))-3.0D0,0,kind=dp), j=1,deg)/)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 10: Chebyshev polynomial deg 20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p(1) = 1D0
    p(2) = 0D0
    p(3) = -200D0
    p(4) = 0D0
    p(5) = 6600D0
    p(6) = 0D0
    p(7) = -84480D0
    p(8) = 0D0
    p(9) = 549120D0
    p(10) = 0D0
    p(11) = -2050048D0
    p(12) = 0D0
    p(13) = 4659200D0
    p(14) = 0D0
    p(15) = -6553600D0
    p(16) = 0D0
    p(17) = 5570560D0
    p(18) = 0D0
    p(19) = -2621440D0
    p(20) = 0D0
    p(21) = 524288D0
    exact_roots = (/ (cmplx(cos((2.0D0*j-1.0D0)*pi/40D0),0,kind=dp), j=1,deg)/)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1))/dble(p(deg+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    ! Poly 11: z^i for i=0,20
    deg = 20
    allocate(exact_roots(deg), p(deg+1), err(deg))
    p = (/ (cmplx(1,0,kind=dp), j=1,deg+1)/)
    exact_roots = (/ (cmplx(cos(2.0D0*j*pi/21.0D0),sin(2.0D0*j*pi/21.0D0),kind=dp), j=1,deg)/)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
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
    allocate(rcoeffs(deg), icoeffs(deg), reigs(deg), ieigs(deg), its(deg), zeros(deg))
    rcoeffs = (/ (dble(p(deg-j+1)), j=1,deg)/)
    icoeffs = (/ (aimag(p(deg-j+1)), j=1,deg)/)
    call zamvw(deg, rcoeffs, icoeffs, reigs, ieigs, its, flag)
    zeros = (/ (cmplx(reigs(j),ieigs(j), kind=dp), j=1,deg)/)
    call sort(zeros, exact_roots, deg)
    err = (/ (abs(zeros(j)-exact_roots(j))/abs(exact_roots(j)), j=1,deg)/)
    write(1, '(ES15.2)') sum(err)/deg
    deallocate(rcoeffs, icoeffs, reigs, ieigs, its, zeros)
    deallocate(exact_roots, p, err)
    close(1)
    call execute_command_line('cp data_files/spec_poly.dat TeX/table_spec_poly.dat')
contains
    !****************************************
    !               init_random_seed        *
    !****************************************
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
    !****************************************
    !               sort                    *
    !****************************************
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