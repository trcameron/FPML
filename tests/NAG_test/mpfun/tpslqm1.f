program tpslqm1

!   Test program for PSLQM1 -- a one-level implementation of the multi-pair 
!   PSLQ algorithm for finding integer relations.  Fortran-90 OMP ARPREC version.

!   David H. Bailey     2004-11-23

!   These parameters are set in the parameter statement below.  Some other
!   parameters are set in other routines.

!     idb = Debug level (0 - 3).
!     itr = Number of trials when KQ = 1, 2 or 3.
!     gam = Gamma, the PSLQ algorithm parameter.  Must be set to at least
!           Sqrt (4/3) = 1.154700538d0.
!     n   = Integer relation vector length.
!     kq  = 0 for the algebraic case [1, al, al^2, ... al^(n-1)], where
!           al = 3^(1/kr) - 2^(1/ks).  Set n = kr * ks + 1 and itr = 1.
!         = 1 for testing algebraic relations of a number read from a file.
!         = 2 for testing additive relations of numbers read from a file.
!         = 3 for testing multiplicative relations of numbers read from a file.
!         = 4 for custom input.
!     kr  = Degree of root of 3 when kq = 0.
!     ks  = Degree of root of 2 when kq = 0.
!     ndp = No. digits of pecision.  Must not exceed mpipl in mpmod90.f.
!           It is recommended that this be set to at least 30 digits more
!           than the expected detection level of the relation.
!     nep = Log10 of epsilon for detections.  Normally set to 10 - ndp or so.
!           Must not exceed the accuracy of input data.
!     nsq = Size of duplicate table in iterdp and itermpw.  Typically nsq = 8.
!     rb  = Log10 of max size (Euclidean norm) of acceptable relation.

use mpmodule
implicit none
double precision d1, d2, drand, gam, rb, rm, second, tm, tm0, tm1
integer i, idb, iq, itr, j, k, kps, kcs, kq, kr, ks, kss, n, ndp, nep, nsq
parameter (gam = 1.154700538d0, idb = 2, itr = 1, kq = 0, kr = 5, ks = 5, &
  n = kr * ks + 1, nsq = 8, ndp = 180, nep = -170, rb = 100.d0)
integer lnm(n)
character*1 chr1(100)
character*40 nam(n), namx
double precision dr(n), r1(n)
type (mp_real) al, t1, t2, r(n), x(n)
external second

!   This is scratch space for pslqm1.  See below.

integer is0(3*n)
double precision s1(n)
type (mp_real) s2(3*n*n+nsq*n+2*n)
save

call mpinit (ndp)
mpoud = 100
mpeps = mpreal (10) ** nep
write (6, 1) itr, gam, n, nsq, kq, kr, ks, rb, ndp, nep
1 format ('PSLQM1 Test Program'/'No. trials = ',i3,3x,'Gamma =',f12.9, &
  3x,'n =',i4,3x,'nsq =',i3,3x,'kq =',i2/'kr =',i2,3x,'ks =',i2,3x, &
  'rb =',1p,d12.4,3x/'Precision level ndp =',i6,' digits',3x, &
  'Epsilon level nep = ',i6)
tm = 0.d0
kps = 0
kcs = 0
if (kq .eq. 1 .or. kq .eq. 2 .or. kq .eq. 3) then
  open (11, file = 'pslq.inp')
  rewind 11
endif

do k = 1, itr
  write (6, 2) k
2 format (/'Start trial', i3)

  if (kq .eq. 0) then

!   This code generates al = 3^(1/kr) - 2^(1/ks).  al is algebraic of degree
!   kr * ks.  Set n = kr * ks + 1 to recover the polynomial satisfied by al.

    al = mpnrtf (mpreal (3.d0), kr) - mpnrtf (mpreal (2.d0), ks)
  elseif (kq .eq. 1) then

!   Read an algebraic constant from a file.

    call mpread (11, al)
  elseif (kq .eq. 2) then

!   Read constants from a file for additive test.

    do i = 1, n
      call mpread (11, al)
      x(i) = al
      write (namx, '(''con'',i3.3)') i
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  elseif (kq .eq. 3) then

!   Read constants from a file for multiplicative test.

    do i = 1, n
      call mpread (11, al)
      x(i) = log (al)
      write (namx, '(''log(con'',i3.3,'')'')') i
      nam(i) = namx(1:11)
      lnm(i) = 11
    enddo
  elseif (kq .eq. 4) then

!   Produce X vector by a custom scheme.

  endif

!   If kq is 0 or 1, generate x = [1, al, al^2, ..., al^(n-1)].

  if (kq .eq. 0 .or. kq .eq. 1) then
    x(1) = 1.d0
    nam(1) = '1'
    lnm(1) = 1

    do i = 2, n
      x(i) = al * x(i-1)
      write (namx, '(''al^'',i3)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  endif

!   Perform relation search.

  tm0 = second ()
  call pslqm1 (idb, gam, n, nsq, rb, x, is0, s1, s2, iq, r)
  tm1 = second ()

!   Output relation, if one was found.

  if (iq .ne. 0) then
    tm = tm + (tm1 - tm0)
    write (6, 4)
4   format ('Recovered relation: 0 =')

    do i = 1, n
      if (r(i) .ne. 0.d0) then
        call mpfform (r(i), 24, 0, chr1)
        write (6, '(''+'',24a1,'' * '',a)') (chr1(j), j = 1, 24), &
          nam(i)(1:lnm(i))
      endif
    enddo

!   Check if original relation was recovered.

    if (kq .eq. 1) then
      d1 = 0.d0
      d2 = 0.d0

      do i = 1, n
        d1 = d1 + abs (dble (r(i)) - r1(i))
        d2 = d2 + abs (dble (r(i)) + r1(i))
      enddo

      if (d1 .eq. 0.d0 .or. d2 .eq. 0.d0) then
        kcs = kcs + 1
        write (6, 5)
5       format ('Original relation recovered.')
      else
        kps = kps + 1
      endif
    endif
  endif
  write (6, 6) tm1 - tm0
6 format ('CPU Time =',f12.4)
enddo

if (itr .gt. 1) then
  write (6, 7) kps, kcs
7 format (/'No. partial successes =',i5/'No. complete successes =',i5)
  kss = kps + kcs
  if (kss .ne. 0) then
    tm = tm / kss
    write (6, 8) tm
8   format ('Ave. CPU time of PS or CS runs =',f10.4)
  endif
endif
stop
end

!------------------------------

!   The following code performs the one-level, multi-pair PSLQ algorithm.
!   David H. Bailey     2004-06-10

subroutine pslqm1 (idb, gam, n, nsq, rb, x, is0, s1, s2, iq, r)

!   This routine allocates the scratch space arrays is0, s1 and s2.  Handling
!   scratch space in this manner is not really necessary in Fortran-90, but 
!   this design facilitates straightforward translation to Fortran-77.

!   Arguments are as follows:
!     idb = Debug flag (0-3).
!     gam = PSLQ gamma parameter, normally set to sqrt(4/3).
!     n   = Dimension of input vector.
!     nsq = Size of tables used in itermpw.
!     rb  = Log10 of max size (Euclidean norm) of acceptable relation (type DP).
!     x   = Input real vector (type mp_real).
!     is0 = Scratch array (type integer -- see size below).
!     s1  = Scratch array (type double precision -- see size below).
!     s2  = Scratch array (type mp_real -- see size below).
!     iq  = Output flag: 0 (unsuccessful) or 1 (successful).
!     r   = Output integer coefficient vector (type mp_real).

use mpmodule
implicit none
integer idb, iq, n, nsq
double precision gam, rb
integer is0(3*n)
double precision s1(n)
type (mp_real) s2(3*n*n+nsq*n+2*n), x(n), r(n)

call pslqm1x (idb, gam, n, nsq, rb, x, is0(1), is0(n+1), is0(2*n+1), &
  s1(1), s2(1), s2(n*n+1), s2(2*n*n+1), s2(3*n*n+1), &
  s2(3*n*n+nsq*n+1), s2(3*n*n+nsq*n+n+1), iq, r)
return
end

subroutine pslqm1x (idb, gam, n, nsq, rb, x, ip, ir, is, dx, b, h, t, &
  syq, s, y, iq, r)

!   The following parameters are set in this routine:
!     ipi = Iteration print interval when idb >= 2.
!     ipm = Iteration check interval for MP iterations.
!     itm = Maximum iteration count.  Run is aborted when this is exceeded.

use mpmodule
implicit none
integer i, idb, imq, ipi, ipm, iq, it, itm, izm, j, j1, n, n1, n2, n3, n4, nsq
parameter (ipi = 25, ipm = 100, itm = 100000)
double precision d1, d2, d3, d4, dplog10, gam, rb, second, tm0, tm1, times(2)
integer ip(n), ir(n), is(n)
double precision dx(n)
type (mp_real) b(n,n), h(n,n), s(n), syq(n,nsq), t(n,n), r(n), x(n), &
  y(n), bound, rn, t1, t2, t3, t4
external bound, dplog10, second

!   Initialize.

if (idb .ge. 2) write (6, 1) n
1 format ('PSLQM1 integer relation detection: n =',i5)
iq = 0
it = 0
imq = 0
rn = 0.d0

do i = 1, 2
  times(i) = 0.d0
enddo

if (idb .ge. 2) write (6, 2) it
2 format ('Iteration',i7,3x,'MP initialization')
tm0 = second ()
call initmp (idb, n, nsq, ip, dx, b, h, s, syq, x, y)
tm1 = second ()
times(1) = tm1 - tm0

!   Perform MP iterations.

if (idb .ge. 2) write (6, 3) it
3 format ('Iteration',i7,3x,'Start MP iterations')

100 continue

it = it + 1
if (idb .eq. 3 .or. idb .ge. 2 .and. mod (it, ipi) .eq. 0) write (6, 4) it
4 format ('Iteration',i7)
tm0 = second ()
call itermp (idb, gam, it, n, nsq, ip, ir, is, dx, b, h, s, syq, t, y, &
  imq, izm)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

if (mod (it, ipm) .eq. 0) then

!   Find min and max absolute value of y vector.

  t1 = 1.d300
  t2 = 0.d0

  do i = 1, n
    t1 = min (t1, abs (y(i)))
    t2 = max (t2, abs (y(i)))
  enddo

  if (idb .ge. 2) then
    call decmd (t1, d1, n1)
    call decmd (t2, d2, n2)
    write (6, 5) it, d1, n1, d2, n2
5   format ('Iteration',i7,3x,'Min, max of y =',f11.6,'D',i5, &
      f11.6,'D',i5)
  endif

!   Compute norm bound.

  t1 = bound (n, h)
  rn = max (rn, t1)
  if (idb .ge. 2) then
    call decmd (t1, d1, n1)
    call decmd (rn, d2, n2)
    write (6, 6) it, d1, n1, d2, n2
6   format ('Iteration',i7,3x,'Norm bound =',f11.6,'D',i5,4x,'Max. bound =', &
      f11.6,'D',i5)
  endif
  if (it .gt. itm) then
    if (idb .ge. 1) write (6, 7) itm
7   format ('Iteration limit exceeded',i7)
    goto 120
  endif
  if (dplog10 (rn) .gt. rb) then
    if (idb .ge. 1) write (6, 8) rb
8   format ('Norm bound limit exceeded.',1p,d15.6)
    goto 120
  endif
endif

if (izm .eq. 0) then
  goto 100
elseif (izm .eq. 1) then
  goto 110
elseif (izm .eq. 2) then
  goto 120
endif

110 continue

!   A relation has been detected.

tm0 = second ()
t1 = 1.d300
t2 = 0.d0
t3 = 0.d0

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
  t3 = t3 + r(i) ** 2
enddo

t3 = sqrt (t3)
call decmd (t3, d3, n3)
d3 = d3 * 10.d0 ** n3

!   Output the final norm bound and other info.

if (idb .ge. 1) then
  t4 = bound (n, h)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  call decmd (t3, d3, n3)
  call decmd (t4, d4, n4)
  write (6, 9) it, d1, n1, d2, n2, d4, n4
9 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',0p,f11.6,'D',i5,f11.6,'D',i5/'Max. bound =',f11.6,'D',i5)
  write (6, 10) j1, d3, n3, d1, n1
10 format ('Index of relation =',i4,3x,'Norm =',f11.6,'D',i5,3x, &
  'Residual =',f11.6,'D',i5)
endif

if (dplog10 (t3) .le. rb) then
  iq = 1
else
  if (idb .ge. 2) write (6, 11)
11 format ('Relation is too large.')
endif

120 continue

if (idb .ge. 2) write (6, 12) times
12 format ('CPU times:'/(5f12.2))

return
end

!------------------------------

!   First-level subroutines.

subroutine initmp (idb, n, nsq, ix, dx, b, h, s, syq, x, y)

!   This initializes MP arrays at the beginning.

use mpmodule
implicit none
integer i, idb, j, n, nsq
double precision d1
integer ix(n)
double precision dx(n)
type (mp_real) b(n,n), h(n,n), s(n), syq(n,nsq), x(n), y(n), t1

if (idb .ge. 3) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, ix, dx, x)
endif

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = 0.d0
  enddo

  b(j,j) = 1.d0
enddo

t1 = 0.d0

!   Compute the x vector, the square root of the partial sum of squares of x,
!   and the y vector, which is the normalized x vector.

do i = n, 1, -1
  t1 = t1 + x(i) ** 2
  s(i) = sqrt (t1)
enddo

t1 = 1.d0 / s(1)

do i = 1, n
  y(i) = t1 * x(i)
  s(i) = t1 * s(i)
enddo

!   Compute the initial h matrix.

!$omp parallel do private (i, j, t1) copyin (mpnw)
do j = 1, n - 1
  do i = 1, j - 1
    h(i,j) = 0.d0
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
    h(i,j) = - y(i) * t1
  enddo
enddo
!$omp end parallel do

!   Zero the syq array.

do j = 1, nsq
  do i = 1, n
    syq(i,j) = 0.d0
  enddo
enddo

if (idb .ge. 3) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine itermp (idb, gam, it, n, nsq, ip, ir, is, dx, b, h, q, syq, t, &
  y, imq, izm)

!   This performs one iteration of the PSLQM algorithm using MP arithmetic.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izm, j, j1, j2, k, mpr, mq, n, &
  n1, nsq, ntl
parameter (ntl = 72)
double precision d1, d2, gam
integer ip(n), ir(n), is(n)
double precision dx(n)
type (mp_real) b(n,n), h(n,n), q(n), syq(n,nsq), t(n,n), y(n)
type (mp_real) t1, t2, t3, t4, teps

teps = 2.d0 ** ntl * mpeps
izm = 0
mpr = nint (0.4d0 * n)

!   Compute q vector = {gam^i * |h(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  q(i) = gam ** i * abs (h(i,i))
enddo

call qsortmp (n - 1, q, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest q(i).

do i = 1, n
  is(i) = 0
enddo

if (imq .eq. 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii .eq. 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of y, and rows of b and h.

!$omp parallel do private (i, j, im, im1, t1) copyin (mpnw)
do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = y(im)
  y(im) = y(im1)
  y(im1) = t1

  do i = 1, n
    t1 = b(im,i)
    b(im,i) = b(im1,i)
    b(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = h(im,i)
    h(im,i) = h(im1,i)
    h(im1,i) = t1
  enddo
enddo
!$omp end parallel do

!   Eliminate the "corners" produced by the above permutation in h.

!$omp parallel do private (i, j, im, im1, t1, t2, t3, t4) copyin (mpnw)
do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im .le. n - 2) then
    t1 = h(im,im)
    t2 = h(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = h(i,im)
      t4 = h(i,im1)
      h(i,im) = t1 * t3 + t2 * t4
      h(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo
!$omp end parallel do

!   Perform reduction on h, using the diagonal scheme.  Multipliers are
!   saved in the t array.

do i = 2, n
!$omp parallel do private (ij, j, k) copyin (mpnw)
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      h(ij,j) = h(ij,j) - t(ij,k) * h(k,j)
    enddo

    t(ij,j) = anint (h(ij,j) / h(j,j))
    h(ij,j) = h(ij,j) - t(ij,j) * h(j,j)
  enddo
!$omp end parallel do
enddo

!   Update y, using the t array.  Find min absolute value of y.

t1 = abs (y(n))
j1 = n

do j = 1, n - 1
  do i = j + 1, n
    y(j) = y(j) + t(i,j) * y(i)
  enddo

  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
enddo

!   Update b, using the t array.

!$omp parallel do private (i, j, k) copyin (mpnw)
do k = 1, n
  do j = 1, n - 1
    do i = j + 1, n
      b(j,k) = b(j,k) + t(i,j) * b(i,k)
    enddo
  enddo
enddo
!$omp end parallel do

!  Find the largest entry of b in the same row as the smallest y.

t2 = 0.d0

do i = 1, n
  t2 = max (t2, abs (b(j1,i)))
enddo

if (t1 .le. t2 * teps) then
  if (idb .ge. 2) then
    call decmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'itermp: Small value in y =',f11.6,'D',i5)
  endif
  if (t1 .le. t2 * mpeps) then
    izm = 1
  else
    if (idb .ge. 1) write (6, 2) it
2   format ('Iteration',i7,3x,'itermp: Precision exhausted.')
    izm = 2
  endif
endif

!   Compare the y vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = 0.d0

  do i = 1, n
    t1 = max (t1, abs (y(i) - syq(i,j)))
  enddo

  if (t1 .le. t2 * teps) then
    if (idb .ge. 2) write (6, 3) it, j
 3  format ('Iteration',i7,3x,'itermp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector y in the table syq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  syq(i,k) = y(i)
enddo

if (idb .ge. 3) then
  write (6, 4)
4 format ('itermp: Updated y:')
!  call matoutmd (1, n, ip, dx, y)
  call matoutmp (1, n, y)
  write (6, 5)
5 format ('itermp: Updated b matrix:')
  call matoutmd (n, n, ip, dx, b)
  write (6, 6)
6 format ('itermp: Updated h matrix:')
  call matoutmd (n, n - 1, ip, dx, h)
endif

return
end

!------------------------------

!   Second- and third-level subroutines.

function bound (n, h)

!   This computes the norm bound using MP arithmetic.

use mpmodule
implicit none
integer i, n
type (mp_real) bound,  h(n,n), t1, t2

t1 = 0.d0

do i = 1, n - 1
  t1 = max (t1, abs (h(i,i)))
enddo

bound = 1.d0 / t1

return
end

function dplog10 (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
double precision da, dplog10, t1
type (mp_real) a

call mpmdc (a%mpr, da, ia)
if (da .le. 0.d0) then
  dplog10 = -1d300
else
  dplog10 = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

subroutine decmd (a, b, ib)

!   For input MP value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (mp_real) a

call mpmdc (a%mpr, da, ia)
if (da .ne. 0.d0) then
  t1 = xlt * ia + log10 (abs (da))
  ib = t1
  if (t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end

function drand ()

!   This routine returns a pseudorandom DP floating number nearly uniformly
!   distributed between 0 and 1 by means of a linear congruential scheme.
!   2^28 pseudorandom numbers with 30 bits each are returned before repeating.

double precision drand, f7, r30, sd, t1, t2, t30
parameter (f7 = 78125.d0, r30 = 0.5d0 ** 30, t30 = 2.d0 ** 30)
save sd
data sd/314159265.d0/

t1 = f7 * sd
t2 = aint (r30 * t1)
sd = t1 - t30 * t2
drand = r30 * sd

return
end

subroutine matoutmd (n1, n2, ix, dx, a)

!   This outputs the MP matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
double precision dx(n2)
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'D',i5))
enddo

return
end

subroutine matoutmp (n1, n2, a)

!   This outputs the MP matrix a.  It may be used in place of calls to matoutmd
!   in the code above if greater accuracy is desired in debug output.

use mpmodule
implicit none
integer i, j, n1, n2
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpwrite (6, a(i,j))
  enddo
enddo

return
end

subroutine qsortmp (n, a, ip)

!   This routine sorts the entries of the N-long MP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.

use mpmodule
implicit none
integer i, iq, it, j, jq, jz, k, l, n
integer ip(n), ik(50), jk(50)
type (mp_real) a(n), s0, s1, s2

do i = 1, n
  ip(i) = i
enddo

if (n .eq. 1) return

k = 1
ik(1) = 1
jk(1) = n

130 i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 .lt. a(ip(l))) goto 160
enddo

i = j
goto 190

160 i = l

do l = j, i, -1
  if (s0 .gt. a(ip(l))) goto 180
enddo

j = i
goto 190

180 j = l
if (i .ge. j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue
if (s0 .ge. a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 k = k - 1
jz = 0
if (j .eq. iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 i = i + 1
if (i .eq. jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz .eq. 0) goto 220
if (j - iq .ge. jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 if (k .gt. 0) goto 130

return
end
