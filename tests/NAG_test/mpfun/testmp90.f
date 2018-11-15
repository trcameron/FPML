program testmp

!   This is the test program for DHB's multiprecision computation package
!   MPFUN.  It exercises most routines and verifies that they are working
!   properly.  IEEE Fortran-90 version.

!   David H. Bailey    2005-03-22

use mpfunmod
double precision ds, dsum, t1, x1
character*1 z
character*80 za
parameter (mx = 6, nx = 2 ** mx, l = 5, n = 8, n4 = nx + 4, &
  ln = l * n, mt = 10 - 6 * nx, np = 20, nt = 64)
dimension aa(l,n), ac(2*l,n), a(n4), al2(n4), b(n4), ds(nt), &
  f1(8), f2(8), ip(np), is1(np), is2(9), nx1(2), pi(n4), x(2*n4), &
  xx(n4,np), x1(2), y(2*n4), z(600)
data aa / &
   1., 0., 40., 0., 0., 1., 0., 24., 0., 0., -1., 0., 15., 0., 0., &
  -1., 0.,  8., 0., 0., 1., 0.,  5., 0., 0.,  1., 0.,  2., 0., 0., &
   0., 0.,  0., 0., 0., 1., 0.,  1., 0., 0./
data ac / &
   1., 0., 40., 0., 0., 5 * 0.,  1., 0., 24., 0., 0., 5 * 0., &
  -1., 0., 15., 0., 0., 5 * 0., -1., 0.,  8., 0., 0., 5 * 0., &
   1., 0.,  5., 0., 0., 5 * 0.,  1., 0.,  2., 0., 0., 5 * 0., &
   0., 0.,  0., 0., 0., 5 * 0.,  1., 0.,  1., 0., 0., 5 * 0./
data ds / &
          0.d0,          0.d0,  532002022.d0,  532002022.d0, &
  535684389.d0,  466695539.d0,  535684389.d0,  466695539.d0, & ! 8
  598587746.d0,  598587746.d0,  544393923.d0,  544393923.d0, &
  550529242.d0,  550529242.d0,  536143302.d0,  516904795.d0, & !16
  536143302.d0,  516904795.d0,  564046428.d0,  564046432.d0, &
  560037982.d0,  546010881.d0,  485070596.d0,  485070596.d0, & !24
  534435267.d0,  534435267.d0,  535883606.d0,  535883606.d0, &
  501186251.d0,  501186250.d0,  471471731.d0,  471471730.d0, & !32
  531664183.d0,  531664183.d0,  531664196.d0,  531664171.d0, &
  573727184.d0,  512972074.d0,  573727185.d0,  512972073.d0, & !40
  510801511.d0,  520623984.d0,  510801523.d0,  520623975.d0, &
  526887598.d0,  526887563.d0,  470800624.d0,  470800624.d0, & !48
  501444259.d0,  501444259.d0,  515452421.d0,  585480619.d0, &
  515452421.d0,  585480292.d0,  561015607.d0,  603898762.d0, & !56
  560752017.d0,  600941033.d0,  529966523.d0,  517231394.d0, &
         98.d0,      3 * 0.d0/
data is1 / &
  20,   1,  11,  10,  16,   5,  17,  12,  15,  13,   9,   4, &
   2,  18,   6,   3,   7,  19,   8,  14/

!   Initialize.  MCR is set to 5 so that the advanced routines can be tested
!   with reasonably short run times.

mpnw = nx
mpwd = 2005 / 7.224719896d0 + 2.d0
mpier = 0
mpmcr = 5
call mpiniq (mpwd, mpnw)
call mpinix (nx)

!   Test the input/output conversion routines.
!   These tests have been disabled, because the routines are now in mpmod90.f.

kt = 1
! za = '10 ^ - 50 x - 3. 14159 26535 89793 23846 26433 83279 50288'
! read (za, '(80A1)') (z(i), i = 1, 80)
! call mpinpc (z, 80, a, mpnw)
! ls = 9
! if (dsum (ls, a) .ne. ds(kt)) call derr (kt, ls, a)
write (6, 1) kt

1 format ('Completed Test',i3)
2  format ('TESTMP Failed On Test No.',i4/'Result: ',60a1)

kt = 2
! za ='10 ^       -50 x -3.14159265358979323846264338327950288'
! call mpoutc (a, z, nn, mpnw)
! nn = min (nn, 55)

! do j = 1, nn
!   if (z(j) .ne. za(j:j)) then
!     write (6, 2) kt, (z(i), i = 1, nn)
!     goto 110
!   endif
! enddo

110  write (6, 1) kt

!   Compute 3. ^ (-13).

kt = 3
f1(1) = 1.
f1(2) = 0.
f1(3) = 3.
f1(4) = 0.
nn = -13
ls = mpnw + 2
call mpnpwr (f1, nn, a, mpnw)
if (dsum (ls, a) .ne. ds(kt)) call derr (kt, ls, a)
write (6, 1) kt

kt = 4
call mpnpwx (f1, nn, a, mpnw)
if (dsum (ls, a) .ne. ds(kt)) call derr (kt, ls, a)
write (6, 1) kt

!   Compute (4 - i) ^ (-25).

kt = 5
f1(1) = 1.
f1(2) = 0.
f1(3) = 4.
f1(4) = 0.
f2(1) = -1.
f2(2) = 0.
f2(3) = 1.
f2(4) = 0.
call mpmmpc (f1, f2, n4, x, mpnw)
nn = -25
ls = mpnw + 2
call mpcpwr (n4, x, nn, y, mpnw)
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt
kt = 6
if (dsum (ls, y(n4+1)) .ne. ds(kt)) call derr (kt, ls, y(n4+1))
write (6, 1) kt

kt = 7
call mpcpwx (n4, x, nn, y, mpnw)
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt
kt = 8
if (dsum (ls, y(n4+1)) .ne. ds(kt)) call derr (kt, ls, y(n4+1))
write (6, 1) kt

!   Compute Sqrt (Sqrt (10)).

kt = 9
f1(1) = 1.
f1(2) = 0.
f1(3) = 10.
f1(4) = 0.
call mpsqrt (f1, a, mpnw)
call mpsqrt (a, b, mpnw)
ls = mpnw + 2
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 10
call mpsqrx (a, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Cbrt (Sqrt (10)).

kt = 11
call mpcbrt (a, b, mpnw)
ls = mpnw + 2
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 12
call mpcbrx (a, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute the 10th root of 10.

kt = 13
call mpnrt (f1, 10, b, mpnw)
ls = mpnw + 2
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 14
call mpnrtx (f1, 10, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute the complex square root of (2., 1.).

kt = 15
f1(1) = 1.
f1(2) = 0.
f1(3) = 2.
f1(4) = 0.
f2(1) = 1.
f2(2) = 0.
f2(3) = 1.
f2(4) = 0.
call mpmmpc (f1, f2, n4, x, mpnw)
call mpcsqt (n4, x, y, mpnw)
ls = mpnw + 1
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt
kt = 16
if (dsum (ls, y(n4+1)) .ne. ds(kt)) call derr (kt, ls, y(n4+1))
write (6, 1) kt

kt = 17
call mpcsqx (n4, x, y, mpnw)
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt
kt = 18
if (dsum (ls, y(n4+1)) .ne. ds(kt)) call derr (kt, ls, y(n4+1))
write (6, 1) kt

!   Compute Pi.

kt = 19
call mppi (b, mpnw)
ls = mpnw + 2
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt
call mpeq (b, pi, mpnw)

kt = 20
call mppix (b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Log (2).

kt = 21
f1(1) = 1.
f1(2) = 0.
f1(3) = 2.
f1(4) = 0.
call mplog (f1, al2, b, mpnw)
ls = mpnw + 2
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt
call mpeq (b, al2, mpnw)

kt = 22
call mplogx (f1, pi, al2, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Log (10).

kt = 23
f1(3) = 10.
call mplog (f1, al2, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 24
call mplogx (f1, pi, al2, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Log (1/4).

kt = 25
call mpdmc (0.25d0, 0, f1)
call mplog (f1, al2, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 26
call mplogx (f1, pi, al2, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Exp (1).

kt = 27
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
call mpexp (f1, al2, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 28
call mpexpx (f1, pi, al2, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Exp (2).

kt = 29
f1(3) = 2.
call mpexp (f1, al2, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 30
call mpexpx (f1, pi, al2, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Exp (-5).

kt = 31
f1(1) = -1.
f1(3) = 5.
call mpexp (f1, al2, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 32
call mpexpx (f1, pi, al2, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Cos and Sin of Pi/4.

kt = 33
call mpmuld (pi, 0.25d0, 0, a, mpnw)
call mpcssn (a, pi, x, y, mpnw)
ls = mpnw + 1
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 34
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt

kt = 35
call mpcssx (a, pi, x, y, mpnw)
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 36
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt

!   Compute Cos and Sin of 39/64 Pi.

kt = 37
call mpmuld (pi, 0.609375d0, 0, a, mpnw)
call mpcssn (a, pi, x, y, mpnw)
ls = mpnw + 1
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 38
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt

kt = 39
call mpcssx (a, pi, x, y, mpnw)
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 40
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt

!   Compute Cos and Sin of -19/64 Pi.

kt = 41
call mpmuld (pi, -0.296875d0, 0, a, mpnw)
call mpcssn (a, pi, x, y, mpnw)
ls = mpnw + 1
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 42
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt

kt = 43
call mpcssx (a, pi, x, y, mpnw)
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 44
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt

!   Compute inverse Cos and Sin of (1, 1).

kt = 45
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
f2(1) = 1.
f2(2) = 0.
f2(3) = 1.
f2(4) = 0.
call mpang (f1, f2, pi, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 46
call mpangx (f1, f2, pi, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute inverse Cos and Sin of (1, 5).

kt = 47
f2(3) = 5.
call mpang (f1, f2, pi, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 48
call mpangx (f1, f2, pi, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute inverse Cos and Sin of (-1, 3).

kt = 49
f1(1) = -1.
f2(3) = 3.
call mpang (f1, f2, pi, b, mpnw)
ls = mpnw + 1
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

kt = 50
call mpangx (f1, f2, pi, b, mpnw)
if (dsum (ls, b) .ne. ds(kt)) call derr (kt, ls, b)
write (6, 1) kt

!   Compute Cosh and Sinh of 0.5.

kt = 51
call mpdmc (0.5d0, 0, f1)
call mpcssh (f1, al2, x, y, mpnw)
ls = mpnw + 1
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 52
if (dsum (ls, y) .ne. ds(kt)) call derr (kt, ls, y)
write (6, 1) kt

kt = 53
call mpcshx (f1, pi, al2, x, y, mpnw)
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

kt = 54
if (abs (dsum (ls, y) - ds(kt)) .gt. 1) call derr (kt, ls, y)
write (6, 1) kt

!   Compute the root near x = 1.42 + 0.69i of the polynomial
!   x^7 + 2 x^5 + 5 x^4 - 8 x^3 - 15 x^2 + 24 x + 40 = 0.

kt = 55
x1(1) = 1.42d0
nx1(1) = 0
x1(2) = 0.69d0
nx1(2) = 0
call mpcpol (n - 1, l, ac, x1, nx1, n4, x, mpnw)
ls = mpnw + 2
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt
kt = 56
if (dsum (ls, x(n4+1)) .ne. ds(kt)) call derr (kt, ls, x(n4+1))
write (6, 1) kt

kt = 57
call mpcplx (n - 1, l, ac, x1, nx1, n4, x, mpnw)
ls = mpnw + 1
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt
kt = 58
if (dsum (ls, x(n4+1)) .ne. ds(kt)) call derr (kt, ls, x(n4+1))
write (6, 1) kt

!   Compute the real root of the above polynomial near x = -1.34.

kt = 59
t1 = -1.34d0
n1 = 0
call mppol (n - 1, l, aa, t1, n1, x, mpnw)
ls = mpnw + 2
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
call mpeq (x, y, mpnw)
write (6, 1) kt

kt = 60
call mppolx (n - 1, l, aa, t1, n1, x, mpnw)
ls = mpnw + 1
if (dsum (ls, x) .ne. ds(kt)) call derr (kt, ls, x)
write (6, 1) kt

!   Sort a pseudo-randomly generated vector.

kt = 62

do j = 1, np
  call mprand (xx(1,j), mpnw)
enddo

call mpsort (np, n4, xx, ip, mpnw)
if (ichk (np, ip, is1) .ne. 0) call ierr (kt, np, ip)
write (6, 1) kt

stop
end

function dsum (n, a)
double precision dsum, s
dimension a(n)

s = 0.d0

do i = 1, n
  s = s + a(i)
enddo

dsum = s
return
end

subroutine derr (n, l, a)
use mpfunmod
double precision dsum
dimension a(l)
character*1 cx(1000)

write (6, 1) n
1 format ('TESTMP Failed On Test No.',i4)
write (6, 2) (a(i), i = 1, l)
2 format ('Result:'/(6f12.0))
! call mpout (6, a, int (7.225 * (l - 2)), cx, mpnw)
write (6, 3) dsum (l, a)
3 format ('Checksum:', f20.0)
return
end

function ichk (n, ia, ib)
dimension ia(n), ib(n)

is = 0

do i = 1, n
  is = is + abs (ia(i) - ib(i))
enddo

ichk = is
return
end

subroutine ierr (n, l, ia)
dimension ia(l)

write (6, 1) n
1 format ('TESTMP Failed On Test No.',i4)
write (6, 2) (ia(i), i = 1, l)
2 format ('Result:'/(8i9))
return
end
