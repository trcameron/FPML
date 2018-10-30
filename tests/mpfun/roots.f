program troots

!   This program finds polynomial roots.  Double precision, double complex,
!   arbitrary precision real and arbitrary precision complex versions are
!   included.
!   David H Bailey   2004-05-27

use mpmodule
implicit none
integer i, n, ndp, nr
parameter (n = 10, ndp = 500)
real*8 d1, d2, d3, da(0:n), second, tm0, tm1
! real*8 dr, dr0
complex*16 dr, dr0
type (mp_real) t1, t2, t3, a(0:n)
! type (mp_real) r0, r
type (mp_complex) r0, r
data da/1.d0, 1.d0, 0.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, 1.d0, 1.d0/
! data dr0/1.176d0/
data dr0/(0.457d0, 0.889d0)/

call mpinit
call mpsetprec (ndp)
mpoud = ndp

write (6, 1) da
1 format ('Arbitrary precision polynomial root finder.'/'Coefficients:'/ &
  (1p,4d19.10))

call dcroot (n, dr0, da, nr, dr)
write (6, 2) dr
2 format ('Double precision root:'/(1p,2d19.10))

!  Copy coefficients in da to a.

do i = 0, n
  a(i) = da(i)
enddo

r0 = dr0

!  Find root of polynomial.

tm0 = second ()
! call rroot (n, r0, a, nr, r)
call croot (n, r0, a, nr, r)
tm1 = second ()
write (6, 3) tm1 - tm0
3 format ('cpu time =',f8.2)
write (6, 4)
4 format ('Arbitrary precision root:')
call mpwrite (6, r)
stop
end

subroutine drroot (n, dr0, da, nr, dr)

!   Finds a single double-precision real root, beginning at dr0.

implicit none
integer i, j, k, n, nit, nr
parameter (nit = 20)
real*8 da(0:n), dad(0:n-1), eps, eps1
real*8 d1, d2, d3, d4, dr, dr0
parameter (eps = 1.d-15)

!   Compute derivative.

do i = 0, n - 1
  dad(i) = (i + 1) * da(i+1)
enddo

d1 = dr0
eps1 = eps * max (abs (d1), 1.d0)

!   Perform up to nit double-precision Newton iterations.

do j = 1, nit

!   Evaluate polynomial at d1.

  d2 = 0.d0

  do i = n, 0, -1
    d2 = d1 * d2 + da(i)
  enddo

!   Evaluate polynomial derivative at d1.

  d3 = 0.d0

  do i = n - 1, 0, -1
    d3 = d1 * d3 + dad(i)
  enddo

  d4 = d1 - d2 / d3
  if (abs (d1 - d4) < eps1) goto 100
  d1 = d4
enddo

write (6, 1)
1 format ('drroot: failed to find double precision real root.')
nr = 0
dr = 0.d0
goto 110

100 continue

nr = 1
dr = d1

110 continue

return
end

subroutine dcroot (n, dr0, da, nr, dr)

!   Finds a single double-precision complex root, beginning at dr0.

implicit none
integer i, j, k, n, nit, nr
parameter (nit = 20)
real*8 da(0:n), dad(0:n-1), eps, eps1
complex*16 d1, d2, d3, d4, dr, dr0
parameter (eps = 1.d-15)

!   Compute derivative.

do i = 0, n - 1
  dad(i) = (i + 1) * da(i+1)
enddo

d1 = dr0
eps1 = eps * max (abs (d1), 1.d0)

!   Perform up to nit double-precision Newton iterations.

do j = 1, nit

!   Evaluate polynomial at d1.

  d2 = 0.d0

  do i = n, 0, -1
    d2 = d1 * d2 + da(i)
  enddo

!   Evaluate polynomial derivative at d1.

  d3 = 0.d0

  do i = n - 1, 0, -1
    d3 = d1 * d3 + dad(i)
  enddo

  d4 = d1 - d2 / d3
  if (abs (d1 - d4) < eps1) goto 100
  d1 = d4
enddo

write (6, 1)
1 format ('dcroot: failed to find double precision complex root.')
nr = 0
dr = 0.d0
goto 110

100 continue

nr = 1
dr = d1

110 continue

return
end

subroutine rroot (n, r0, a, nr, r)

!   Finds a single arbitrary-precision real root, beginning at r0.

use mpmodule
implicit none
integer i, j, k, n, nit, nr, nwp, nws
parameter (nit = 100)
type (mp_real) a(0:n), ad(0:n-1), eps1
type (mp_real) t1, t2, t3, t4, r, r0

!   Compute derivative.

do i = 0, n - 1
  ad(i) = (i + 1) * a(i+1)
enddo

call mpgetprecwords (nws)
nwp = 3
call mpsetprecwords (nwp)
t1 = r0
t2 = 0.d0
if (r0 /= t2) then
  eps1 = 2.d0 ** (-96) * abs (t1)
else
  eps1 = 2.d0 ** (-96)
endif

!   Perform up to nit double-precision Newton iterations.

do j = 1, nit

!   Evaluate polynomial at t1.

  t2 = 0.d0

  do i = n, 0, -1
    t2 = t1 * t2 + a(i)
  enddo

!   Evaluate polynomial derivative at t1.

  t3 = 0.d0

  do i = n - 1, 0, -1
    t3 = t1 * t3 + ad(i)
  enddo

  t4 = t1 - t2 / t3
  if (nwp == 3) then
    if (abs (t1 - t4) < eps1) then
      nwp = min (2 * nwp - 1, nws)
      call mpsetprecwords (nwp)
    endif
  elseif (nwp < nws) then
    nwp = min (2 * nwp - 1, nws)
    call mpsetprecwords (nwp)
  else
    t1 = t4
    goto 100
  endif
  t1 = t4
enddo

write (6, 1)
1 format ('rroot: failed to find arbitrary precision real root.')
nr = 0
r = 0.d0
goto 110

100 continue

nr = 1
r = t1

110 continue

return
end

subroutine croot (n, r0, a, nr, r)

!   Finds a single arbitrary-precision complex root, beginning at r0.

use mpmodule
implicit none
integer i, j, k, n, nit, nr, nwp, nws
parameter (nit = 100)
type (mp_real) a(0:n), ad(0:n-1), eps1
type(mp_complex) t1, t2, t3, t4, r, r0

!   Compute derivative.

do i = 0, n - 1
  ad(i) = (i + 1) * a(i+1)
enddo

call mpgetprecwords (nws)
nwp = 3
call mpsetprecwords (nwp)
t1 = r0
t2 = (0.d0, 0.d0)
if (r0 /= t2) then
  eps1 = 2.d0 ** (-96) * abs (t1)
else
  eps1 = 2.d0 ** (-96)
endif

!   Perform up to nit double-precision Newton iterations.

do j = 1, nit

!   Evaluate polynomial at t1.

  t2 = (0.d0, 0.d0)

  do i = n, 0, -1
    t2 = t1 * t2 + a(i)
  enddo

!   Evaluate polynomial derivative at t1.

  t3 = (0.d0, 0.d0)

  do i = n - 1, 0, -1
    t3 = t1 * t3 + ad(i)
  enddo

  t4 = t1 - t2 / t3
  if (nwp == 3) then
    if (abs (t1 - t4) < eps1) then
      nwp = min (2 * nwp - 1, nws)
      call mpsetprecwords (nwp)
    endif
  elseif (nwp < nws) then
    nwp = min (2 * nwp - 1, nws)
    call mpsetprecwords (nwp)
  else
    t1 = t4
    goto 100
  endif
  t1 = t4
enddo

write (6, 1)
1 format ('croot: failed to find multiple precision complex root.')
nr = 0
r = (0.d0, 0.d0)
goto 110

100 continue

nr = 1
r = t1

110 continue

return
end


