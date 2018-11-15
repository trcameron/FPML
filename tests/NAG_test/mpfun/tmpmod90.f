program tmod

!  This is the test progrm for the F-90 based MP translation modules.

!  David H. Bailey   2004-07-12

use mpmodule
implicit type (mp_real) (a-h, o-z)
type (mp_integer) ia, ib, ic
type (mp_complex) c, d, e
parameter (n = 25)
dimension a(n), b(n)

call mpinit (1000)

!   Character-to-MP assignment, generic MPREAL function, pi and e.

x = '1.234567890 1234567890 1234567890 D-100'
ee = exp (mpreal (1.d0))
call mpwrite (6, x, mppic, ee)
s = 0.d0

!   Loop with subscripted MP variables.

do i = 1, n
  a(i) = 2 * i + 1
  b(i) = 2.d0 * a(i) * (a(i) + 1.d0)
  s = s + b(i) ** 2
enddo
call mpwrite (6, s)

!   Expression with mixed MPI and MPR entities.

ia = s
ib = 262144
s = (s + 327.25d0) * mod (ia, 4 * ib)
call mpwrite (6, s)

!   A complex square root reference.

e = sqrt (mpcmpl (2.d0 * s, s))
call mpwrite (6, e)

!   External and intrinsic MP function references in expressions.

s = dot (n, a, b)
t = 2.d0 * sqrt (s) ** 2
call mpwrite (6, s, t)

s = s / 1048576.d0
t = s + 2.d0 * log (s)
x = 3 + nint (t) * 5
call mpwrite (6, s, t, x)

!   A deeply nested expression with function references.

x = (s + (2 * (s - 5) + 3 * (t - 5))) * exp (cos (log (s)))
call mpwrite (6, x)

!   A special MP subroutine call (computes both cos and sin of S).

call mpcssnf (s, x, y)
t = 1.d0 - (x ** 2 + y ** 2)

!   IF-THEN-ELSE construct involving MP variables.

if (abs (t) .lt. mpeps) then
  call mpwrite (6, t)
else
  call mpwrite (6, mpeps)
endif

stop
end

function dot (n, a, b)

!   MP function subprogram.

use mpmodule
type (mp_real) a(n), b(n), dot, s

s = 0.d0

do i = 1, n
  s = s + a(i) * b(i)
enddo

dot = s
return
end
