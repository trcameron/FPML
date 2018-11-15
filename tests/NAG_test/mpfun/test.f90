program test
use mpfunmod
use mpmodule
implicit none
integer i, ndp, nchx
parameter (ndp = 100000)
type (mp_real) t1, t2
character*1 chx(ndp+200)

call mpinit (ndp)
mpoud = ndp
t1 = mppic / 10.d0
t2 = cos (t1)**2 + sin(t1)**2
write (6, *) 'cos(t1)^2 + sin(t1)^2 ='
call mpwrite (6, t2)
write (6, *) 'pi to ', ndp, ' digits ='

 call mpwrite (6, mppic)

!   This code does the equivalent of the previous line with mpwrite, but is much
!   faster in execution (note however this requires the character*1 line above):

call mpoutcx (mppic%mpr, chx, nchx, mpnwx)
write (6, '(78a1)') (chx(i), i = 1, nchx)

stop
end 
