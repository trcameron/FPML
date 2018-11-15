program fixmodx

!   This converts MPFUN calls in file mpmod90.f to the advanced equivalents.  
!   In other words, calls to MPMUL are converted to call MPMULX, etc.
!   This program is used as follows to produce a converted file mpmod90x.f.
!   For example, using gfortran, one would type:

!   gfortran -ffree-form -O2 -o fixmodx fixmodx.f
!   ./fixmodx < mpmod90.f > mpmod90x.f

!   In addition to the changes produced by this program, the line

!   call mpinix (mpnw)

!   (without the ! of course) must be inserted in subroutine mpinit after
!   the line  mpnw = mpwds + 1 (this is one of the first executable lines)
!   in the output file mpmod90x.f.  Also, of course, you should change the
!   value of mpipl in the parameter statement
!   parameter (mpipl = 2000, mpiou = 56, mpiep = 10 - mpipl)
!   (which is one of the first non-comment lines in the file mpmod90x.f)
!   to be a higher value -- otherwise there is no point in using this
!   converted file.

!   Users should also note that calls to MPWRITE are not efficient when the 
!   output precision level is above about 10,000 digits or so.  For very high
!   precision levels such as this, call the advanced routine MPOUTXX from the 
!   MPFUN package for binary-to-decimal conversion, and then use your own
!   WRITE statement for output (MPOUTXX produces a array of character*1 data).
!   For example, the code

!   integer ndp
!   parameter (ndp = 100000)
!   call mpinit (ndp)
!   zz = mppic
!   call mpwrite (6, zz)

!   can be replaced by the followinog:

!   integer i, ndp, nchx
!   parameter (ndp = 100000)
!   character*1 chx(ndp+200)
!   call mpinit (ndp)
!   zz = mppic
!   call mpoutcx (zz%mpr, chx, nchx, mpnwx)
!   write (6, '(78a1)') (chx(i), i = 1, nchx)

!   In the above code, note that the array chx must be allocated with at
!   at least ndp+200 character*1 cells, or memory overwrite errors can occur.

!   David H. Bailey   2010-07-16

character*80 lin, linx

ilog = 0

100 read (5, '(a)', end = 200) lin
linx = lin

do i = 1, 40
  if (lin(i:i+4) .eq. 'call ') then
    j = i + 5
    if (lin(j:j+5) .eq. 'mpang ') linx(j:j+5) = 'mpangx'
    if (lin(j:j+5) .eq. 'mpcdiv') linx(j:j+5) = 'mpcdvx'
    if (lin(j:j+5) .eq. 'mpcmul') linx(j:j+5) = 'mpcmlx'
    if (lin(j:j+5) .eq. 'mpcpwr') linx(j:j+5) = 'mpcpwx'
    if (lin(j:j+5) .eq. 'mpcssh') then
      linx(j:j+5) = 'mpcshx'
      ix = j + 6 + index (lin(j+7:), ' ')
      linx(ix+1:ix+11) = 'mppic%mpr, '
      linx(ix+12:) = lin(ix+1:)
    endif
    if (lin(j:j+5) .eq. 'mpcsqr') linx(j:j+5) = 'mpcsqx'
    if (lin(j:j+5) .eq. 'mpcssn') linx(j:j+5) = 'mpcssx'
    if (lin(j:j+5) .eq. 'mpdiv ') linx(j:j+5) = 'mpdivx'
    if (lin(j:j+5) .eq. 'mpexp ') then
      linx(j:j+5) = 'mpexpx'
      ix = j + 6 + index (lin(j+7:), ' ')
      linx(ix+1:ix+11) = 'mppic%mpr, '
      linx(ix+12:) = lin(ix+1:)
    endif
    if (lin(j:j+5) .eq. 'mplog ') then
      ilog = ilog + 1
      linx(j:j+5) = 'mplogx'
      ix = j + 6 + index (lin(j+7:), ' ')
      if (ilog .le. 2) then
        linx(ix+1:ix+6) = 't1, '
        linx(ix+7:) = lin(ix+1:)
      else
        linx(ix+1:ix+11) = 'mppic%mpr, '
        linx(ix+12:) = lin(ix+1:)
      endif
    endif
    if (lin(j:j+5) .eq. 'mpmul ') linx(j:j+5) = 'mpmulx'
    if (lin(j:j+5) .eq. 'mpnpwr') linx(j:j+5) = 'mpnpwx'
    if (lin(j:j+5) .eq. 'mpnrt ') linx(j:j+5) = 'mpnrtx'
    if (lin(j:j+5) .eq. 'mppi (') linx(j:j+5) = 'mppix('
    if (lin(j:j+5) .eq. 'mpsqrt') linx(j:j+5) = 'mpsqrx'
  endif
enddo

ln = lnblk (linx)
write (6, '(a)') linx(1:ln)
goto 100

200 stop
end

function lnblk (lin)
character*(*) lin

!   This finds the index of the last non-blank character in LIN.

ln = len (lin)

do i = ln, 1, -1
  if (lin(i:i) .ne. ' ') goto 110
enddo

i = 0
110 lnblk = i

return
end
