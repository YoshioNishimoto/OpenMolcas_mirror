!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!********************************************************
!** Public-domain library routines used by casvb only. **
!********************************************************
!**********************
!** EISPACK ROUTINES **
!**********************

subroutine tql2(nm,n,d,e,z,ierr)
! this subroutine is a translation of the algol procedure tql2,
! num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
! wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
! this subroutine finds the eigenvalues and eigenvectors
! of a symmetric tridiagonal matrix by the ql method.
! the eigenvectors of a full symmetric matrix can also
! be found if  tred2  has been used to reduce this
! full matrix to tridiagonal form.
!
! on input
!
!    nm must be set to the row dimension of two-dimensional
!      array parameters as declared in the calling program
!      dimension statement.
!
!    n is the order of the matrix.
!
!    d contains the diagonal elements of the input matrix.
!
!    e contains the subdiagonal elements of the input matrix
!      in its last n-1 positions.  e(1) is arbitrary.
!
!    z contains the transformation matrix produced in the
!      reduction by  tred2, if performed.  if the eigenvectors
!      of the tridiagonal matrix are desired, z must contain
!      the identity matrix.
!
!  on output
!
!    d contains the eigenvalues in ascending order.  if an
!      error exit is made, the eigenvalues are correct but
!      unordered for indices 1,2,...,ierr-1.
!
!    e has been destroyed.
!
!    z contains orthonormal eigenvectors of the symmetric
!      tridiagonal (or full) matrix.  if an error exit is made,
!      z contains the eigenvectors associated with the stored
!      eigenvalues.
!
!    ierr is set to
!      zero       for normal return,
!      j          if the j-th eigenvalue has not been
!                 determined after 30 iterations.
!
! calls pythag for  sqrt(a*a + b*b) .
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nm, n
real(kind=wp), intent(inout) :: d(n), e(n), z(nm,n)
integer(kind=iwp), intent(out) :: ierr
integer(kind=iwp) :: i, ii, j, k, l, l1, l2, m, mml
real(kind=wp) :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2
real(kind=wp), external :: pythag

ierr = 0
c3 = Zero ! dummy initialize
s2 = Zero ! dummy initialize
if (n == 1) return

e(1:n-1) = e(2:)

f = Zero
tst1 = Zero
e(n) = Zero

do l=1,n
  j = 0
  h = abs(d(l))+abs(e(l))
  if (tst1 < h) tst1 = h
  ! .......... look for small sub-diagonal element ..........
  do m=l,n
    tst2 = tst1+abs(e(m))
    if (tst2 == tst1) exit
    ! .......... e(n) is always zero, so there is no exit through the bottom of the loop ..........
  end do

  if (m /= l) then
    do
      if (j == 30) then
        ! .......... set error -- no convergence to an eigenvalue after 30 iterations ..........
        ierr = l
        return
      end if
      j = j+1
      ! .......... form shift ..........
      l1 = l+1
      l2 = l1+1
      g = d(l)
      p = (d(l1)-g)/(Two*e(l))
      r = pythag(p,One)
      d(l) = e(l)/(p+sign(r,p))
      d(l1) = e(l)*(p+sign(r,p))
      dl1 = d(l1)
      h = g-d(l)

      d(l2:) = d(l2:)-h

      f = f+h
      ! .......... ql transformation ..........
      p = d(m)
      c = One
      c2 = c
      el1 = e(l1)
      s = Zero
      mml = m-l
      ! .......... for i=m-1 step -1 until l do -- ..........
      do ii=1,mml
        c3 = c2
        c2 = c
        s2 = s
        i = m-ii
        g = c*e(i)
        h = c*p
        r = pythag(p,e(i))
        e(i+1) = s*r
        s = e(i)/r
        c = p/r
        p = c*d(i)-s*g
        d(i+1) = h+s*(c*g+s*d(i))
        ! .......... form vector ..........
        do k=1,n
          h = z(k,i+1)
          z(k,i+1) = s*z(k,i)+c*h
          z(k,i) = c*z(k,i)-s*h
        end do

      end do

      p = -s*s2*c3*el1*e(l)/dl1
      e(l) = s*p
      d(l) = c*p
      tst2 = tst1+abs(e(l))
      if (tst2 <= tst1) exit
    end do
  end if
  d(l) = d(l)+f
end do
! .......... order eigenvalues and eigenvectors ..........
do ii=2,n
  i = ii-1
  k = i
  p = d(i)

  do j=ii,n
    if (d(j) >= p) cycle
    k = j
    p = d(j)
  end do

  if (k == i) cycle
  d(k) = d(i)
  d(i) = p

  do j=1,n
    p = z(j,i)
    z(j,i) = z(j,k)
    z(j,k) = p
  end do

end do

return

end subroutine tql2
