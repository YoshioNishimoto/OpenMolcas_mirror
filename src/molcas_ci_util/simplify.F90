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

subroutine SIMPLIFY(FRAC)

implicit real*8(A-H,O-Z)
integer A, B, FRAC
dimension FRAC(2)

if (FRAC(1) == 0) return
! Find GCD of numerator and denominator
A = FRAC(1)
B = FRAC(2)
60 continue
if (B /= 0) then
  ITEMP = B
  B = mod(A,B)
  A = ITEMP
  goto 60
end if
FRAC(1) = FRAC(1)/A
FRAC(2) = FRAC(2)/A

end subroutine SIMPLIFY