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

subroutine DCORR(JREFX,AREF,ICSPCK,INTSYM,INDX,DMO)

use mrci_global, only: ENP, IPRINT, IRC, LN, Lu_27, NREF
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: JREFX(*), ICSPCK(*), INTSYM(*), INDX(*)
real(kind=wp) :: AREF(*), DMO(*)
integer(kind=iwp) :: I, IAD27, II1, IJ, IK, INDA, IOC
real(kind=wp) :: FAC, TSUM
integer(kind=iwp), external :: ICUNP
!Statement function
integer(kind=iwp) :: JO, L
JO(L) = ICUNP(ICSPCK,L)

! CORRECTION TO DENSITY MATRIX IN ACPF CASE.
if (IPRINT >= 7) write(u6,*) ' ENP IN DENS =',ENP
FAC = One-(One/ENP)
IAD27 = 0
call dDAFILE(Lu_27,2,AREF,NREF,IAD27)
IK = 0
do INDA=1,IRC(1)
  II1 = (INDA-1)*LN
  if (JREFX(INDA) /= 0) then
    IK = IK+1
    TSUM = AREF(IK)*AREF(IK)*FAC
    IJ = 0
    do I=1,LN
      IOC = (1+JO(II1+I))/2
      IJ = IJ+I
      DMO(IJ) = DMO(IJ)+IOC*TSUM
    end do
  end if
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(INTSYM)
  call Unused_integer_array(INDX)
end if

end subroutine DCORR