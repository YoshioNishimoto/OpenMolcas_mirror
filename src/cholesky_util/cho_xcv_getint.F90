!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_XCV_GetInt(irc,ListCD,l_ListCD,ListSP,l_ListSP,NVT,l_nVT,xInt,l_Int)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: calculate integrals (CD|J) where CD belongs to the shell
!          pairs in ListCD and J runs over all Cholesky vectors.
!
! NOTE: this will only work if nQual and iQuAB are properly set to
! correspond to parent diagonals. Also, mySP should be be set to a
! trivial array [i.e. mySP(i)=i].

use Cholesky, only: iOff_col, nDim_Batch, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_ListCD, ListCD(l_ListCD), l_ListSP, ListSP(l_ListSP), l_NVT, NVT(l_NVT), l_Int
real(kind=wp), intent(out) :: xInt(l_Int)
integer(kind=iwp) :: iCD, iSP, iSym, n

! Set return code
irc = 0

! Offsets to symmetry blocks of the integrals
n = 0
do iSym=1,nSym
  iOff_Col(iSym) = n
  n = n+nDim_Batch(iSym)*NVT(iSym)
end do

! Check allocation of xInt
if (n > l_Int) then
  irc = 1
  return
end if

! Calculate integrals
xInt(1:n) = Zero
do iSP=1,l_ListSP
  do iCD=1,l_ListCD
    call Cho_MCA_CalcInt_4(xInt,n,ListCD(iCD),ListSP(iSP))
  end do
end do

end subroutine Cho_XCV_GetInt
