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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

integer function RecNo(itype,iRoot)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute the Record number of a vector                            *
!                                                                      *
!     calling arguments:                                               *
!     itype   : integer                                                *
!               vector type: 1 = H_diag                                *
!                            2 = CI_vec                                *
!                            3 = Sig_vec                               *
!     iRoot   : integer                                                *
!               root number                                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
#include "rasdim.fh"
#include "davctl.fh"

RecNo = 0
if (itype == 1) then
  H_diag_RecNo = 1
  RecNo = H_diag_RecNo
else if (itype == 2) then
  CI_vec_RecNo = 1+PageNo(iRoot)
  RecNo = CI_vec_RecNo
else if (itype == 3) then
  Sig_vec_RecNo = 1+nkeep+PageNo(iRoot)
  RecNo = Sig_vec_RecNo
else if (itype == 4) then
  tmp_CI_vec_RecNo = 1+2*nKeep+iRoot
  RecNo = tmp_CI_vec_RecNo
else if (itype == 5) then
  tmp_Sig_vec_RecNo = 1+2*nKeep+n_Roots+iRoot
  RecNo = tmp_Sig_vec_RecNo
else
  write(6,*) 'RecNo: itype does not match'
  write(6,*) 'itype = ',itype
  call Abend()
end if

return

end function RecNo