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
! Stuff for Cholesky vector subtraction screening (subscr):
Module ChoSubScr
Implicit None
Private
Public:: Cho_SScreen, SSTau, SubScrStat, DSubScr, ip_DSPNm, l_DSPNm, SSNorm
Logical Cho_SScreen
Real*8  SSTau, SubScrStat(2)
Integer ip_DSPNm
Integer l_DSPNm
Character(LEN=3) SSNorm

Real*8, Allocatable:: DSubScr(:)
End Module ChoSubScr
