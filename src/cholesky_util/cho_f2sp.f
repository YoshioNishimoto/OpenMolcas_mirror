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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************
      Integer Function Cho_F2SP(iSP)
!
!     Thomas Bondo Pedersen, March 2006.
!
!     Purpose: return reduced shell pair index for full shell pair index
!              iSP. If not found, 0 is returned.
!              Note: nnShl_SP is used to avoid problems in parallel runs
!              when swapping nnShl and nnShl_G. If properly set,
!              nnShl_SP = nnShl_G.
!
      use ChoArr, only: iSP2F
      use ChoSP, only: nnShl_SP
      Implicit None
      Integer iSP

      Integer jSP

      Cho_F2SP = 0
      Do jSP = 1,nnShl_SP
         If (iSP2F(jSP) .eq. iSP) Then
            Cho_F2SP = jSP
            Return
         End If
      End Do

      End
