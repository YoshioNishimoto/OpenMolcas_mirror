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
      SubRoutine Cho_P_SetPass(Diag,Sync,DiaSh,iSySh,iLoc,Conv,nPotSh)
!
!     Purpose: check convergence and, if not converged, set up next
!              next integral pass
!              Diag is the local diagonal; the global diagonal is
!              synchronized if Sync=.True. Note that DiaSh and iSySh
!              must be allocated with dimension nnShl_G (global number
!              of shell pairs in 1st reduced set). iLoc is the location
!              to use in the reduced set index arrays. On exit,
!              Conv=.True. if converged and nPotSh is the number of
!              shell pairs whose max. diagonal element is larger than
!              the decomposition threshold.
!
      use ChoSwp, only: Diag_G
      Implicit None
      Real*8  Diag(*)
      Logical Sync, Conv
      Real*8  DiaSh(*)
      Integer iSySh(*)
      Integer iLoc, nPotSh
#include "cho_para_info.fh"
#include "choglob.fh"

      If (Cho_Real_Par) Then

!        Sync diagonal if requested.
!        ---------------------------

         If (Sync) Then
            Call Cho_P_SyncDiag(Diag,iLoc)
         End If

!        Swap local and global index arrays and set next integral pass
!        original serial routine.
!        -------------------------------------------------------------

         Call Cho_P_IndxSwp()
         Call Cho_SetPass(Diag_G,DiaSh,iSySh,iLoc,Conv,nPotSh)
         Call Cho_P_IndxSwp()

      Else

         Call Cho_SetPass(Diag,DiaSh,iSySh,iLoc,Conv,nPotSh)

      End If

      End
