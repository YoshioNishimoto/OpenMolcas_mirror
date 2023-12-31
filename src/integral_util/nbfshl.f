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
!***********************************************************************
!                                                                      *
!  Integer Function NbfShl    returns # of bf for shell,symmetry       *
!                                                                      *
!***********************************************************************
!----------------------------------------------------------------------
      Integer Function nbfshl(iSkal,irp)
      use iSD_data, only: iSD
      use SOAO_Info, only: iAOtSO
!----------------------------------------------------------------------
      Implicit None
      Integer iSkal, irp

      Integer iAO, iCmp, i
!
!  returns number of basis functions for given shell and symmetry
!
      nbfshl=0
      iAO    = iSD( 7,iSkal)
      iCmp   = iSD( 2,iSkal)
!     loop over components of shell...
      Do i=1, iCmp
         If (iAOtSO(iAO+i,irp)>0) nbfshl = nbfshl + iSD(3,iSkal)
      End Do

      return
      End Function nbfshl
