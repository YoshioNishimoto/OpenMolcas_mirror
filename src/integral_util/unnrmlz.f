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
      SubRoutine UnNrmlz(Exp,nPrim,Coeff,nCntrc,iAng)
      use Constants, only: TwoP34, Two, Three, Four
      Implicit None
      Integer nPrim, nCntrc, iAng
      Real*8 Exp(nPrim), Coeff(nPrim,nCntrc)

      Integer i, j
!
      Do  i = 1, nCntrc
         Do  j = 1, nPrim
            Coeff(j,i) = Coeff(j,i) /( TwoP34 *
     &                (Four*Exp(j))**((Two*DBLE(iAng)+Three)/Four))
         End Do
      End Do

      End SubRoutine UnNrmlz
