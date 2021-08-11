************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine NACInt(xyz,nCent,H12,Bf,lWrite_,Label,dBf,ldB,lIter)
      use Slapaf_Info, only: NAC
      Implicit Real*8  (a-h,o-z)
#include "real.fh"
#include "constants.fh"
      Real*8   Bf(3,nCent), xyz(3,nCent), dBf(3*nCent,3*nCent)
      Logical lWrite_, ldB
      Character*8 Label
*
*
      H12=Zero
      If (lWrite_) Then
         Write (6,'(2A,F18.8,A,F18.8,A)')
     &             Label,' : H12               = ',
     &             H12, ' hartree '
      End If
*
*---- Compute the WDC B-matrix
*
      Do iCent = 1, nCent
         Fact=DBLE(iDeg(xyz(1,iCent)))
         Do iCar = 1, 3
            Bf(iCar,iCent)=Fact*NAC(iCar,iCent,lIter)
         End Do
      End Do
#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'NACInt, lIter:',lIter
      Call RecPrt('NAC',' ',NAC(1,1,lIter),3,nCent)
      Call RecPrt('Bf',' ',Bf,3,nCent)
#endif
*
*---- Compute the cartesian derivative of the B-Matrix.
*
      If (ldB) Then
         dBf(:,:)=Zero
*        Do i = 1, 3*nCent
*           dBf(i,i)=0.1D-1
*        End Do
C        Call RecPrt('dBf','(9F9.1)',dBf,3*nCent,3*nCent)
*
      End If
*
      Return
      End
