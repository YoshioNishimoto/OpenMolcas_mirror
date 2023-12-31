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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine vAssmbl(Rnxyz,Axyz,la,Rxyz,lr,Bxyz,lb,nZeta,HerW,nHer,
     &                  Temp)
!***********************************************************************
!                                                                      *
! Object: to assemble the cartesian components of the multipole moment *
!         matrix within the framework of the Gauss-Hermite quadrature. *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************
      use Constants
      Implicit None
      Integer la,lr,lb,nZeta,nHer
      Real*8 Rnxyz(nZeta*3,0:la,0:lb,0:lr), HerW(nHer),
     &       Axyz(nZeta*3,nHer,0:la),
     &       Rxyz(nZeta*3,nHer,0:lr),
     &       Bxyz(nZeta*3,nHer,0:lb), Temp(nZeta*3,nHer)
#ifdef _DEBUGPRINT_
      Character(LEN=80) Label
#endif
      Integer ia, ib, iHer, iZCar, ir
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' In vAssmbl:HerW',' ',HerW,1,nHer)
      Call RecPrt(' In vAssmbl:Axyz',' ',Axyz,nZeta*3,nHer*(la+1))
      Call RecPrt(' In vAssmbl:Bxyz',' ',Bxyz,nZeta*3,nHer*(lb+1))
      Call RecPrt(' In vAssmbl:Rxyz',' ',Rxyz,nZeta*3,nHer*(lr+1))
#endif
!
!
      rnxyz(:,:,:,:)=Zero
      Do 100 ia = 0, la
         Do 110 ib = 0, lb
            Do 111 iHer = 1, nHer
               Do 112 iZCar = 1, 3*nZeta
                  Temp(iZCar,iHer) =  Axyz(iZCar,iHer,ia)*
     &                                Bxyz(iZCar,iHer,ib)*HerW(iHer)
 112           Continue
 111        Continue
            Do 120 ir = 0, lr
!
!     Generate the cartesian components of the multipole moment
!     matrix as a sum of the value of the integrand, evaluated
!     at a root, times a weight.
!
               Do 30 iHer = 1, nHer
                  Do 10 iZCar = 1, 3*nZeta
                     Rnxyz(iZCar,ia,ib,ir) = Rnxyz(iZCar,ia,ib,ir) +
     &                             Temp(iZCar,iHer)*
     &                             Rxyz(iZCar,iHer,ir)
 10               Continue
 30            Continue
!
#ifdef _DEBUGPRINT_
               Write (Label,'(A,I2,A,I2,A,I2,A)')
     &         ' In vAssmbl: Rnxyz(',ia,',',ib,',',ir,')'
               Call RecPrt(Label,' ',Rnxyz(1,ia,ib,ir),nZeta,3)
#endif
 120        Continue
 110     Continue
 100  Continue
!
      End SubRoutine vAssmbl
