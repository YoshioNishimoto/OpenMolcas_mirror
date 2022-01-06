************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2006, Per Ake Malmqvist                                *
*               2010, Grigory A. Shamov                                *
************************************************************************
      Subroutine CWIG(mGrid,Coeff,iSpin,F_xc)
************************************************************************
*                                                                      *
* Object: To compute the Powerful Wigner correlation functional        *
* Functional Repository (http://www.cse.clrc.ac.uk/qcg/dft)            *
*    Steward PA, Gill PW. J.Chem.Soc. Fadaday Trans. 91 (1995) 4337    *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. June 2006                    *
************************************************************************
      use nq_Grid, only: Rho
      use nq_Grid, only: vRho
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20

* IDORD=Order of derivatives to request from XPBE:
      idord=1


      if (ispin.eq.1) then
* ispin=1 means spin zero.
* T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(1,iGrid))
         if(rhoa.lt.Ta) goto 110

         call cWIG_(idord,rhoa,rhoa,F,dFdrhoa,dFdrhob,
     &          d2Fdra2,d2Fdrb2,d2Fdrab)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*(F)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
C     F is ok, but dFdrho comes in two pieces ? nope
* Note: For xpbe, dFdgammaab is zero.
 110     continue
        end do
      else
* ispin .ne. 1, use both alpha and beta components.
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(1,iGrid))
         rhob=max(1.0D-24,Rho(2,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) goto 210
         call cWIG_(idord,rhoa,rhob,F,dFdrhoa,dFdrhob,
     &          d2Fdra2,d2Fdrb2,d2Fdrab)

         F_xc(iGrid)=F_xc(iGrid)+Coeff*(F)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
         vRho(2,iGrid)=vRho(2,iGrid)+Coeff*dFdrhob

 210     continue
        end do
      end if

      Return
      End

      Subroutine cWIG_(idord,rho_a,rho_b, F
     &                ,dFdra,dFdrb,
     &                 d2Fdra2,d2Frb2,d2Fdrab)
************************************************************************
*                                                                      *
* Object: The Powerful Wigner Functional                               *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author: Grigory A Shamov, U. of Manitoba. Winter 2009           *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      parameter(one=1.0d0)
      parameter(four=4.0d0)
      parameter(eight=4.0d0)

      parameter(d= 0.349d0)
      parameter(a= 0.04918d0)

      rhoa=max(1.0D-24,rho_a)
      rhob=max(1.0D-24,rho_b)


      temp1=rhob+rhoa
      temp2=one/temp1
      temp3=d/temp1**(1.0d+0/3.0d+0)+1
      temp4=one/temp3
      temp5=(-4.0d+0)*a*d*rhoa*rhob/(3.0d+0*temp1**(7.0d+0/
     &   3.0d+0)*temp3**2)
      temp6=4*a*rhoa*rhob*temp4/temp1**2
      f = -4*a*rhoa*rhob*temp2*temp4


      if (idord .lt. 1) goto 99

      dFdra = -4*a*rhob*temp2*temp4+temp6+temp5
      dFdrb = -4*a*rhoa*temp2*temp4+temp6+temp5

      if ( idord .lt. 2) goto 99

      Call WarningMessage(2,
     & 'Second derivatives werent implemented for Wigner functional')
      call Abend()

99    continue

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(d2Fdra2)
         Call Unused_real(d2Frb2)
         Call Unused_real(d2Fdrab)
      End If
      End
