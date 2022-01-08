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
! Copyright (C) 2006, Per Ake Malmqvist                                *
!               2009, Grigory A. Shamov                                *
!***********************************************************************
      Subroutine xB86(mGrid,Coeff,iSpin,F_xc)
!***********************************************************************
!                                                                      *
! Object: Becke86 exchange functional, from Appendix of Becke's B07 ref*
!        Becke, Johnson, J.Chem.Phys 127, 124108 (2007)                *
!        derived with Maxima                                           *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!      Author:template by Per Ake Malmqvist,                           *
!             University of Lund, SWEDEN. June 2006                    *
!             Grigory A Shamov, U of Manitoba, winter 2009             *
!***********************************************************************
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vRho, vSigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20


! IDORD=Order of derivatives to request from XPBE:
      idord=1
!
      Rho_Min=T_X*1.0D-2
!
      if (ispin.eq.1) then
! ispin=1 means spin zero.

! T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         rhoa=rho(1,iGrid)
         if(rhoa.lt.Ta) goto 110
         sigmaaa=Sigma(1,iGrid)

         call xB86_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*(2.0D0*Fa)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
! Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
! Note: For xpbe, dFdgammaab is zero.
 110     continue
        end do

      else
! ispin .ne. 1, use both alpha and beta components.

        do iGrid=1,mgrid
         rhoa=Max(Rho_Min,rho(1,iGrid))
         rhob=Max(Rho_Min,rho(2,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) goto 210
         sigmaaa=Sigma(1,iGrid)
         call xB86_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         sigmabb=Sigma(3,iGrid)
         call xB86_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid)=F_xc(iGrid)+Coeff*(Fa+Fb)
         vRho(1,iGrid)=vRho(1,iGrid)+Coeff*dFdrhoa
         vRho(2,iGrid)=vRho(2,iGrid)+Coeff*dFdrhob
! Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
! Note: For xpbe, dFdgammaab is zero.
         vSigma(1,iGrid)=vSigma(1,iGrid)+Coeff*dFdgammaaa
         vSigma(3,iGrid)=vSigma(3,iGrid)+Coeff*dFdgammabb
 210     continue
        end do

      end if

      Return
      End

      subroutine xB86_(idord,rho_s,gamma_s,B88,dB88dr,dB88dg,d2B88dr2,d2B88drdg,d2B88dg2)
      implicit real*8 (a-h,o-z)
      parameter(third=1.0d0/3.0d0)
      parameter(four3=4.0d0/3.0d0)
      parameter(seven3=7.0d0/3.0d0)
      parameter(b=0.0036d0)
      parameter(g=0.0040d0)

!     parameter(xldacff=0.930525736349100025D0)

      rho=rho_s+1.0D-16
      gamma=gamma_s+1.0D-16

      r43 = rho**four3
! lda part:
!     xlda=-xldacff*r43

      B88 = -b*gamma/(g*gamma/r43 + r43)

      if(idord.lt.1) goto 99

      temp1=g*gamma/rho**(8.0E+0/3.0E+0)+1

      dB88dr = 4.0E+0*b*gamma/(3.0E+0*rho**(7.0E+0/3.0E+0)*temp1)+ (-8.0E+0)*b*g*gamma**2/(3.0E+0*rho**5*temp1**2)

      temp1=g*gamma/rho**(8.0E+0/3.0E+0)+1

      dB88dg = b*g*gamma/(rho**4*temp1**2) - b/(rho**(4.0E+0/3.0E+0)*temp1)

      if(idord.lt.2) goto 99

      T1=g**2
      T2=g*gamma/rho**(8.d+0/3.d+0)+1
      T3=1/T2**3
      T4=gamma**2
      T5=1/T2**2
      T6=1/T2

      d2B88Dr2 = (-2.7999999999999997d+1)*b*gamma*T6/(9.d+0*rho**(1.d+1/3.d+0))+1.52d+2*b* &
         g*T4*T5/(9.d+0*rho**6)+(-1.28d+2)*b*T1*gamma**3*T3/(9.d+0*rho**                   &
         (2.6d+1/3.d+0))

      d2B88DrDg = 4.d+0*b*T6/(3.d+0*rho**(7.d+0/3.d+0))+(-2.d+1)*b*g*gamma*T5/(3.d+0*rho**5) &
      +1.6d+1*b*T1*T4*T3/(3.d+0*rho**(2.2999999999999998d+1/3.d+0))

      d2B88Dg2 = 2*b*g*T5/rho**4-2*b*T1*gamma*T3/rho**(2.d+1/3.d+0)


  99  continue

      return
      end
