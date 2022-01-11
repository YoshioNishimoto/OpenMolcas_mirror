#define _NEWCODE_
#ifdef _NEWCODE_
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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
      Subroutine xS12g(mGrid,Coeff,nD,F_xc)
      use xc_f03_lib_m
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vSigma
      use libxc
      implicit none
      integer :: mGrid, nD, nRho
      Real*8 :: F_xc(mGrid)
      Real*8 :: Coeff

      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! PBE correlation
      integer*4, parameter :: func_id = 495

      nRho=SIZE(Rho,1)
      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=2.00D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=4.00D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=0.50D0*vSigma(:,1:mGrid)
      End If

      call libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

      call xc_f03_func_end(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=0.50D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=0.25D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=2.00D0*vSigma(:,1:mGrid)
      End If

      Return

      End Subroutine xS12g
      Subroutine xS12h(mGrid,Coeff,nD,F_xc)
      use xc_f03_lib_m
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vSigma
      use libxc
      implicit none
      integer :: mGrid, nD, nRho
      Real*8 :: F_xc(mGrid)
      Real*8 :: Coeff

      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! PBE correlation
      integer*4, parameter :: func_id = 496

      nRho=SIZE(Rho,1)
      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=2.00D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=4.00D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=0.50D0*vSigma(:,1:mGrid)
      End If

      call libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

      call xc_f03_func_end(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=0.50D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=0.25D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=2.00D0*vSigma(:,1:mGrid)
      End If

      Return

    End Subroutine xS12h
#else
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
      Subroutine xS12g(mGrid,Coeff,iSpin,F_xc)
      use nq_Grid, only: vSigma
      use nq_Grid, only: vRho
      Implicit None
      Integer mGrid, iSpin
      Real*8 Coeff, F_xc(mGrid)
      integer, parameter :: gh_switch=1
      Call xS12gh(mGrid,Coeff,iSpin,F_xc,gh_switch)
      End Subroutine xS12g
      Subroutine xS12h(mGrid,Coeff,iSpin,F_xc)
      use nq_Grid, only: vSigma
      use nq_Grid, only: vRho
      Implicit None
      Integer mGrid, iSpin
      Real*8 Coeff, F_xc(mGrid)
      integer, parameter :: gh_switch=2
      Call xS12gh(mGrid,Coeff,iSpin,F_xc,gh_switch)
      End Subroutine xS12h
      Subroutine xS12gh(mGrid,Coeff,iSpin,F_xc,gh_switch)
!***********************************************************************
!                                                                      *
! Object S12g from Marcel Swart                                        *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
!             University of Lund, SWEDEN. June 2006  (B88 code)        *
!             Grigory A. Shamov, U of Manitoba, Dec 2009               *
!***********************************************************************
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vRho, vSigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20
      integer gh_switch

! IDORD=Order of derivatives to request from XPBE:
      idord=1
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

         call xS12g_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2,gh_switch)
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

         call xS12g_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2,gh_switch)

         sigmabb=Sigma(3,iGrid)
         call xS12g_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2,gh_switch)

!         write(6,*) rhoa,rhob,sigmaaa,sigmabb,Fa,Fb

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

      subroutine xS12g_(idord,rho_s,gamma_s,B88,dB88dr,dB88dg,d2B88dr2,d2B88drdg,d2B88dg2,gh_switch)
      implicit real*8 (a-h,o-z)
      parameter(third=1.0d0/3.0d0)
      parameter(four3=4.0d0/3.0d0)
      parameter(seven3=7.0d0/3.0d0)
!      parameter(bcoef=0.0042d0)
!     parameter(xldacff=0.930525736349100025D0)
      integer gh_switch
      real*8 gamma,rho,gdenom, hdenom

      parameter(b=1.0d0/137.0d0)

! initialized to Zero to avoid nasty complains from compilers
        rA = 0.0d0
        rK = 0.0d0
        rB = 0.0d0
        rC = 0.0d0
        rD = 0.0d0
        rE = 0.0d0

      if(gh_switch.eq.1) then
! GGA non-hybrid parameter set
        rA = 1.03842032d0
        rK = 0.757d0
        rB = 1.0d0 + rK - rA
        rC = 0.00403198d0
        rD = 0.00104596d0
        rE = 0.00594635d0
      elseif(gh_switch.eq.2) then
! GGA hybrid parameter set
        rA = 1.02543951d0
        rK = 0.757d0
        rB = 1.0d0 + rK - rA
        rC = 0.00761554d0
        rD = 0.00211063d0
        rE = 0.00604672d0
      endif

      C = -(1.5d0)*(0.75d0/acos(-1.0d0))**(third)
!     C = -0.9305257277176876


!      rho = min(rho_s , 1.0D-16 )
         rho=rho_s
!      gamma = min(gamma_s, 1.0D-16 )
      gamma=gamma_s
      rho13 = rho**third
      rho43 = rho**four3
! lda part:
!     xlda=-xldacff*rho43
! Note: Use x=sqrt(gamma)/rho**four3
      x = sqrt(gamma)/rho43
      x2 = x*x

!      hgi = 0.50d0/gamma

      gdenom = 1.0d0 + rC*x2 + rD*x2*x2
      hdenom = 1.0d0 + rE*x2
      ums = 1.0d0-(1.0d0/gdenom)
      vms = 1.0d0-(1.0d0/hdenom)
      g = rB*ums*vms
!
      dudx = (2.0d0*rC*x + 4.0d0*rD*x2*x)/(gdenom**2.0d0)
      dvdx = 2.0d0*rE*x/(hdenom**2.0d0)
      dg = C*rB*(dudx*vms + ums*dvdx)


      B88 =  C*rho43*(rA + g)

      if(idord.lt.1) goto 99

      dB88Dr = rA*(4d0/3d0)*rho13*C + (4d0/3d0)*rho13*(C*g-x*dg)

      t = dg / dsqrt(gamma)
      dB88Dg = t * 0.5d0

      if(idord.lt.2) goto 99

      write(6,*) 'S12g 2nd derivs not programmed'
      Call Abend()
!      d2B88Dr2 = (-1.d+1)*b*rho**((-8.d+0)/3.d+0)*gamma**(3.d+0/4.d+0)/9
!     1   .d+0
!      d2B88DrDg = b*rho**((-5.d+0)/3.d+0)*gamma**((-1.d+0)/4.d+0
!     2   )/2.d+0
!      d2B88Dg2 = 3.d+0*b*rho**((-2.d+0)/3.d+0)*gamma**((-5.d+0
!     3   )/4.d+0)/1.6d+1


  99  continue
      return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(d2B88dr2)
         Call Unused_real(d2B88drdg)
         Call Unused_real(d2B88dg2)
      End If
      end
#endif