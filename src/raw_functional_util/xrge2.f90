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
      Subroutine xRGE2(mGrid,Coeff,nD,F_xc)
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

      ! Becke exchange 88 exchange
      integer*4, parameter :: func_id = 142

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

    End Subroutine xRGE2
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
!               2010, Grigory A. Shamov                                *
!***********************************************************************
      Subroutine XRGE2(mGrid,Coeff,iSpin,F_xc)
!***********************************************************************
!                                                                      *
! Object: To compute the exchange part of RGE2, by                     *
!  Ruzsinszky, Csonka, Scuseria, JCTC 2009, 5, 763-769                 *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
!             University of Lund, SWEDEN. June 2006                    *
!             Grigory A Shamov, U of Manitoba, spring 2010             *
!***********************************************************************
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vRho, vSigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20

! IDORD=Order of derivatives to request from XPBE:
      idord=1


      if (ispin.eq.1) then
! ispin=1 means spin zero.
! T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(1,iGrid))
         if(rhoa.lt.Ta) goto 110
         sigmaaa=Sigma(1,iGrid)

         call testRGE2_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)
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
         rhoa=max(1.0D-24,Rho(1,iGrid))
         rhob=max(1.0D-24,Rho(2,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) goto 210
         sigmaaa=Sigma(1,iGrid)
         call testRGE2_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         sigmabb=Sigma(3,iGrid)
         call testRGE2_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

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

! Cmu to be modified to 10/81 to make it PBEsol
!      Data Ckp, Cmu / 0.804D0, 0.12345679012346D0/


      Subroutine testXPBE_(idord,rho_s,sigma_s,f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2)
!***********************************************************************
!                                                                      *
! Object: To compute the PBE exchange functional from MOLPRO manual    *
!         using Maxima derivatives and optimization/fortran generator  *
!  Ref: PBE,"Generalized gradient approximation made simple,"          *
!      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!      Author: Grigory Shamov, U of Manitoba, Feb-March 2010           *
!***********************************************************************
      Implicit None
      Integer idord
      real*8 rho_s,sigma_s, rho, gamma, k, mu, Tpi, Cx, one
      real*8 f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2
      parameter(Tpi=3.141592653589793d0,Cx=0.9305257363491d0, one=1.0d0 )

      Data k, mu / 0.804D0, 0.2195149727645171D0/
!      Data CkF / 3.0936677262801359310D0/
!      Data CeX /-0.73855876638202240588D0/
      real*8 T1,T2,T3,T4,T5,T6 ,T7, T8
!
      rho=max(1.0D-24,rho_s)
      gamma=max(1.0D-24,sigma_s)


      T1=rho**(4.d+0/3.d+0)
      T2=one/k
      T3=1.6455307846020562d-2*T2*mu*gamma/rho**(8.d+0/3.d+0)+1.d+0
      T4=-1.d+0*k/T3+k+1.d+0
      T5=one/rho**(7.d+0/3.d+0)
      T6=one/T3**2
      T7=mu**2
      T8=one/T3**3

      F = -9.305257363491001d-1*T1*T4

      dFdr = 4.083223320071841d-2*mu*T5*gamma*T6-1.2407009817988002d+0*rho**(1.d+0/3.d+0)*T4

      dFdg = -1.5312087450269407d-2*mu*T6/T1

      d2Fdr2 = -4.135669939329333d-1*T4/rho**(2.d+0/3.d+0)-4.0832233200718443d-2*mu*gamma*  &
         T6/rho**(1.d+1/3.d+0)+3.583503825911055d-3*T2*T7*gamma**2*T8/rho**6

      d2Fdrdg = 2.0416116600359208d-2*mu*T5*T6-1.343813934716646d-3*T2*T7*gamma*T8/rho**5

      d2Fdg2 = 5.03930225518742d-4*T2*T7*T8/rho**4

      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(idord)
      End

      Subroutine testRGE2_(idord,rho_s,sigma_s,f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2)
!***********************************************************************
!                                                                      *
! Object: To compute the RGE exchange functional. Procedure: the PBE   *
!  exchange functional from MOLPRO manual with Fx changed to RGE, the  *
!  mu parameter changed to be as of PBEsol. Code was made              *
!         using Maxima derivatives and optimization/fortran generator  *
!                                                                      *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!      Author: Grigory Shamov, U of Manitoba, Feb-March 2010           *
!***********************************************************************

      Implicit None
      Integer idord
      real*8 rho_s,sigma_s, rho, gamma, k, mu, Tpi, Cx, one
      real*8 f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2
      parameter(Tpi=3.141592653589793d0,Cx=0.9305257363491d0, one=1.0d0 )

! Cmu to be modified to 10/81 to make it PBEsol
      Data k, mu / 0.804D0, 0.12345679012346D0/

!      Data k, mu / 0.804D0, 0.2195149727645171D0/
!      Data CkF / 3.0936677262801359310D0/
!      Data CeX /-0.73855876638202240588D0/
      real*8 T1, T2, T3, T4, T5, T6
      real*8 T7, T8, T9, T10, T11
      real*8 T12, T13, T14, T15, T16
!
      rho=max(1.0D-24,rho_s)
      gamma=max(1.0D-24,sigma_s)

      T1= rho**(4.d+0/3.d+0)
      T2=one/k
      T3=one/rho**(8.d+0/3.d+0)
      T4=one/k**2
      T5=  mu**2
      T6=one/rho**(1.6d+1/3.d+0)
      T7=gamma**2
      T8=2.7077715630730587d-4*T4*T5*T6*T7+1.6455307846020562d-2*T2*mu*T3*gamma+1.d+0
      T9=  -1.d+0*k/T8+k+1.d+0
      T10=one/rho**(1.1000000000000001d+1/3.d+0)
      T11=one/rho**(1.9d+1/3.d+0)
      T12=-1.4441448336389645d-3*T4*T5*T11*T7-4.388082092272149d-2*T2*mu*T10*gamma
      T13=one/T8**2
      T14=rho**(1.d+0/3.d+0)
      T15=5.415543126146117d-4*T4*T5*T6*gamma+1.6455307846020562d-2*T2*mu*T3
      T16=one/T8**3

      F = -9.305257363491001d-1*T1*T9

      dFdr = -1.2407009817988002d+0*T14*T9-9.305257363491001d-1*k*T1*T12*T13
      dFdg = -9.305257363491001d-1*k*T1*T15*T13

      d2Fdr2 = -4.135669939329333d-1*T9/rho**(2.d+0/3.d+0)-2.4814019635975995d+0*k                          &
         *T14*T12*T13-9.305257363491001d-1*k*T1*(9.146250613046776d-3*T4                                    &
         *T5*T7/rho**(2.2000000000000003d+1/3.d+0)+1.6089634338331216d-1                                    &
         *T2*mu*gamma/rho**(1.3999999999999999d+1/3.d+0))*T13+1.8610514726982003d+0*k*T1*T12**2*T16

      d2Fdrdg = -1.2407009817988002d+0*k*                                                                   &
         T14*T15*T13-9.305257363491001d-1*k*T1*(-2.888289667277929d-3*T4                                    &
         *T5*T11*gamma-4.388082092272149d-2*T2*mu*T10)*T13+1.8610514726982003d+0*k*T1*T15*T12*T16

      d2Fdg2 = 1.8610514726982003d+0*k*T1*T15**2*T16-5.03930225518742d-4*T2*T5*T13/rho**4

      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(idord)
      End
#endif