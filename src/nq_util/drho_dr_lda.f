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
* Copyright (C) 2000,2002, Roland Lindh                                *
************************************************************************
#define _NEWCODE_
#ifdef _NEWCODE_
      Subroutine dRho_dR_LDA(nD,dRho_dR,ndRho_dR,mGrid,
     &                       list_s,nlist_s,
     &                       QabAO,ipTabAO,mAO,nTabAO,
     &                       nGrad_Eff,list_g,
     &                       Fact,ndc,
     &                       list_bas,Index,nIndex)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
*                                                                      *
*             Modified May-June 2002 in Tokyo for DFT gradient by RL.  *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      use nq_Grid, only: Grid_AO, Dens_AO, Ind_Grd, TabAO
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "disp.fh"
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer On, Off
      Parameter (On=1, Off=0)
      Integer list_s(2,nlist_s), list_g(3,nlist_s),
     &        ipTabAO(nlist_s), list_bas(2,nlist_s), Index(nIndex)
      Real*8 QabAO(nTabAO), Fact(ndc**2),
     &       dRho_dR(ndRho_dR,mGrid,nGrad_Eff)
      Integer IndGrd_Eff(3,2)
      Integer ipD(2)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      crap = QabAO(1)
      iCrap=ipTAbAO(1)
      dRho_dR(:,:,:)=Zero
      nAO = SIZE(Dens_AO,1)
      Dens_AO(:,:,:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
      iOff = 0
      Do ilist_s=1,nlist_s
         iS      =list_s(1,ilist_s)
         iCmp       =iSD( 2,iS)
         iBas_Eff  = list_bas(1,ilist_s)
         iBas       =iSD( 3,iS)
         kDCRE   =list_s(2,ilist_s)
         iShell     =iSD(11,iS)
         mdci       =iSD(10,iS)
         index_i    =list_bas(2,ilist_s)
         nFunc_i=iBas*iCmp
         n_iBas=iBas_Eff*iCmp

*
         Call ICopy(3,list_g(1,ilist_s),1,IndGrd_Eff(1,1),1)
         Indx = List_g(1,ilist_s)
         Ind_Grd(1,iOff+1:iOff+n_iBas)=Indx
         Indy = List_g(2,ilist_s)
         Ind_Grd(2,iOff+1:iOff+n_iBas)=Indy
         Indz = List_g(3,ilist_s)
         Ind_Grd(3,iOff+1:iOff+n_iBas)=Indz
*
         jOff = 0
         Do jlist_s=1,ilist_s
            jS      =list_s(1,jlist_s)
            kDCRR   =list_s(2,jlist_s)
            jCmp       =iSD( 2,jS)
            jBas       =iSD( 3,jS)
            jBas_Eff   =list_bas(1,jlist_s)
            mdcj       =iSD(10,jS)
            jShell     =iSD(11,jS)
            index_j    =list_bas(2,jlist_s)
            nFuncji=jBas*jCmp
            n_jBas=jBas_Eff*jCmp
*
            Call ICopy(3,list_g(1,jlist_s),1,IndGrd_Eff(1,2),1)
*
            mDij=iBas*jBas*iCmp*jCmp
*
*---------- Get the Density
*
            ijS=iTri(iShell,jShell)
            ipTmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ipTmp,nD)
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*)
               Write (6,*) 'iS,jS=',iS,jS
               Write (6,*) 'mDCRij,mDij=',mDCRij,mDij
               Write (6,*) 'ipDij,ipDSij,ipDDij=',ipDij,ipDSij,ipDDij
            End If
#endif
*
            ij = (mdcj-1)*ndc + mdci
*
            iER=iEOr(kDCRE,kDCRR)
            lDCRER=NrOpr(iER)
*
            ip_D_a=ipDij+lDCRER*mDij
            ip_D_b=ip_D_a
            If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
            ipD(1)=ip_D_a
            ipD(2)=ip_D_b
*
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*) 'dRho_dR_LDA'
               nBB = iBas*jBas
               nCC = iCmp*jCmp
               Write (6,*) 'iShell,jShell=',iShell,jShell
               Write (6,*) 'kDCRE,kDCRR=',kDCRE,kDCRR
               Write (6,*) 'iER,lDCRER=',iER,lDCRER
               Call RecPrt('DAij',' ',DeDe(ip_D_a),nBB,nCC)
               If (nD.ne.1)
     &            Call RecPrt('DBij',' ',DeDe(ip_D_b),nBB,nCC)
            End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
            Do iD = 1, nD
               Do j_R = 1, n_jBas
                  jCB = Index(index_j-1+j_R)    ! Real index
                  j_A = j_R + jOff            ! Absolute index
*
                  Do i_R = 1, n_iBas
                     iCB = Index(index_i-1+i_R)
                     i_A = i_R + iOff
*
                     ij_D = (jCB-1)*nFunc_i + iCB - 1
                     DAij =DeDe(ipD(iD)+ij_D)*Fact(ij)
                     Dens_AO(i_A,j_A,iD) = DAij
                     Dens_AO(j_A,i_A,iD) = DAij
*
                  End Do          ! iCB
               End Do             ! jCB
            End Do
*                                                                      *
************************************************************************
*                                                                      *
            jOff = jOff + jBas_Eff*jCmp
         End Do                      ! jlist_s
         iOff = iOff + iBas_Eff*iCmp
      End Do                         ! ilist_s
*     Write (*,*) 'ddot',ddot_(SIZE(Dens_AO),Dens_AO,1,One,0)
*                                                                      *
************************************************************************
*                                                                      *
*   not that this is redundant and already done in mk_Rho!!!!
      Call DGEMM_('N','N',mAO*mGrid,nAO*nD,nAO,
     &            One,TabAO,mAO*mGrid,
     &                Dens_AO,nAO,
     &            Zero,Grid_AO,mAO*mGrid)
*     Write (*,*) 'gdot',ddot_(SIZE(Grid_AO),Grid_AO,1,One,0)

      Do iD = 1, nD
         Do iAO = 1, nAO
*
*---------- Loop over cartesian components
*
            Do iCar = 1, 3

               Ind_xyz=Ind_Grd(iCar,iAO)
               j = iCar + 1

               If (Ind_xyz/=0) Then
                  Do iGrid = 1, mGrid
                     dRho_dR(iD,iGrid,Ind_xyz)=dRho_dR(iD,iGrid,Ind_xyz)
     &                             + Two * Grid_AO(1,iGrid,iAO,iD)
     &                             * TabAO(j,iGrid,iAO)
                  End Do
               End If

            End Do
         End Do
      End Do

*     Do iAO = 1, nAO
*        Do iCar = 1, 3
*           Ind_xyz = Ind_Grd(iCar,iAO)
*           If (Ind_xyz==9) Then
*              Write (*,*) 'iAO,iCar,Ind_xyz=',iAO,iCar,Ind_xyz
*              Write (*,*) Grid_AO(:,mGrid,iAO,1)
*           End If
*        End Do
*     End Do
*     Write (6,*) dRho_dR(1,mGrid,9)

*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('dRho_dR_LDA: dRho_dR',' ',dRho_dR,
     &                       ndRho_dR*mGrid,nGrad_Eff)
*
#endif
*     Stop 123
      Return
      End
#else
      Subroutine dRho_dR_LDA(nD,dRho_dR,ndRho_dR,mGrid,
     &                       list_s,nlist_s,
     &                       TabAO,ipTabAO,mAO,nTabAO,
     &                       nGrad_Eff,list_g,
     &                       Fact,ndc,
     &                       list_bas,Index,nIndex)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
*                                                                      *
*             Modified May-June 2002 in Tokyo for DFT gradient by RL.  *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "disp.fh"
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer On, Off
      Parameter (On=1, Off=0)
      Integer list_s(2,nlist_s), list_g(3,nlist_s),
     &        ipTabAO(nlist_s), list_bas(2,nlist_s), Index(nIndex)
      Real*8 TabAO(nTabAO), Fact(ndc**2),
     &       dRho_dR(ndRho_dR,mGrid,nGrad_Eff)
      Integer IndGrd_Eff(3,2)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (Debug) Then
         Call RecPrt('dRho_dR_LDA:Dens',' ',Dens,nDens,nD)
         Write (6,*) 'mAO=',mAO
         Write (6,*) 'mGrid=',mGrid
         Write (6,*) 'nTabAO=',nTabAO
         Write (6,*) 'nlist_s=',nlist_s
         Write (6,*) 'nD=',nD
         Do iList_s = 1, nList_s
            Write (6,*) 'iList_s=',iList_s
            iS = list_s(1,ilist_s)
            iCmp  = iSD( 2,iS)
            iBas_Eff  = list_bas(1,ilist_s)
            mTabAO=iBas_Eff*iCmp
            mdci  = iSD(10,iS)
            Call RecPrt('dRho_dR_LDA: TabAO',' ',
     &                  TabAO(ipTabAO(iList_s)),mAO,mGrid*mTabAO)
         End Do
      End If
#endif
*
      dRho_dR(:,:,:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
      Do ilist_s=1,nlist_s
         iS      =list_s(1,ilist_s)
         iCmp       =iSD( 2,iS)
         iBas_Eff  = list_bas(1,ilist_s)
         iBas       =iSD( 3,iS)
         kDCRE   =list_s(2,ilist_s)
         iShell     =iSD(11,iS)
         mdci       =iSD(10,iS)
         index_i    =list_bas(2,ilist_s)
         nFunc_i=iBas*iCmp
         n_iBas=iBas_Eff*iCmp

*
         Call ICopy(3,list_g(1,ilist_s),1,IndGrd_Eff(1,1),1)
*
         Do jlist_s=1,ilist_s
            jS      =list_s(1,jlist_s)
            kDCRR   =list_s(2,jlist_s)
            jCmp       =iSD( 2,jS)
            jBas       =iSD( 3,jS)
            jBas_Eff   =list_bas(1,jlist_s)
            mdcj       =iSD(10,jS)
            jShell     =iSD(11,jS)
            index_j    =list_bas(2,jlist_s)
            nFuncji=jBas*jCmp
            n_jBas=jBas_Eff*jCmp
*
            Call ICopy(3,list_g(1,jlist_s),1,IndGrd_Eff(1,2),1)
*
            mDij=iBas*jBas*iCmp*jCmp
*
*---------- Get the Density
*
            ijS=iTri(iShell,jShell)
            ipTmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ipTmp,nD)
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*)
               Write (6,*) 'iS,jS=',iS,jS
               Write (6,*) 'mDCRij,mDij=',mDCRij,mDij
               Write (6,*) 'ipDij,ipDSij,ipDDij=',ipDij,ipDSij,ipDDij
            End If
#endif
*
            ij = (mdcj-1)*ndc + mdci
*
            Fij=Two
            If (ilist_s.eq.jlist_s) Fij=One
*
            iER=iEOr(kDCRE,kDCRR)
            lDCRER=NrOpr(iER)
*
            ip_D_a=ipDij+lDCRER*mDij
            ip_D_b=ip_D_a
            If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij

*
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*) 'dRho_dR_LDA'
               nBB = iBas*jBas
               nCC = iCmp*jCmp
               Write (6,*) 'iShell,jShell=',iShell,jShell
               Write (6,*) 'kDCRE,kDCRR=',kDCRE,kDCRR
               Write (6,*) 'iER,lDCRER=',iER,lDCRER
               Call RecPrt('DAij',' ',DeDe(ip_D_a),nBB,nCC)
               If (nD.ne.1)
     &            Call RecPrt('DBij',' ',DeDe(ip_D_b),nBB,nCC)
            End If
#endif
            If (nD.eq.1) Then
               Call Do_Rho2da(dRho_dR,   mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Fij,IndGrd_Eff,
     &                       Index(index_i),Index(index_j))
            Else
               Call Do_Rho2d_(dRho_dR,   mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Fij,IndGrd_Eff,
     &                       Index(index_i),Index(index_j))
            End If
*
         End Do                      ! jlist_s
      End Do                         ! ilist_s
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('dRho_dR_LDA: dRho_dR',' ',dRho_dR,
     &                       ndRho_dR*mGrid,nGrad_Eff)
*
#endif
      Return
      End
      Subroutine Do_Rho2da(dRho_dR,   mGrid,nGrad_Eff,
     &                     DAij,          mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(   mGrid,nGrad_Eff), DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
      Integer Ind(3)
      Data Ind/2,3,4/
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB=Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
*
*---------- Loop over cartesian components
*
            Do iCar = 1, 3
               Ind_1=IndGrd_Eff(iCar,1)
               Ind_2=IndGrd_Eff(iCar,2)
               j = Ind(iCar)
               If (Ind_1.ne.0.and.Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(iGrid,Ind_1)=dRho_dR(iGrid,Ind_1)
     &                                   +Prod_1*DAij_
                     dRho_dR(iGrid,Ind_2)=dRho_dR(iGrid,Ind_2)
     &                                   +Prod_2*DAij_
                  End Do ! iGrid
               Else If (Ind_1.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     dRho_dR(iGrid,Ind_1)=dRho_dR(iGrid,Ind_1)
     &                                   +Prod_1*DAij_
                  End Do ! iGrid
               Else If (Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(iGrid,Ind_2)=dRho_dR(iGrid,Ind_2)
     &                                   +Prod_2*DAij_
                  End Do ! iGrid
               End If
*
            End Do  ! iCar
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
      Subroutine Do_Rho2d_(dRho_dR,   mGrid,nGrad_Eff,
     &                     DAij,DBij,     mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(   2,mGrid,nGrad_Eff),
     &       DAij(iBas*iCmp,jBas*jCmp), DBij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
      Integer Ind(3)
      Data Ind/2,3,4/
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB=Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            DBij_=DBij(iCB,jCB)*Fact
*
*---------- Loop over cartesian components
*
            Do iCar = 1, 3
               Ind_1=IndGrd_Eff(iCar,1)
               Ind_2=IndGrd_Eff(iCar,2)
               j = Ind(iCar)
               If (Ind_1.ne.0.and.Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,Ind_1)=dRho_dR(1,iGrid,Ind_1)
     &                                     +Prod_1*DAij_
                     dRho_dR(2,iGrid,Ind_1)=dRho_dR(2,iGrid,Ind_1)
     &                                     +Prod_1*DBij_
                     dRho_dR(1,iGrid,Ind_2)=dRho_dR(1,iGrid,Ind_2)
     &                                     +Prod_2*DAij_
                     dRho_dR(2,iGrid,Ind_2)=dRho_dR(2,iGrid,Ind_2)
     &                                     +Prod_2*DBij_
                  End Do ! iGrid
               Else If (Ind_1.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,Ind_1)=dRho_dR(1,iGrid,Ind_1)
     &                                     +Prod_1*DAij_
                     dRho_dR(2,iGrid,Ind_1)=dRho_dR(2,iGrid,Ind_1)
     &                                     +Prod_1*DBij_
                  End Do ! iGrid
               Else If (Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,Ind_2)=dRho_dR(1,iGrid,Ind_2)
     &                                     +Prod_2*DAij_
                     dRho_dR(2,iGrid,Ind_2)=dRho_dR(2,iGrid,Ind_2)
     &                                     +Prod_2*DBij_
                  End Do ! iGrid
               End If
*
            End Do ! iCar
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
#endif
