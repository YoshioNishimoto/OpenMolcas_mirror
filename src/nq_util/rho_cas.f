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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine Rho_CAS(Dens,nDens,nD,Rho,nRho,mGrid,
     &                   list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                   Fact,mdc,list_bas,Index,nIndex)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer list_s(2,nlist_s), ipTabAO(nlist_s), list_bas(2,nlist_s),
     &        Index(nIndex)
      Real*8 Dens(nDens,nD), Rho(nRho,mGrid), Fact(mdc**2),
     &       TabAO(nTabAO)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*define _TIME_
#ifdef _TIME_
#endif
#ifdef _DEBUGPRINT_
      If (Debug) Then
         Call RecPrt('Rho_CAS:Dens',' ',Dens,nDens,nD)
         Write (6,*) 'mAO=',mAO
         Write (6,*) 'mGrid=',mGrid
         Write (6,*) 'nTabAO=',nTabAO
         Write (6,*) 'nlist_s=',nlist_s
         Do iList_s = 1, nList_s
            Write (6,*) 'iList_s=',iList_s
            iSkal = list_s(1,ilist_s)
            iCmp  = iSD( 2,iSkal)
            iBas  = iSD( 3,iSkal)
            mTabAO=iBas*iCmp
            Call RecPrt('Rho_CAS: TabAO',' ',
     &                  TabAO(ipTabAO(iList_s)),mAO,mGrid*mTabAO)
         End Do
      End If
#endif
*
      Call FZero(Rho,nRho*mGrid)
*                                                                      *
************************************************************************
*                                                                      *
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         kDCRE=list_s(2,ilist_s)
         index_i=list_bas(2,ilist_s)
         mdci  = iSD(10,iSkal)
         iShell= iSD(11,iSkal)
         nFunc_i=iBas*iCmp
*
         mDij=nFunc_i*nFunc_i
*
*------- Get the Density
*
         ijS=iTri(iShell,iShell)
         ip_Tmp=ipDijs
         Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)
*
         ij = (mdci-1)*mdc + mdci
*
         iER=iEOr(kDCRE,kDCRE)
         lDCRER=NrOpr(iER)
*
         ip_D_a=ipDij+lDCRER*mDij
         ip_D_b=ip_D_a
         If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
*
         If (nD.eq.1) Then
            Call Do_Rho7a_d(Rho,nRho,mGrid,
     &                      DeDe(ip_D_a),mAO,TabAO(ipTabAO(iList_s)),
     &                      iBas,iBas_Eff,iCmp,
     &                      Fact(ij),Index(index_i))
         Else
            Call Do_Rho7_d(Rho,nRho,mGrid,
     &                     DeDe(ip_D_a),DeDe(ip_D_b),mAO,
     &                     TabAO(ipTabAO(iList_s)),
     &                     iBas,iBas_Eff,iCmp,
     &                     Fact(ij),Index(index_i))
         End If

*
         Do jlist_s=1,ilist_s-1
            jSkal = list_s(1,jlist_s)
            kDCRR=list_s(2,jlist_s)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            index_j =list_bas(2,jlist_s)
            mdcj  = iSD(10,jSkal)
            jShell= iSD(11,jSkal)
            nFunc_j=jBas*jCmp
*
            mDij=nFunc_i*nFunc_j
*
*---------- Get the Density
*
            ijS=iTri(iShell,jShell)
            ip_Tmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*)
               Write (6,*) 'iS,jS=',iSkal,jSkal
               Write (6,*) 'mDCRij,mDij=',mDCRij,mDij
               Write (6,*) 'ipDij,ipDSij,ipDDij=',ipDij,ipDSij,ipDDij
            End If
#endif
*
            ij = (mdcj-1)*mdc + mdci
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
               Write (6,*) 'Rho_CAS'
               nBB = iBas*jBas
               nCC = iCmp*jCmp
               Write (6,*) 'iShell,jshell=', iShell,jshell
               Write (6,*) 'kDCRE,kDCRR=', kDCRE,kDCRR
               Write (6,*) 'iER,lDCRER=',iER,lDCRER
               Call RecPrt('DAij',' ',DeDe(ip_D_a),nBB,nCC)
               If (nD.ne.1)
     &            Call RecPrt('DBij',' ',DeDe(ip_D_b),nBB,nCC)
            End If
#endif
*
            If (nD.eq.1) Then
               If (iShell.ge.jShell) Then
               Call Do_Rho7a(Rho,nRho,mGrid,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Two,Index(index_i),Index(index_j))
               Else
               Call Do_Rho7a(Rho,nRho,mGrid,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Two,Index(index_i),Index(index_j))
               End If
            Else
               If (iShell.ge.jShell) Then
               Call Do_Rho7_(Rho,nRho,mGrid,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Two,Index(index_i),Index(index_j))
               Else
               Call Do_Rho7_(Rho,nRho,mGrid,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Two,Index(index_i),Index(index_j))
               End If
            End If
*
         End Do                      ! jlist_s
      End Do                         ! ilist_s
*
#ifdef _DEBUGPRINT_
      If (Debug) Then
c        Do iGrid=1,mGrid_Eff
c           Write (*,*) (Rho(iRho,iGrid),iRho=1,nRho)
c        End Do
         Call RecPrt('Rho_CAS: Rho',' ',Rho,nRho,mGrid)
      End If
*
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Dens)
#endif
#ifdef _TIME_
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nSym)
      End
      Subroutine Do_Rho7a(Rho,nRho,mGrid,
     &                    DAij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid), DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer Index_i(iBas_Eff*iCmp), Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
#endif
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB = Index_j(jCB_Eff)
*

         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
               Prod_55=TabAO1(5,iGrid,iCB_Eff)*TabAO2(5,iGrid,jCB_Eff)
               Prod_88=TabAO1(8,iGrid,iCB_Eff)*TabAO2(8,iGrid,jCB_Eff)
               Prod_10=TabAO1(10,iGrid,iCB_Eff)*TabAO2(10,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(5,iGrid)= Rho(5,iGrid)
     &                     + (Prod_55+Prod_88+Prod_10)*DAij_
            End Do    ! iGrid
*
         End Do          ! iCB
      End Do             ! jCB
*
#ifdef _TIME_
#endif
      Return
      End
      Subroutine Do_Rho7_(Rho,nRho,mGrid,
     &                    DAij,DBij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid),
     &       DAij(iBas*iCmp,jBas*jCmp), DBij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer Index_i(iBas_Eff*iCmp), Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB = Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            DBij_=DBij(iCB,jCB)*Fact
            Dij_ =Half*(Abs(DAij_)+Abs(DBij_))
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
               Prod_55=TabAO1(5,iGrid,iCB_Eff)*TabAO2(5,iGrid,jCB_Eff)
               Prod_88=TabAO1(8,iGrid,iCB_Eff)*TabAO2(8,iGrid,jCB_Eff)
               Prod_10=TabAO1(10,iGrid,iCB_Eff)*TabAO2(10,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) +     Prod_11      *DBij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(5,iGrid)=Rho(5,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(6,iGrid)=Rho(6,iGrid) + (Prod_21+Prod_12)*DBij_
               Rho(7,iGrid)=Rho(7,iGrid) + (Prod_31+Prod_13)*DBij_
               Rho(8,iGrid)=Rho(8,iGrid) + (Prod_41+Prod_14)*DBij_
               Rho( 9,iGrid)= Rho( 9,iGrid)
     &                      + (Prod_55+Prod_88+Prod_10)*DAij_
               Rho(10,iGrid)= Rho(10,iGrid)
     &                      + (Prod_55+Prod_88+Prod_10)*DBij_
            End Do    ! iGrid
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
      Subroutine Do_Rho7a_d(Rho,nRho,mGrid,
     &                    DAii,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                    Fact,Index_i)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid), DAii(iBas*iCmp,iBas*iCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp)
      Integer Index_i(iBas_Eff*iCmp)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
#endif
      Do jCB_Eff = 1, iBas_Eff*iCmp
         jCB=Index_i(jCB_Eff)
*
         DAii_=DAii(jCB,jCB)*Fact
         Do iGrid = 1, mGrid
            Prod_11=TabAO1(1,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_21=TabAO1(2,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_31=TabAO1(3,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_41=TabAO1(4,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_55=TabAO1(5,iGrid,jCB_Eff)*TabAO1(5,iGrid,jCB_Eff)
            Prod_88=TabAO1(8,iGrid,jCB_Eff)*TabAO1(8,iGrid,jCB_Eff)
            Prod_10=TabAO1(10,iGrid,jCB_Eff)*TabAO1(10,iGrid,jCB_Eff)
*
            Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11*DAii_
            Rho(2,iGrid)=Rho(2,iGrid) + Two*Prod_21*DAii_
            Rho(3,iGrid)=Rho(3,iGrid) + Two*Prod_31*DAii_
            Rho(4,iGrid)=Rho(4,iGrid) + Two*Prod_41*DAii_
            Rho(5,iGrid)= Rho(5,iGrid)
     &                  + (Prod_55+Prod_88+Prod_10)*DAii_
         End Do    ! iGrid
*
         Do iCB_Eff = 1, jCB_Eff-1
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAii(iCB,jCB)*Fact*Two
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
               Prod_55=TabAO1(5,iGrid,iCB_Eff)*TabAO1(5,iGrid,jCB_Eff)
               Prod_88=TabAO1(8,iGrid,iCB_Eff)*TabAO1(8,iGrid,jCB_Eff)
               Prod_10=TabAO1(10,iGrid,iCB_Eff)*TabAO1(10,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(5,iGrid)= Rho(5,iGrid)
     &                     + (Prod_55+Prod_88+Prod_10)*DAij_
            End Do    ! iGrid
*
         End Do          ! iCB
      End Do             ! jCB
*
#ifdef _TIME_
#endif
      Return
      End
      Subroutine Do_Rho7_d(Rho,nRho,mGrid,
     &                     DAii,DBii,
     &                     mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                     Fact,Index_i)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid),
     &       DAii(iBas*iCmp,iBas*iCmp), DBii(iBas*iCmp,iBas*iCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp)
      Integer Index_i(iBas_Eff*iCmp)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, iBas_Eff*iCmp
         jCB=Index_i(jCB_Eff)
*
         DAii_=DAii(jCB,jCB)*Fact
         DBii_=DBii(jCB,jCB)*Fact
         Dii_ =Half*(Abs(DAii_)+Abs(DBii_))
         Do iGrid = 1, mGrid
            Prod_11=TabAO1(1,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_21=TabAO1(2,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_31=TabAO1(3,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_41=TabAO1(4,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
            Prod_55=TabAO1(5,iGrid,jCB_Eff)*TabAO1(5,iGrid,jCB_Eff)
            Prod_88=TabAO1(8,iGrid,jCB_Eff)*TabAO1(8,iGrid,jCB_Eff)
            Prod_10=TabAO1(10,iGrid,jCB_Eff)*TabAO1(10,iGrid,jCB_Eff)
*
            Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11*DAii_
            Rho(2,iGrid)=Rho(2,iGrid) +     Prod_11*DBii_
            Rho(3,iGrid)=Rho(3,iGrid) + Two*Prod_21*DAii_
            Rho(4,iGrid)=Rho(4,iGrid) + Two*Prod_31*DAii_
            Rho(5,iGrid)=Rho(5,iGrid) + Two*Prod_41*DAii_
            Rho(6,iGrid)=Rho(6,iGrid) + Two*Prod_21*DBii_
            Rho(7,iGrid)=Rho(7,iGrid) + Two*Prod_31*DBii_
            Rho(8,iGrid)=Rho(8,iGrid) + Two*Prod_41*DBii_
            Rho( 9,iGrid)= Rho( 9,iGrid)
     &                   + (Prod_55+Prod_88+Prod_10)*DAii_
            Rho(10,iGrid)= Rho(10,iGrid)
     &                   + (Prod_55+Prod_88+Prod_10)*DBii_
         End Do    ! iGrid
*
         Do iCB_Eff = 1, jCB_Eff-1
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAii(iCB,jCB)*Fact*Two
            DBij_=DBii(iCB,jCB)*Fact*Two
            Dij_ =Half*(Abs(DAij_)+Abs(DBij_))
*
            Do iGrid = 1, mGrid
               Prod_11= TabAO1(1,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_21= TabAO1(2,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_12= TabAO1(1,iGrid,iCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_31= TabAO1(3,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_13= TabAO1(1,iGrid,iCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_41= TabAO1(4,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_14= TabAO1(1,iGrid,iCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
               Prod_55=TabAO1(5,iGrid,iCB_Eff)*TabAO1(5,iGrid,jCB_Eff)
               Prod_88=TabAO1(8,iGrid,iCB_Eff)*TabAO1(8,iGrid,jCB_Eff)
               Prod_10=TabAO1(10,iGrid,iCB_Eff)*TabAO1(10,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) +     Prod_11      *DBij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(5,iGrid)=Rho(5,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(6,iGrid)=Rho(6,iGrid) + (Prod_21+Prod_12)*DBij_
               Rho(7,iGrid)=Rho(7,iGrid) + (Prod_31+Prod_13)*DBij_
               Rho(8,iGrid)=Rho(8,iGrid) + (Prod_41+Prod_14)*DBij_
               Rho( 9,iGrid)= Rho( 9,iGrid)
     &                      + (Prod_55+Prod_88+Prod_10)*DAij_
               Rho(10,iGrid)= Rho(10,iGrid)
     &                      + (Prod_55+Prod_88+Prod_10)*DBij_
            End Do    ! iGrid
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
