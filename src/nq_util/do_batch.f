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
* Copyright (C) 2000,2021, Roland Lindh                                *
*               2021, Jie Bao                                          *
************************************************************************
      Subroutine Do_Batch(Kernel,Func,mGrid,
     &                    list_s,nlist_s,List_Exp,List_Bas,
     &                    Index,nIndex,
     &                    FckInt,nFckDim,nFckInt,
     &                    ipTabAO,mAO,nSym,nD,
     &                    nP2_ontop,Do_Mo,
     &                    TabMO,TabSO,nMOs,
     &                    Do_Grad,Grad,nGrad,ndRho_dR,nGrad_Eff,iNQ,
     &                    EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use iSD_data
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Phase_Info
      use KSDFT_Info
      use nq_Grid, only: Grid, Weights, Rho, GradRho, Sigma, nRho
      use nq_Grid, only: vRho, vSigma, vTau, vLapl
      use nq_Grid, only: l_CASDFT, TabAO, TabAO_Pack, dRho_dR
      use nq_Grid, only: F_xc, F_xca, F_xcb
      use nq_Grid, only: Fact, Tmp, SOs, Angular, Mem
      use nq_Grid, only: D1UnZip, P2UnZip
      use nq_pdft
      use nq_MO, only: DoIt, CMO, D1MO, P2_ontop
      use Grid_On_Disk
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "SysDef.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "debug.fh"
#include "ksdft.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
#include "pamint.fh"
      Integer list_s(2,nlist_s),List_Exp(nlist_s),
     &        ipTabAO(nlist_s+1,2),Index(nIndex),
     &        List_Bas(2,nlist_s)
      Real*8 A(3), RA(3), Grad(nGrad), FckInt(nFckInt,nFckDim),
     &       TabMO(mAO,mGrid,nMOs),TabSO(mAO,mGrid,nMOs),
     &       PDFTPot1(nPot1),PDFTFocI(nPot1),PDFTFocA(nPot1)
      Logical Do_Grad,Do_Mo,Unpack
      Logical l_tanhr
      Real*8 P2_ontop_d(nP2_ontop,nGrad_Eff,mGrid)
      Real*8,DIMENSION(:),ALLOCATABLE::P2MOCube,P2MOCubex,P2MOCubey,
     &                                 P2MOCubez,MOs,MOx,MOy,MOz
*     MOs,MOx,MOy and MOz are for active MOs.
*     MOas is for all MOs.
      Integer nPMO3p
      Real*8 EG_OT(nTmpPUVX)
      Real*8, Allocatable:: RhoI(:,:), RhoA(:,:)
      Real*8, Allocatable:: TmpCMO(:)
      Integer, Allocatable:: TDoIt(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      nTabAO=Size(TabAO)
      nCMO  =Size(CMO)
      T_Rho=T_X*1.0D-4
      l_tanhr=.false.

      CALL PDFTMemAlloc(mGrid,nOrbt)
      If ( Functional_Type.eq.CASDFT_Type .or.
     &     l_casdft ) Then
         mRho = nP2_ontop
         Call mma_allocate(RhoI,mRho,mGrid,Label='RhoI')
         Call mma_allocate(RhoA,mRho,mGrid,Label='RhoA')
         RhoI(:,:)=Zero
         RhoA(:,:)=Zero
      Else
         mRho=-1
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the AOs on the grid points.                             *
*                                                                      *
************************************************************************
*
      TabAO(:,:,:)=Zero
      UnPack=.False.
      If (NQ_Direct.eq.Off .and. (Grid_Status.eq.Use_Old .and.
     &      .Not.Do_Grad         .and.
     &    Functional_Type.eq.Old_Functional_Type)) Then
*
*------- Retrieve the AOs from disc
*
         Call iDaFile(Lu_Grid,2,ipTabAO,2*(nlist_s+1),iDisk_Grid)
         mTabAO=ipTabAO(nlist_s+1,2)-1
         Call dDaFile(Lu_Grid,2,TabAO,mTabAO,iDisk_Grid)
         Unpack=Packing.eq.On
*
      Else
*
*------- Generate the values of the AOs on the grid
*
         Mem(:)=Zero
         ipxyz=1
*
         iOff = 1
         Do ilist_s=1,nlist_s
            ish=list_s(1,ilist_s)

            iShll = iSD( 0,iSh)
            iAng  = iSD( 1,iSh)
            iCmp  = iSD( 2,iSh)
            iBas  = iSD( 3,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            iPrim = iSD( 5,iSh)
            iPrim_Eff=List_Exp(ilist_s)
            iAO   = iSD( 7,iSh)
            mdci  = iSD(10,iSh)
            iShll = iSD(0,iSh)
            iCnttp= iSD(13,iSh)
            iCnt  = iSD(14,iSh)
            A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
            nDrv     = mRad - 1
            nForm    = 0
            Do iDrv  = 0, nDrv
               nForm = nForm + nElem(iDrv)
            End Do
            nTerm    = 2**nDrv
            nxyz     = mGrid*3*(iAng+mRad)
!           nRadial  = iBas_Eff*mGrid*mRad
            ipRadial = ipxyz + nxyz
*
            iR=list_s(2,ilist_s)
*
            ipx=iPhase(1,iR)
            ipy=iPhase(2,iR)
            ipz=iPhase(3,iR)
            px=DBLE(iPhase(1,iR))
            py=DBLE(iPhase(2,iR))
            pz=DBLE(iPhase(3,iR))
            RA(1) = px*A(1)
            RA(2) = py*A(2)
            RA(3) = pz*A(3)
            iSym=NrOpr(iR)
*
*---------- Evaluate AOs at RA
*
            ipTabAO(iList_s,1)=iOff
*                                                                      *
            Call AOEval(iAng,mGrid,Grid,Mem(ipxyz),RA,
     &                  Shells(iShll)%Transf,
     &                  RSph(ipSph(iAng)),nElem(iAng),iCmp,
     &                  Angular,nTerm,nForm,T_X,mRad,
     &                  iPrim,iPrim_Eff,Shells(iShll)%Exp,
     &                  Mem(ipRadial),iBas_Eff,
     &                  Shells(iShll)%pCff(1,iBas-iBas_Eff+1),
     &                  TabAO_Pack(iOff:),
     &                  mAO,px,py,pz,ipx,ipy,ipz)
            iOff = iOff + mAO*mGrid*iBas_Eff*iCmp
*
         End Do
         ipTabAO(nList_s+1,1)=iOff
*
*        AOs are packed and written to disk.
*
         If (NQ_Direct.eq.Off .and. (Grid_Status.eq.Regenerate .and.
     &       .Not.Do_Grad)) Then
*
            If (Packing.eq.On) Then
               Unpack=.True.
*
*------------- Pack before they are put on disc
*
               jOff = 1
               Do ilist_s=1,nlist_s
                  ish=list_s(1,ilist_s)
                  iCmp  = iSD( 2,iSh)
                  iBas_Eff = List_Bas(1,ilist_s)
                  nData=mAO*mGrid*iBas_Eff*iCmp
*
*                 Check if we should store any AOs at all!
*
                  iOff = ipTabAO(ilist_s,1)
                  If (nData.gt.SIZE(Tmp)) Then
                     Call WarningMessage(2,'nData.gt.SIZE(Tmp)')
                     Call Abend()
                  End If
                  call dcopy_(nData,TabAO_Pack(iOff:),1,Tmp,1)
                  Call PkR8(0,nData,nByte,Tmp,TabAO_Pack(jOff:))
                  mData = (nByte+RtoB-1)/RtoB
                  If (mData.gt.nData) Then
                     Call WarningMessage(2,'mData.gt.nData')
                     Write (6,*) 'nData=',nData
                     Write (6,*) 'nData=',nData
                     Call Abend()
                  End If
                  ipTabAO(iList_s,2)=nByte
                  jOff = jOff + mData
               End Do
               ipTabAO(nList_s+1,2)=jOff
            Else
               ipTabAO(nList_s+1,2)=ipTabAO(nList_s+1,1)
            End If
*
            Call iDaFile(Lu_Grid,1,ipTabAO,2*(nlist_s+1),iDisk_Grid)
            mTabAO=ipTabAO(nList_s+1,2)-1
            Call dDaFile(Lu_Grid,1,TabAO,mTabAO,iDisk_Grid)
*
         End If
*
      End If
*
*---- Unpack AOs
*
      If (Unpack) Then
*
         jOff = ipTabAO(nlist_s+1,2)
         Do ilist_s=nlist_s,1,-1
            ish=list_s(1,ilist_s)
            iCmp  = iSD( 2,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            nData=mAO*mGrid*iBas_Eff*iCmp
            nByte=ipTabAO(ilist_s,2)
*
            iOff = ipTabAO(ilist_s,1)
            If (nByte.gt.0) Then
               mData = (nByte+RtoB-1)/RtoB
               jOff = jOff - mData
               If (mData.gt.SIZE(Tmp)) Then
                  Call WarningMessage(2,'mData.gt.SIZE(Tmp)')
                  Call Abend()
               End If
               Call UpkR8(0,nData,nByte,TabAO_Pack(jOff:),Tmp)
               call dcopy_(nData,Tmp,1,TabAO_Pack(iOff:),1)
            Else
               mData=0
               TabAO_Pack(1:nData)=Zero
            End If
         End Do
*
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Evaluate some MOs on the grid                                    *
*                                                                      *
************************************************************************
*                                                                      *
      If (Do_MO) Then
         Call FZero(TabMO,mAO*mGrid*nMOs)
         Call FZero(TabSO,mAO*mGrid*nMOs)
*
         Call mma_Allocate(TmpCMO,nCMO,Label='TmpCMO')
         Call mma_Allocate(TDoIt,nMOs,Label='TDoIt')
         Do ilist_s=1,nlist_s
            ish=list_s(1,ilist_s)
            iCmp  = iSD( 2,iSh)
            iBas  = iSD( 3,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            iPrim = iSD( 5,iSh)
            iAO   = iSD( 7,iSh)
            mdci  = iSD(10,iSh)
*
*---------- Allocate memory for SO and MO
*
            kAO   = iCmp*iBas*mGrid
            nDeg  = nSym/dc(mdci)%nStab
            nSO   = kAO*nDeg*mAO
            Call FZero(SOs,nSO)
*
            iR=list_s(2,ilist_s)
            iSym=NrOpr(iR)
*
*---------- Distribute contributions of AOs if this particular shell
*           on to the SOs of this shell. The SOs are only stored
*           temporarily!
*
            Call SOAdpt_NQ(TabAO_Pack(ipTabAO(iList_s,1):),mAO,mGrid,
     &                     iBas,iBas_Eff,iCmp,iSym,SOs,nDeg,iAO)
*
            Call  SODist2(SOs,mAO,mGrid,iBas,iCmp,nDeg,TabSO,
     &                    nMOs,iAO,TmpCMO,nCMO,TDoIt)
*
            Call  SODist(SOs,mAO,mGrid,iBas,iCmp,nDeg,TabMO,
     &                   nMOs,iAO,CMO,nCMO,DoIt)
*
         End Do
         Call mma_deAllocate(TDoIt)
         Call mma_deAllocate(TmpCMO)
      End If
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*---- Compute Rho, Grad Rho, Tau, Laplacian, and the Sigma vectors.
*     In case of gradient calculations compute Cartesian derivatives
*     of Rho, Grad Rho, Tau, and the Laplacian.
*                                                                      *
      Call Mk_Rho(list_s,nlist_s,Fact,ndc,list_bas,Index,nIndex,Do_Grad)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      If (l_casdft) then
         Dens_t1=Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
         Dens_a1=Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
         Dens_b1=Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)

         nPMO3p=1
         IF (lft.and.lGGA) nPMO3p=mGrid*NASHT

         CALL mma_allocate(P2MOCube,mGrid*NASHT)
         CALL mma_allocate(P2MOCubex,nPMO3p)
         CALL mma_allocate(P2MOCubey,nPMO3p)
         CALL mma_allocate(P2MOCubez,nPMO3p)
         CALL mma_allocate(MOs,mGrid*NASHT)
         CALL mma_allocate(MOx,mGrid*NASHT)
         CALL mma_allocate(MOy,mGrid*NASHT)
         CALL mma_allocate(MOz,mGrid*NASHT)

         CALL CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,
     &                     nPMO3p,
     &                     MOs,MOx,MOy,MOz,TabMO,P2Unzip,mAO,mGrid,nMOs,
     &                     Do_Grad)
         Call Fzero(P2_ontop,nP2_ontop*mGrid)

         If (.not.Do_Grad) then !regular MO-based run
            Call Do_PI2(D1MO,SIZE(D1MO),TabMO,mAO,mGrid,
     &                  nMOs,P2_ontop,nP2_ontop,RhoI,
     &                  RhoA,mRho,Do_Grad,
     &                  P2MOCube,MOs,MOx,MOy,MOz)
         Else !AO-based run for gradients
!           nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
            P2_ontop_d(:,:,:) = 0
            Call  Do_Pi2grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &                       P2_ontop,nP2_ontop,nGrad_Eff,
     &                       list_s,nlist_s,list_bas,
     &                       D1MO,SIZE(D1MO),TabMO,P2_ontop_d,
     &                       RhoI,RhoA,mRho,nMOs,CMO,
     &                       nCMO,TabSO,nsym,lft,
     &                       P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,
     &                       nPMO3p,MOs,MOx,MOy,MOz)
         End If

         CALL TranslateDens(P2_OnTop,dRho_dr,P2_OnTop_d,
     &                       l_tanhr,nRho,mGrid,nP2_OnTop,
     &                       ndRho_dR,nGrad_Eff,Do_Grad)

         CALL mma_deallocate(P2MOCube)
         CALL mma_deallocate(P2MOCubex)
         CALL mma_deallocate(P2MOCubey)
         CALL mma_deallocate(P2MOCubez)
         CALL mma_deallocate(MOs)
         CALL mma_deallocate(MOx)
         CALL mma_deallocate(MOy)
         CALL mma_deallocate(MOz)

         If (lGGA) Then
            If (nD.eq.1) Then
               Do iGrid=1, mGrid
                  Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                          +GradRho(2,iGrid)**2
     &                          +GradRho(3,iGrid)**2
               End Do
            Else
               Do iGrid=1, mGrid
                  Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                          +GradRho(2,iGrid)**2
     &                          +GradRho(3,iGrid)**2
                  Sigma(2,iGrid)=GradRho(1,iGrid)*GradRho(4,iGrid)
     &                          +GradRho(2,iGrid)*GradRho(5,iGrid)
     &                          +GradRho(3,iGrid)*GradRho(6,iGrid)
                  Sigma(3,iGrid)=GradRho(4,iGrid)**2
     &                          +GradRho(5,iGrid)**2
     &                          +GradRho(6,iGrid)**2
               End Do
            End If
         End If

*        Integrate out the number of electrons
         Dens_t2=Dens_t2+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
         Dens_a2=Dens_a2+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
         Dens_b2=Dens_b2+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)

      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Integrate out the number of electrons, |grad|, and tau
*
      If (Functional_type.eq.LDA_Type) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
      Else If (Functional_type.eq.GGA_type) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD,T_Rho)
      Else If (Functional_type.eq.meta_GGA_type1) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD,T_Rho)
         Tau_I =Tau_I +Compute_Tau (Weights,mGrid,nD,T_Rho)
      Else If (Functional_type.eq.meta_GGA_type2) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD,T_Rho)
         Tau_I =Tau_I +Compute_Tau (Weights,mGrid,nD,T_Rho)
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- Evaluate the functional on the grid                              *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      vRho(:,1:mGrid)=Zero
      If (Allocated(vSigma)) vSigma(:,1:mGrid)=Zero
      If (Allocated(vTau)) vTau(:,1:mGrid)=Zero
      If (Allocated(vLapl)) vLapl(:,1:mGrid)=Zero
      F_xc(1:mGrid)=Zero
      If (l_casdft) Then
         F_xca(1:mGrid)=Zero
         F_xcb(1:mGrid)=Zero
      End If
*
*1)   evaluate the energy density, the derivative of the functional with
*     respect to rho and grad rho.
*
      Call Kernel(mGrid,nD)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
*     Integrate the energy of the functional
*
      Func=Func+DDot_(mGrid,Weights,1,F_xc,1)
      If (l_casdft) Then
         Funcaa=Funcaa+DDot_(mGrid,Weights,1,F_xca,1)
         Funcbb=Funcbb+DDot_(mGrid,Weights,1,F_xcb,1)
         Funccc=Func-Funcaa-Funcbb
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Do_Grad) Then
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the DFT contribution to the gradient                  *
*                                                                      *
************************************************************************
*                                                                      *
         Call DFT_Grad(Grad,nGrad,nD,Grid,mGrid,dRho_dR,ndRho_dR,
     &                nGrad_Eff,Weights,iNQ)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         If (l_casdft) Then
*                                                                      *
************************************************************************
*                                                                      *
*------- For MC-PDFT optionally compute stuff for the CP-MC-PDFT       *
*                                                                      *
************************************************************************
*                                                                      *
           If (do_pdftPot) then
              CALL mma_allocate(MOs ,mGrid*NASHT)
              CALL TransferMO(MOas,TabMO,mAO,mGrid,nMOs,1)
              IF (lft.and.lGGA) THEN
                 CALL TransferMO(MOax,TabMO,mAO,mGrid,nMOs,2)
                 CALL TransferMO(MOay,TabMO,mAO,mGrid,nMOs,3)
                 CALL TransferMO(MOaz,TabMO,mAO,mGrid,nMOs,4)
              END IF
              CALL TransActMO(MOs, TabMO,mAO,mGrid,nMOs)
              Call Calc_Pot1(PDFTPot1,TabMO,mAO,mGrid,nMOs,P2_ontop,
     &                       nP2_ontop,MOas)
              Call Calc_Pot2(EG_OT,mGrid,P2_ontop,nP2_ontop)
              Call PDFTFock(PDFTFocI,PDFTFocA,D1Unzip,mGrid,MOs)
              CALL mma_deallocate(MOs)
           End If
*                                                                      *
************************************************************************
*                                                                      *
         Else
*                                                                      *
************************************************************************
*                                                                      *
*------- Compute the DFT contribution to the Fock matrix               *
*                                                                      *
************************************************************************
*                                                                      *
           Call DFT_Int(list_s,nlist_s,FckInt,nFckInt,nD,Fact,ndc,
     &                  list_bas)
*                                                                      *
************************************************************************
*                                                                      *
         End If
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      CALL PDFTMemDeAlloc()

      If (Allocated(RhoI)) Then
         Call mma_deallocate(RhoI)
         Call mma_deallocate(RhoA)
      End If

      Return
      End
