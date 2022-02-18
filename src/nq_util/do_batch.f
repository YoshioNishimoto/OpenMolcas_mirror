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
     &                    Index,nIndex,FckInt,nFckDim,nFckInt,
     &                    mAO,nD,nP2_ontop,Do_Mo,TabMO,TabSO,nMOs,
     &                    Do_Grad,Grad,nGrad,ndRho_dR,nGrad_Eff,iNQ,
     &                    EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use iSD_data
      use SOAO_Info, only: iAOtSO
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Phase_Info
      use KSDFT_Info
      use nq_Grid, only: Grid, Weights, Rho, nRho
      use nq_Grid, only: GradRho, Sigma
      use nq_Grid, only: l_CASDFT, TabAO, TabAO_Pack, dRho_dR
      use nq_Grid, only: F_xc, F_xca, F_xcb, kAO, Grid_AO
      use nq_Grid, only: Fact, Angular, Mem
      use nq_Grid, only: D1UnZip, P2UnZip
      use nq_Grid, only: Dens_AO, iBfn_Index
      use nq_pdft
      use nq_MO, only: CMO, D1MO, P2_ontop
      use Grid_On_Disk
      use nq_Info
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "SysDef.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "debug.fh"
#include "ksdft.fh"
#include "nsd.fh"
#include "setup.fh"
#include "pamint.fh"
      Integer list_s(2,nlist_s),List_Exp(nlist_s),Index(nIndex),
     &        List_Bas(2,nlist_s)
      Real*8 A(3), RA(3), Grad(nGrad), FckInt(nFckInt,nFckDim),
     &       TabMO(mAO,mGrid,nMOs),TabSO(mAO,mGrid,nMOs),
     &       PDFTPot1(nPot1),PDFTFocI(nPot1),PDFTFocA(nPot1)
      Logical Do_Grad,Do_Mo
      Logical l_tanhr
      Real*8 P2_ontop_d(nP2_ontop,nGrad_Eff,mGrid)
      Real*8,DIMENSION(:),ALLOCATABLE::P2MOCube,P2MOCubex,P2MOCubey,
     &                                 P2MOCubez,MOs,MOx,MOy,MOz
*     MOs,MOx,MOy and MOz are for active MOs.
*     MOas is for all MOs.
      Integer nPMO3p
      Real*8 EG_OT(nTmpPUVX)
      Real*8, Allocatable:: RhoI(:,:), RhoA(:,:)
      Real*8, Allocatable:: TabAO_Tmp(:)
      Integer :: TabAO_Size(2)
      Integer, Allocatable :: Tmp_Index(:,:)

*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      nCMO  =Size(CMO)
      l_tanhr=.false.

      If (l_casdft ) Then
         CALL PDFTMemAlloc(mGrid,nOrbt)
         mRho = nP2_ontop
         Call mma_allocate(RhoI,mRho,mGrid,Label='RhoI')
         Call mma_allocate(RhoA,mRho,mGrid,Label='RhoA')
         RhoI(:,:)=Zero
         RhoA(:,:)=Zero
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Set up an indexation translation between the running index of
*     the AOIntegrals and the actual basis function index
*
      nBfn=0
      Do iList_s = 1, nList_s
         iBas_Eff=List_Bas(1,ilist_s)
         iSkal    =list_s(1,ilist_s)
         iCmp  = iSD( 2,iSkal)
         nBfn=nBfn+iBas_Eff*iCmp
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the AOs on the grid points.                             *
*                                                                      *
************************************************************************
*
      TabAO(:,:,:)=Zero
      TabAO_Size(:)=0
*                                                                      *
************************************************************************
*                                                                      *
!     Compute AO's or retrive from disk
!
      If (NQ_Direct.eq.Off .and. (Grid_Status.eq.Use_Old .and.
     &      .Not.Do_Grad         .and.
     &    Functional_Type.eq.Old_Functional_Type)) Then
*                                                                      *
************************************************************************
*                                                                      *
*------- Retrieve (and unpack) the AOs from disc
*
         Call iDaFile(Lu_Grid,2,TabAO_Size,2,iDisk_Grid)
         If (TabAO_Size(1)==0) Then
            Call Terminate()
            Return
         End If
         nBfn=TabAO_Size(1)
         Call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
         Call iDaFile(Lu_Grid,2,iBfn_Index,Size(iBfn_Index),iDisk_Grid)

         nByte=TabAO_Size(2)
         If (Packing.eq.On) Then
            mTabAO = (nByte+RtoB-1)/RtoB
         Else
            mTabAO = nByte
         End If
         Call dDaFile(Lu_Grid,2,TabAO,mTabAO,iDisk_Grid)
*
         If (Packing.eq.On) Then
            nData=Size(TabAO)
            Call mma_Allocate(TabAO_Tmp,nData,Label='TabAO_Tmp')
            nByte=TabAO_Size(2)
            Call UpkR8(0,nData,nByte,TabAO_Pack,TabAO_Tmp)
            TabAO_Pack(:)=TabAO_Tmp(:)
            Call mma_deAllocate(TabAO_Tmp)
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
!        Compute the AO's
!
         Call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
         iBfn_Index(:,:)=0
*
*------- Generate the values of the AOs on the grid
*
         Mem(:)=Zero
         ipxyz=1
*
*#define _ANALYSIS_
#ifdef _ANALYSIS_
      Thr=T_Y
      Write (6,*)
      Write (6,*) ' Sparsity analysis of AO blocks'
      mlist_s=0
#endif
!        iOff  = 1
         iBfn  = 0
         iBfn_s= 0
         iBfn_e= 0
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
!
!           Set up the unsifted version of iBfn_Index
!
            iAdd = iBas-iBas_Eff
            iBfn_s = iBfn + 1
            Do i1 = 1, iCmp
               iSO1 = iAOtSO(iAO+i1,0) ! just used when nIrrep=1
               Do i2 = 1, iBas_Eff
                  IndAO1 = i2 + iAdd
                  Indi = iSO1 + IndAO1 -1

                  iBfn = iBfn + 1
                  iBfn_Index(1,iBfn) = Indi
                  iBfn_Index(2,iBfn) = ilist_s
                  iBfn_Index(3,iBfn) = i1
                  iBfn_Index(4,iBfn) = i2
                  iBfn_Index(5,iBfn) = mdci
                  iBfn_Index(6,iBfn) = IndAO1
               End Do
            End Do
            iBfn_e = iBfn

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
*
*---------- Evaluate AOs at RA
*
            Call AOEval(iAng,mGrid,Grid,Mem(ipxyz),RA,
     &                  Shells(iShll)%Transf,
     &                  RSph(ipSph(iAng)),nElem(iAng),iCmp,
     &                  Angular,nTerm,nForm,T_Y,mRad,
     &                  iPrim,iPrim_Eff,Shells(iShll)%Exp,
     &                  Mem(ipRadial),iBas_Eff,
     &                  Shells(iShll)%pCff(1,iBas-iBas_Eff+1),
     &                  TabAO(:,:,iBfn_s:),
     &                  mAO,px,py,pz,ipx,ipy,ipz)
#ifdef _ANALYSIS_
            ix = iDAMax_(mAO*mGrid*iBas_Eff*iCmp,TabAO_Pack(iOff),1)
            TMax = Abs(TabAO_Pack(iOff-1+ix))
            If (TMax<Thr) Then
               mlist_s=mlist_s+1
               Write (6,*) ' ilist_s: ',ilist_s
               Write (6,*) ' TMax:    ',TMax
            End If
#endif
!
!           At this time eliminate individual basis functions which have
!           an insignificant contribution to any of the grid points we
!           are processing at this stage.
!
            Thr=T_Y
            iSkip=0
            kBfn = iBfn_s - 1
            Do jBfn = iBfn_s, iBfn_e
               Call Spectre(SMax)
               If (SMax<Thr) Then
                 iSkip=iSkip+1
               Else
                 kBfn = kBfn + 1
                 If (kBfn/=jBfn) Then
                    TabAO(:,:,kBfn)=TabAO(:,:,jBfn)
                    iBfn_Index(:,kBfn) = iBfn_Index(:,jBfn)
                 End If
               End If
            End Do
            iBfn = kBfn

!           iOff = iBfn*mAO*mGrid + 1
*
*                                                                      *
************************************************************************
*                                                                      *
         End Do
*                                                                      *
************************************************************************
*                                                                      *
         ! reduced the size of the table to be exactly that of the
         ! number of functions that have non-zero contributions.
         If (iBfn/=nBfn) Then
            If (iBfn==0) Then
               TabAO_Size(:)=0
               Call iDaFile(Lu_Grid,1,TabAO_Size,2,iDisk_Grid)
               Call Terminate()
               Return
             End If
            Call mma_allocate(Tmp_Index,6,iBfn,Label='Tmp_Index')
            Tmp_Index(:,1:iBfn) = iBfn_Index(:,1:iBfn)
            Call mma_deallocate(iBfn_Index)
            nBfn=iBfn
            Call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
            iBfn_Index(:,:) = Tmp_Index(:,:)
            Call mma_deallocate(Tmp_Index)
         End If
         TabAO_Size(1)=nBfn

#ifdef _ANALYSIS_
         Write (6,*) ' % AO blocks that can be eliminated: ',
     &             1.0D2*DBLE(mlist_s)/DBLE(nlist_s)
         Write (6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Allocate(Dens_AO,nBfn,nBfn,nD,Label='Dens_AO')
      Call mma_Allocate(Grid_AO,kAO,mGrid,nBfn,nD,Label='Grid_AO')
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Evaluate some MOs on the grid                                    *
*                                                                      *
************************************************************************
*                                                                      *
      If (Do_MO) Then

         ! First, symmatry adapt the AOs
         TabSO(:,:,:)=Zero
         jlist_s=0
         Call mk_SOs(TabSO,mAO,mGrid,nMOs,List_s,List_Bas,nList_s,
     &               jlist_s)

         ! Second, transform SOs to MOs
         TabMO(:,:,:)=Zero
         Call mk_MOs(TabSO,mAO,mGrid,TabMO,nMOs,CMO,nCMO)
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
         Dens_t1=Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,0)
         Dens_a1=Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,1)
         Dens_b1=Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,2)

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
            Call  Do_Pi2grad(mAO,mGrid,P2_ontop,nP2_ontop,nGrad_Eff,
     &                       list_s,nlist_s,list_bas,
     &                       D1MO,SIZE(D1MO),TabMO,P2_ontop_d,
     &                       RhoI,RhoA,mRho,nMOs,CMO,
     &                       nCMO,TabSO,lft,
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

*        Integrate out the number of electrons
         Dens_t2=Dens_t2+Comp_d(Weights,mGrid,Rho,nRho,nD,0)
         Dens_a2=Dens_a2+Comp_d(Weights,mGrid,Rho,nRho,nD,1)
         Dens_b2=Dens_b2+Comp_d(Weights,mGrid,Rho,nRho,nD,2)

      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Allocated(Sigma)) Then
         If (Size(Sigma,1)==1) Then
            Do iGrid=1, mGrid
               Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                       +GradRho(2,iGrid)**2
     &                       +GradRho(3,iGrid)**2
            End Do
         Else
            Do iGrid=1, mGrid
               Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                       +GradRho(2,iGrid)**2
     &                       +GradRho(3,iGrid)**2
               Sigma(2,iGrid)=GradRho(1,iGrid)*GradRho(4,iGrid)
     &                       +GradRho(2,iGrid)*GradRho(5,iGrid)
     &                       +GradRho(3,iGrid)*GradRho(6,iGrid)
               Sigma(3,iGrid)=GradRho(4,iGrid)**2
     &                       +GradRho(5,iGrid)**2
     &                       +GradRho(6,iGrid)**2
            End Do
         End If
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Integrate out the number of electrons, |grad|, and tau
*
      Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD)
      Select Case (Functional_type)
      Case (LDA_Type)
      Case (GGA_type)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD)
      Case (meta_GGA_type1,meta_GGA_type2)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD)
         Tau_I =Tau_I +Compute_Tau (Weights,mGrid,nD)
      End Select
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- Evaluate the functional on the grid                              *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     evaluate the energy density, the derivative of the functional with
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
           Call DFT_Int(list_s,nlist_s,FckInt,nFckInt,nD,Fact,ndc)
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
*     AOs on the grid are (packed and) written to disk.
*
      If (NQ_Direct.eq.Off .and. (Grid_Status.eq.Regenerate .and.
     &    .Not.Do_Grad)) Then
*
         TabAO_Size(1)=nBfn
         If (Packing.eq.On) Then
*
*---------- Pack before they are put on disc
*
            nData=mAO*mGrid*nBfn
            Call mma_Allocate(TabAO_Tmp,nData,Label='TabAO_Tmp')
            TabAO_Tmp(1:nData)=TabAO_Pack(1:nData)
            Call PkR8(0,nData,nByte,TabAO_Tmp,TabAO_Pack)
            mData = (nByte+RtoB-1)/RtoB
            If (mData.gt.nData) Then
               Call WarningMessage(2,'mData.gt.nData')
               Write (6,*) 'nData=',nData
               Write (6,*) 'nData=',nData
               Call Abend()
            End If
            TabAO_Size(2)=nByte
            Call mma_deAllocate(TabAO_Tmp)
         Else
            mData=mAO*mGrid*nBfn
            TabAO_Size(2)=mData
         End If
*
         Call iDaFile(Lu_Grid,1,TabAO_Size,2,iDisk_Grid)
         Call iDaFile(Lu_Grid,1,iBfn_Index,Size(iBfn_Index),iDisk_Grid)
         mTabAO=mData
         Call dDaFile(Lu_Grid,1,TabAO,mTabAO,iDisk_Grid)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Terminate()

      Return

      Contains
        Subroutine Terminate()
           If (l_casdft) CALL PDFTMemDeAlloc()

           If (Allocated(RhoI)) Then
              Call mma_deallocate(RhoI)
              Call mma_deallocate(RhoA)
           End If
           If (Allocated(iBfn_Index)) Call mma_deAllocate(iBfn_Index)
           If (Allocated(Grid_AO)) Call mma_deAllocate(Grid_AO)
           If (Allocated(Dens_AO)) Call mma_deAllocate(Dens_AO)
        End Subroutine Terminate

        Subroutine Spectre(SMax)
           Integer iGrid, iAO
           Real*8 SMax
           SMax=Zero
           Do iGrid = 1, mGrid
             Do iAO = 1, mAO
                SMax = Max (SMax,
     &                 Abs(Weights(iGrid)*TabAO(iAO,iGrid,jBfn)))
             End Do
           End Do
!          If (SMax<Thr) Write (6,*) SMax, TMax
        End Subroutine Spectre

      End Subroutine Do_Batch
