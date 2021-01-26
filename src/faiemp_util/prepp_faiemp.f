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
* Copyright (C) Ben Swerts                                             *
************************************************************************
      SubRoutine PrepP_FAIEMP(nBas_Valence, nBT, nBVT)
************************************************************************
*                                                                      *
* Object: to set up the handling of the 2nd order density matrix for   *
*         the calculation of the 2-electron FAIEMP derivatives         *
*                                                                      *
*     Author: Ben Swerts                                               *
*                                                                      *
* Based on PrepP                                                       *
************************************************************************
      use aces_stuff, only: Gamma_On
      use pso_stuff
      use Basis_Info
      use Sizes_of_Seward, only: S
      use Symmetry_Info, only: nIrrep
      Implicit None
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "etwas.fh"
#include "nac.fh"
      Integer nBas_Valence(0:7),nBT,nBVT,nFro(0:7)
      Character*8 RlxLbl,Method, KSDFT*16
      Logical lPrint
      Integer i,iBas,iGo,iIrrep,ij,ipt,ipTmp1
      Integer iSpin,jBas,nAct,nDens_Valence,nsa,nTst
      Integer iRout,iPrint,iComp
      Real*8  Get_ExFac,CoefX,CoefR
      External Get_ExFac
      Character*60 Fmt
      Real*8, Allocatable:: D1AV(:)

*
*...  Prologue
      iRout = 205
      iPrint = nPrint(iRout)
      lPrint=.True.
      iD0Lbl=1
      iComp=1
*
      nDens = nBT
      nDens_Valence = nBVT
*
      lsa=.False.
      Gamma_On=.False.
      lPSO=.false.
*
*...  Get the method label
      Call Get_cArray('Relax Method',Method,8)
      nCMo = S%n2Tot
      mCMo = S%n2Tot
      If (Method.eq. 'KS-DFT  ' .or.
     &    Method.eq. 'CASDFT  ' ) Then
         Call Get_iScalar('Multiplicity',iSpin)
         Call Get_cArray('DFT functional',KSDFT,16)
         Call Get_dScalar('DFT exch coeff',CoefX)
         Call Get_dScalar('DFT corr coeff',CoefR)
         ExFac=Get_ExFac(KSDFT)
         CoulFac=One
      Else
         iSpin=0
         ExFac=One
         CoulFac=One
      End If
*
*...  Check the wave function type
*
*                                                                      *
************************************************************************
*                                                                      *
      If ( Method.eq.'RHF-SCF ' .or.
     &     Method.eq.'UHF-SCF ' .or.
     &    (Method.eq.'KS-DFT  '.and.iSpin.eq.1) .or.
     &     Method.eq.'ROHF    ' ) then
         If (lPrint) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: ',Method
            If (Method.eq.'KS-DFT  ') Then
               Write (6,'(2A)') ' Functional type:   ',KSDFT
               Fmt = '(1X,A26,20X,F18.6)'
               Write(6,Fmt)'Exchange scaling factor',CoefX
               Write(6,Fmt)'Correlation scaling factor',CoefR
            End If
            Write (6,*)
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'Corr. WF' ) then
         If (lPrint) Then
            Write (6,*)
            Write (6,*)
     &         ' Wavefunction type: an Aces 2 correlated wavefunction'
            Write (6,*)
         End If
         Gamma_On=.True.
         Call Aces_Gamma()
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'RASSCF  ' .or.
     &          Method.eq.'CASSCF  ' .or.
     &          Method.eq.'CASDFT  ') then
*
         Call Get_iArray('nAsh',nAsh,nIrrep)
         nAct = 0
         Do iIrrep = 0, nIrrep-1
            nAct = nAct + nAsh(iIrrep)
         End Do
         If (nAct.gt.0) lPSO=.true.
*
         nDSO = nDens
         mIrrep=nIrrep
         Call ICopy(nIrrep,nBas,1,mBas,1)
         If (lPrint) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: ', Method
            If (Method.eq.'CASDFT  ')
     &         Write (6,'(2A)') ' Functional type:   ',KSDFT
            Write (6,*)
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'CASSCFSA' .or.
     &          Method.eq.'RASSCFSA' ) then
         Call Get_iArray('nAsh',nAsh,nIrrep)
         nAct = 0
         Do iIrrep = 0, nIrrep-1
            nAct = nAct + nAsh(iIrrep)
         End Do
         If (nAct.gt.0) lPSO=.true.
         nDSO = nDens
         Call Get_iScalar('SA ready',iGo)
         If (iGO.eq.1) lSA=.true.
         mIrrep=nIrrep
         Call ICopy(nIrrep,nBas,1,mBas,1)
         If (lPrint.and.lSA) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: State average ',
     &                         Method(1:6)
            Write (6,*)
         Else If (lPrint) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: ', Method
         End If
         Method='RASSCF  '
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Write (6,*)
         Write (6,*) ' Wavefunction type:',Method
         Write (6,*) ' Illegal type of wave function!'
         Write (6,*) ' ALASKA can not continue'
         Write (6,*)
         Call Quit_OnUserError()
      End If
*
*...  Read the (non) variational 1st order density matrix
*...  density matrix in AO/SO basis
         nsa=1
         If (lsa) nsa=4
         mDens=nsa
         Call mma_allocate(D0  ,nDens,mDens,Label='D0')
         Call mma_allocate(DVar,nDens,mDens,Label='DVar')
         D0  (:,:)=Zero
         DVar(:,:)=Zero
         Call Get_D1ao(D0,nDens)
         Call Get_D1ao_Var(DVar,nDens)
*
         Call ReIndexFrag(D0,nDens,nDens_Valence,nBas,
     &                    nBas_Valence, nIrrep)
         Call ReIndexFrag(DVar,nDens,nDens_Valence,nBas,
     &                    nBas_Valence, nIrrep)
         Call AddFragDens(D0,nDens,nDens_Valence,nBas_Valence)
         Call AddFragDens(DVar,nDens,nDens_Valence,nBas_Valence)
*
         Call mma_allocate(DS,nDens,Label='DS')
         Call mma_allocate(DSVar,nDens,Label='DSVar')
         DS(:)=Zero
         DSVar(:)=Zero
         If(Method.eq.'UHF-SCF ' .or.
     &      Method.eq.'ROHF    ' .or.
     &      Method.eq.'Corr. WF'      ) Then
           Call Get_D1sao(DS,nDens)
           Call Get_D1sao_Var(DSVar,nDens)
         End If
*
*     Unfold density matrix
*
      ij = -1
      Do 10 iIrrep = 0, nIrrep-1
         Do 11 iBas = 1, nBas(iIrrep)
            Do 12 jBas = 1, iBas-1
               ij = ij + 1
               D0   (1+ij,1) = Half*D0   (1+ij,1)
               DVar (1+ij,1) = Half*DVar (1+ij,1)
               DS   (1+ij)   = Half*DS   (1+ij)
               DSVar(1+ij)   = Half*DSVar(1+ij)
 12         Continue
            ij = ij + 1
 11      Continue
 10   Continue
      If (iPrint.ge.99) Then
         RlxLbl='D1AO    '
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],D0)
         RlxLbl='D1AO-Var'
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DVar)
         RlxLbl='DSAO    '
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DS)
         RlxLbl='DSAO-Var'
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DSVar)
      End If
*
*...  Get the MO-coefficients
         If (Method.eq.'UHF-SCF ' .or.
     &       Method.eq.'ROHF    ' .or.
     &       Method.eq.'Corr. WF'      ) Then
            nsa=2
         Else
            nsa=1
            If (lsa) nsa=2
         End If
         kCMO=nsa
         Call mma_allocate(CMO,mCMO,kCMO,Label='CMO')
         Call Get_CMO(CMO(:,1),mCMO)
         If (iPrint.ge.99) Then
            ipTmp1 = 1
            Do iIrrep = 0, nIrrep-1
               Call RecPrt(' CMO''s',' ',
     &                     CMO(ipTmp1,1),nBas_Valence(iIrrep),
     &                     nBas_Valence(iIrrep))
               ipTmp1 = ipTmp1 + nBas_Valence(iIrrep)**2
            End Do
         End If
*
*
*...  Get additional information in the case of a RASSCF wave function
*...  Get the number of inactive, active and frozen orbitals
         If (.not.lpso) Goto 1000
         Call Get_iScalar('nSym',i)
         Call Get_iArray('nIsh',nIsh,i)
         Call Get_iArray('nAsh',nAsh,i)
         Call Get_iArray('nFro',nFro,i)
         If (iPrint.ge.99) Then
            Write (6,*) ' nISh=',nISh
            Write (6,*) ' nASh=',nASh
            Write (6,*) ' nFro=',nFro
         End If
         nAct = 0
         nTst = 0
         Do iIrrep = 0, nIrrep-1
            nAct = nAct + nAsh(iIrrep)
            nTst = nTst + nFro(iIrrep)
         End Do
         If (nTst.ne.0) Then
            Write (6,*)
            Write (6,*) ' No frozen orbitals are allowed!'
            Write (6,*) ' ALASKA can not continue'
            Write (6,*)
            Call Quit_OnUserError()
         End If
*
*...  Get the one body density for the active orbitals
*     (not needed for SA-CASSCF)
         nG1 = nAct*(nAct+1)/2
         nsa=1
         If (lsa) nsa=0
         mG1=nsa
         Call mma_allocate(G1,nG1,mG1,Label='G1')
         If (nsa.gt.0) Then
            Call Get_D1MO(G1(:,1),nG1)
            If (iPrint.ge.99) Call TriPrt(' G1',' ',G1(:,1),nAct)
         End If
*
*...  Get the two body density for the active orbitals
         nG2 = nG1*(nG1+1)/2
         nsa=1
         if (lsa) nsa=2
         mG2=nsa
         Call mma_allocate(G2,nG2,mG2,Label='G2')
         Call Get_P2MO(G2(:,1),nG2)
         If (iPrint.ge.99) Call TriPrt(' G2',' ',G2(1,1),nG1)
         If (lsa) Then

*  CMO1 Ordinary CMO's
*
*  CMO2 CMO*Kappa
*
           Call Get_LCMO(CMO(:,2),mCMO)
           If (iPrint.ge.99) Then
            ipTmp1 = 1
            Do iIrrep = 0, nIrrep-1
               Call RecPrt('LCMO''s',' ',
     &                     CMO(ipTmp1,2),nBas_Valence(iIrrep),
     &                     nBas_Valence(iIrrep))
               ipTmp1 = ipTmp1 + nBas_Valence(iIrrep)**2
            End Do
           End If
*
* P are stored as
*                            _                     _
*   P1=<i|e_pqrs|i> + sum_i <i|e_pqrs|i>+<i|e_pqrs|i>
*   P2=sum_i <i|e_pqrs|i>
*
           Call Get_PLMO(G2(:,2),nG2)
           Call Daxpy_(nG2,1.0d0,G2(:,2),1,G2(:,1),1)
           If(iPrint.ge.99)Call TriPrt(' G2L',' ',G2(:,2),nG1)
           If(iPrint.ge.99)Call TriPrt(' G2T',' ',G2(:,1),nG1)
*
           Call Get_D2AV(G2(:,2),nG2)
           If (iPrint.ge.99) Call TriPrt('G2A',' ',G2(:,2),nG2)
*
*
*  Densities are stored as:
*
*       ipd0 AO:
*
*       D1 = inactive diagonal density matrix
*                                _                 _
*       D2 = <i|E_pq|i> + sum_i <i|E_pq|i>+<i|E_pq|i> + sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> - 1/2 D1
*
*       D3 = sum_i <i|E_pq|i> (active)
*
*       D4 = sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> (inactive)
*
*       G1 = <i|e_ab|i>
*       G2 = sum i <i|e_ab|i>
           Call Getmem('TMP','ALLO','REAL',ipt,2*ndens)
           Call Get_D1I(CMO(1,1),D0(1,1),Work(ipT),
     &                nish,nBas_Valence,nIrrep)
           Call Getmem('TMP','FREE','REAL',ipt,ndens)

           Call dcopy_(nDens_Valence,DVar,1,D0(1,2),1)
           If (.not.isNAC) call daxpy_(ndens,-Half,D0(1,1),1,D0(1,2),1)
           If (iprint.gt.90)Call PrMtrx('D0',[iD0Lbl],iComp,[1],D0)
*
*   This is necessary for the kap-lag
*
           nG1 = nAct*(NAct+1)/2
           Call mma_Allocate(D1AV,nG1,Label='D1AV')
           Call Get_D1AV(D1AV,nG1)
           Call Get_D1A(CMO(1,1),D1AV,D0(1,3),
     &                 nIrrep,nBas_Valence,nish,nash,nDens_Valence)
           Call mma_deallocate(D1AV)
*
           Call Get_DLAO(D0(1,4),nDens)
         End If
         If (iPrint.ge.99) Call TriPrt(' G2',' ',G2(1,1),nG1)
*
*...  Close 'RELAX' file
1000     Continue
*
*...  Epilogue, end
      Return
      End

      SubRoutine ReIndexFrag(Array, nDens, nDens_Valence, nBas,
     &                       nBas_Valence,nIrrep)
************************************************************************
*                                                                      *
* Input: Array(nDens) filled up to Array(nDens_Valence)                *
* Output: Reindexed Array so fragment densities can be put at the      *
*         proper places.                                               *
* Only needed in case of symmetry.                                     *
*                                                                      *
************************************************************************
      Implicit None
#include "real.fh"
      Integer nDens, nDens_Valence,nIrrep
      Real*8  Array(nDens)
      Integer nBas(0:7), nBas_Valence(0:7)
      Integer indexLarge,indexSmall,iIrrep

      If(nIrrep.eq.1) return

      indexLarge = nDens + 1
      indexSmall = nDens_Valence + 1
      Do iIrrep = nIrrep-1, 0, -1
* calculate the position in the hypothetical Array(nDens_Valence)
* and in the needed Array(nDens)
        indexLarge = indexLarge - nBas(iIrrep)
        indexSmall = indexSmall - nBas_Valence(iIrrep)
* move the data
        call dcopy_(nBas_Valence(iIrrep), Array(indexSmall), 1,
     &                                   Array(indexLarge), 1)
        call dcopy_(nBas_Valence(iIrrep),[Zero],0,Array(indexSmall),1)
      End Do
      Return
      End

      SubRoutine AddFragDens(Array, nDens, nDens_Valence, nBas_Valence)
************************************************************************
*                                                                      *
* Input: Array(Size) filled with valence density at proper positions   *
*        for a density including fragments.                            *
* Output: Updated with fragment densities at their proper positions.   *
*                                                                      *
************************************************************************
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iOper
      Implicit None
#include "real.fh"
#include "WrkSpc.fh"
      Integer nDens,nDens_Valence
      Real*8  Array(nDens)
      Integer nBas_Valence(0:7)
      Logical EnergyWeight
      Integer iPrint,maxDens,iCnttp,ipFragDensAO,iDpos,ipFragDensSO
      Integer i,j,iCnt,iFpos,iFD,mdc,iIrrep,nBasC
      Real*8  rDummy(1)

      If(nIrrep.ne.1) Then
        write(6,*) 'AddFragDens: Symmetry not implemented yet'
        Call Abend()
      End If
*
      iPrint=0
*
* Each fragment needs it''s (symmetrized) density matrix added along the diagonal
* This density matrix first has to be constructed from the MO coefficients
* so allocate space for the largest possible density matrix
      maxDens = 0
      Do iCnttp = 1, nCnttp
        If (dbsc(iCnttp)%nFragType.gt.0) maxDens = Max(maxDens,
     &                        dbsc(iCnttp)%nFragDens
     &                      *(dbsc(iCnttp)%nFragDens+1)/2)
      End Do
      Call GetMem('FragDSO','Allo','Real',ipFragDensSO,maxDens)
c     If(nIrrep.ne.1) Then
c       Call GetMem('FragDAO','Allo','Real',ipFragDensAO,maxDens)
c     Else
        ipFragDensAO = ipFragDensSO
c     End If

      iDpos = 1 ! position in the total density matrix
      Do iIrrep = 0, nIrrep - 1
        nBasC = nBas_Valence(iIrrep)
        iDpos = iDpos + nBasC*(nBasC+1)/2
        mdc = 0
        Do 1000 iCnttp = 1, nCnttp
          If(dbsc(iCnttp)%nFragType.le.0) Then
            mdc = mdc + dbsc(iCnttp)%nCntr
            Go To 1000
          End If

* construct the density matrix
          EnergyWeight = .false.
          Call MakeDens(dbsc(iCnttp)%nFragDens,
     &                  dbsc(iCnttp)%nFragEner,
     &                  dbsc(iCnttp)%FragCoef,
     &                  rDummy,
     &                  EnergyWeight,Work(ipFragDensAO))
* create the symmetry adapted version if necessary
* (fragment densities are always calculated without symmetry)
C         If(nIrrep.ne.1) Call SymmDens(Work(ipFragDensAO),
C    &      Work(ipFragDensSO))
          If(iPrint.ge.99) Call TriPrt('Fragment density',' ',
     &      Work(ipFragDensSO),dbsc(iCnttp)%nFragDens)
          Do iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
* only add fragment densities that are active in this irrep
* => the following procedure still has to be verified thoroughly
*    but appears to be working
            If(iAnd(dc(mdc)%iChCnt,iIrrep).eq.iOper(iIrrep)) Then
* add it at the correct location in the large custom density matrix
              iFpos = 1
c              ! position in fragment density matrix
              Do i = 1, dbsc(iCnttp)%nFragDens
                iDpos = iDpos + nBasC
                Do j = 0, i-1
                  Array(iDpos + j) = Work(ipFragDensSO + iFpos + j - 1)
                End Do
                iDpos = iDpos + i
                iFpos = iFpos + i
              End Do
              nBasC = nBasC + dbsc(iCnttp)%nFragDens
            End If
          End Do
 1000   Continue
      End Do
      If(iPrint.ge.19) Then
        iFD = 1
        Do iIrrep = 0, nIrrep - 1
          Call TriPrt('Combined density',' ',Array(iFD),nBas(iIrrep))
          iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
        End Do
      End If
      Call GetMem('FragDSO','Free','Real',ipFragDensSO,maxDens)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nDens_Valence)
      End
