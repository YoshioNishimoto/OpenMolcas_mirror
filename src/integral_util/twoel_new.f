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
! Copyright (C) 1990,1991,1993, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine TwoEl_Sym_New(iS_,jS_,kS_,lS_,
     &           Coor,
     &           iAnga,iCmp,iShell,iShll,iAO,iAOst,
     &           NoInts,iStb,jStb,kStb,lStb,
     &           nAlpha,iPrInc, nBeta,jPrInc,
     &           nGamma,kPrInc,nDelta,lPrInc,
     &           Data1,mData1,nData1,Data2,mData2,nData2,
     &           IJeqKL,kOp,
     &           Dij,mDij,mDCRij,Dkl,mDkl,mDCRkl,Dik,mDik,mDCRik,
     &           Dil,mDil,mDCRil,Djk,mDjk,mDCRjk,Djl,mDjl,mDCRjl,
     &           Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl,
     &           FckTmp,nFT,Zeta,ZInv,IndZet,Kappab,P,nZeta,
     &           Eta,EInv,IndEta,Kappcd,Q,nEta,
     &           SOInt,nSOInt,Wrk,nWork2,
     &           Shijij,nHRRAB,nHRRCD,Aux,nAux)
!***********************************************************************
!                                                                      *
! Object: to generate the SO integrals for four fixed centers and      *
!         fixed basis set types.                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to use Schwartz inequality for pre-  *
!          screening, July 1991.                                       *
!          Modified for direct SCF, January '93                        *
!***********************************************************************
      use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Phase_Info
      use Gateway_Info, only: ThrInt, CutInt
      use Symmetry_Info, only: nIrrep
      use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
      use Int_Options, only: ExFac, Thize, W2Disc, IntOnly=>PreSch
      use Int_Options, only: Disc_Mx, Disc, Quad_ijkl
      use k2_arrays, only: TwoHam=>pFq, Dens=>pDq
      use Breit, only: nOrdOp, nComp
      use Constants
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
#include "print.fh"
#include "twoswi.fh"
      Real*8 SOInt(iBasi*jBasj*kBask*lBasl,nSOInt)
      Real*8 Coor(3,4), CoorM(3,4), CoorAC(3,2),
     &       Data1(mData1,nData1),Data2(mData2,nData2),
     &       Zeta(nZeta), ZInv(nZeta), KappAB(nZeta), P(nZeta,3),
     &       Eta(nEta),   EInv(nEta),  KappCD(nEta),  Q(nEta,3),
     &       Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj),
     &       Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl),
     &       Wrk(nWork2), QInd(2), Aux(nAux), FckTmp(nFT),
     &       Dij(mDij,mDCRij),Dkl(mDkl,mDCRkl),Dik(mDik,mDCRik),
     &       Dil(mDil,mDCRil),Djk(mDjk,mDCRjk),Djl(mDjl,mDCRjl)
      Integer iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iStabN(0:7),
     &        iStabM(0:7), IndZet(nZeta), IndEta(nEta),
     &        iAO(4), iAnga(4), iCmp(4),
     &        iShell(4), iShll(4), kOp(4), iAOst(4), jOp(6), iWR(2)
      Logical NoPInts, Shijij, AeqB, CeqD, AeqC, ABeqCD,
     &        EQ, Copy, NoCopy,Do_TnsCtl,
     &        IJeqKL,IeqK,JeqL,
     &        lEmpty, Prescreen_On_Int_Only, DoCoul, DoExch,
     &        Scrij, Scrkl, Scrik, Scril, Scrjk, Scrjl,
     &        Batch_On_Disk,
     &        NoInts, DoAOBatch, All_Spherical
#ifdef _DEBUGPRINT_
      Character ChOper(0:7)*3
      Data ChOper/' E ',' x ',' y ',' xy',' z ',' xz',' yz','xyz'/
#endif
      Data Copy/.True./, NoCopy/.False./, jOp/0,0,0,0,0,0/
#include "SysDef.fh"
      External EQ, lEmpty

      Interface
      Integer Function iGet(A,n)
      Integer :: n
      Real*8, Target :: A(*)
      End Function iGet
      End Interface
!
!     Declaration of statement functions to compute canonical index
!
      Integer :: ixyz, nabSz
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
!
      If (nOrdOp/=0) Then
         Write (6,*) 'Breit two-electron integrals not implemented yet'
         Write (6,*) 'Symmetry adaptation different since the operator'
         Write (6,*) 'is not symmetric.'
      End If
      Call TwoEl_Sym_New_Internal(Data1,Data2)
!
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iS_)
         Call Unused_integer(jS_)
         Call Unused_integer(kS_)
         Call Unused_integer(lS_)
         Call Unused_integer(iPrInc)
         Call Unused_integer(kPrInc)
         Call Unused_real_array(Dens)
      End If
!
!     This is to allow type punning without an explicit interface
      Contains
      Subroutine TwoEl_Sym_New_Internal(Data1,Data2)
      Real*8, Target :: Data1(mData1,nData1),Data2(mData2,nData2)
      Integer, Pointer :: iData1(:),iData2(:)
!
#ifdef _DEBUGPRINT_
      iRout = 12
      iPrint = nPrint(iRout)
#endif
!
      All_Spherical=Shells(iShll(1))%Prjct.and.
     &              Shells(iShll(2))%Prjct.and.
     &              Shells(iShll(3))%Prjct.and.
     &              Shells(iShll(4))%Prjct
      QInd(1)=Quad_ijkl
      RST_triplet=One
!
      Do_tnsctl=.False.
      kabcd=0
!
#ifdef _DEBUGPRINT_
      Do iData = 1, nData1
         Call RecPrt('Twoel: Data1',' ',
     &                Data1(1,iData),nAlpha*nBeta,nDArray)
      End Do
      Do iData = 1, nData2
         Call RecPrt('Twoel: Data2',' ',
     &               Data2(1,iData),nGamma*nDelta,nDArray)
      End Do
      Call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
      Call RecPrt('Coeff2',' ',Coeff2,nBeta,jBasj)
      Call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
      Call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
      Call RecPrt('Coor',' ',Coor,3,4)
#endif
!
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      iSmAng=la+lb+lc+ld
      LmbdT=0
      nijkl = iBasi*jBasj*kBask*lBasl*nComp
      nab = iCmp(1)*iCmp(2)
      ncd = iCmp(3)*iCmp(4)
      nabcd = nab*ncd
      nInts =nijkl*nabcd
      SOInt(:,:) = Zero
      ipAOInt=1
      iW3=1+nInts
      iW4=1
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Find the Double Coset Representatives for center A and B
!
      Call DCR(LmbdR,dc(iStb)%iStab,dc(iStb)%nStab,
     &               dc(jStb)%iStab,dc(jStb)%nStab,iDCRR,nDCRR)
      u = DBLE(dc(iStb)%nStab)
      v = DBLE(dc(jStb)%nStab)
#ifdef _DEBUGPRINT_
      If (iPrint.ge.9) Write (6,'(20A)') ' {R}=(',
     &      (ChOper(iDCRR(i)),',',i=0,nDCRR-1),')'
#endif
!
!-----Find stabilizer for center A and B
!
      Call Inter(dc(iStb)%iStab,dc(iStb)%nStab,
     &           dc(jStb)%iStab,dc(jStb)%nStab,iStabM,lStabM)
!
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Find the Double Coset Representatives for center C and D.
!
      Call DCR(LmbdS,dc(kStb)%iStab,dc(kStb)%nStab,
     &               dc(lStb)%iStab,dc(lStb)%nStab,iDCRS,nDCRS)
      w = DBLE(dc(kStb)%nStab)
      x = DBLE(dc(lStb)%nStab)
#ifdef _DEBUGPRINT_
      If (iPrint.ge.9) Write (6,'(20A)') ' {S}=(',
     &      (ChOper(iDCRS(i)),',',i=0,nDCRS-1),')'
#endif
!
!-----Find stabilizer for center C and D
!
      Call Inter(dc(kStb)%iStab,dc(kStb)%nStab,
     &           dc(lStb)%iStab,dc(lStb)%nStab,iStabN,lStabN)
!                                                                      *
!***********************************************************************
!                                                                      *
!
!-----Find the Double Coset Representatives for the two charge
!     distributions.
!
      Call DCR(LmbdT,iStabM,lStabM,iStabN,lStabN,iDCRT,nDCRT)
!                                                                      *
!***********************************************************************
!                                                                      *
      kOp(1)=NrOpr(0)
      call dcopy_(3,Coor(1,1),1,CoorM(1,1),1)
      Do 100 lDCRR = 0, nDCRR-1
         kOp(2)=NrOpr(iDCRR(lDCRR))
         Call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
         AeqB = EQ(CoorM(1,1),CoorM(1,2))
!
         lDCR1=NrOpr(iDCRR(lDCRR))+1
!
         vijij=Data1(ip_abMax(nZeta),lDCR1)
!
!switch (to generate better start orbitals...)
         If (NDDO .AND. .NOT.AeqB) Go To 100
!switch
         MxDCRS = nDCRS-1
         Do 200 lDCRS = 0, MxDCRS
            RS_doublet=DBLE(lDCRS*nDCRR+lDCRR+1)
            call dcopy_(3,Coor(1,3),1,CoorM(1,3),1)
            Call OA(iDCRS(lDCRS),Coor(1:3,4),CoorM(1:3,4))
            CeqD = EQ(Coor(1,3),CoorM(1,4))
!
!switch (to generate better start orbitals...)
            If (NDDO .AND. .NOT.CeqD) Go To 200
!
            lDCR2=NrOpr(iDCRS(lDCRS))+1
!
!-----------Pickup estimated largest integral value (AO)
!
            vijkl = vijij * Data2(ip_abMax( nEta),lDCR2)
            Do 300 lDCRT = 0, nDCRT-1
               ipAOInt=1
               iW3=1+nInts
               RS_doublet=DBLE(lDCRS*nDCRR+lDCRR+1)
               RST_triplet=DBLE(lDCRT*nDCRR*nDCRS)+RS_doublet
               QInd(2)=RST_triplet
!              Write (*,*) QInd(1), QInd(2)
               iDCRTS=iEor(iDCRT(lDCRT),iDCRS(lDCRS))
               Call OA(iDCRTS,Coor(1:3,4),CoorM(1:3,4))
               Call OA(iDCRT(lDCRT),Coor(1:3,3),CoorM(1:3,3))
               AeqC = EQ(CoorM(1,1),CoorM(1,3))
               ABeqCD = AeqB .and. CeqD .and. AeqC
               If (ABeqCD .and. Mod(iSmAng,2).eq.1) Go To 300
!lwj
!lwj           For Spherical Gaussians, batches like
!lwj           (DS|SS), (FP|SS) and (FS|PS) vanish as well
               If (ABeqCD .and. All_Spherical .and.
     &                2*Max(la,lb,lc,ld).gt.iSmAng) Go To 300
!lwj
!
#ifdef _DEBUGPRINT_
               If (iPrint.ge.9) Write (6,'(6A)')
     &         ' R=',ChOper(iDCRR(lDCRR)),
     &         ', S=',ChOper(iDCRS(lDCRS)),
     &         ', T=',ChOper(iDCRT(lDCRT))
#endif
!
               kOp(3) = NrOpr(iDCRT(lDCRT))
               kOp(4) = NrOpr(iEor(iDCRT(lDCRT),iDCRS(lDCRS)))
!
               ix1 = 1
               iy1 = 1
               iz1 = 1
               ix2 = iPhase(1,iDCRT(lDCRT))
               iy2 = iPhase(2,iDCRT(lDCRT))
               iz2 = iPhase(3,iDCRT(lDCRT))
               lDCRE_=0
               lDCRT_=iDCRT(lDCRT)
!
!------------- Find index to desymmetrized Dij, Dkl, Dik, Dil, Djk, and
!              Djl. Some care has to be taken here. Assume that there
!              are two operators, T and S which generates the center
!              pairs A,T(B) and A,S(B). If these pairs are symmetry
!              related we will only
!
               If (DoFock) Then
!--------------Dij
               iR = iDCRR(lDCRR)
               jOp(1) = NrOpr(iR) + 1
!--------------Dkl
               iS = iDCRS(lDCRS)
               jOp(2)= NrOpr(iS) + 1
!--------------Dik
               iT  = iDCRT(lDCRT)
               jOp(3)= NrOpr(iT) + 1
!--------------Dil
               iTS = iEor(iT,iS)
               jOp(4)= NrOpr(iTS) + 1
!--------------Djk
               iRT = iEor(iR,iT)
               jOp(5)= NrOpr(iRT) + 1
!--------------Djl
               iRTS= iEor(iRT,iS)
               jOp(6)= NrOpr(iRTS) + 1
!
               End If ! DoFock
!                                                                      *
!***********************************************************************
!                                                                      *
!--------------Prescreening at group level
!
!--------------Select if batch should be written on the disc at the
!              first iteration and then later on read.
!
               Batch_On_Disk = (vijkl.gt.Thize) .and.
     &                (Disc+DBLE(nInts+2+2/RtoI).le.Disc_Mx)
!
!--------------Set prescreening level
!
!              IntOnly = T  prescreening on integral value only
!              IntOnly = F  prescreening on integral & density matrix
!
!              If integral batch is written to disc, no prescreening
!              since prescreening on integrals only was done in k2 loop.
!
               Prescreen_On_Int_Only = IntOnly
               If (DoIntegrals) Prescreen_On_Int_Only=.True.
               If (Batch_On_Disk) Prescreen_On_Int_Only = .True.
!
!--------------Prescreening based on the 1st order density in AO
!              basis (contracted). Observe that this is a rough
!              estimate of if there will be any contribution due to
!              this integral batch.
!
!--------------Special care here for RS, RT and ST degeneracy.
!
!............. Get maximum density elements in AO basis.
!
               If (DoFock) Then
                  vij = Dij(mDij,jOp(1))
                  vkl = Dkl(mDkl,jOp(2))
                  vik = Dik(mDik,jOp(3))
                  vil = Dil(mDil,jOp(4))
                  vjk = Djk(mDjk,jOp(5))
                  vjl = Djl(mDjl,jOp(6))
!                 Coulomb contributions
                  Scrkl = vij*vijkl.ge.ThrInt
                  Scrij = vkl*vijkl.ge.ThrInt
                  DoCoul  = Scrij .or. Scrkl
!
!                 Exchange contributions
                  vik = vik/Four
                  vil = vil/Four
                  vjk = vjk/Four
                  vjl = vjl/Four
                  Scrjl = vik*vijkl.ge.ThrInt
                  Scrjk = vil*vijkl.ge.ThrInt
                  Scril = vjk*vijkl.ge.ThrInt
                  Scrik = vjl*vijkl.ge.ThrInt
                  DoExch = Scrjl .or. Scrjk .or. Scril .or. Scrik
                  If (FckNoClmb) DoCoul=.False.
                  If (FckNoExch) DoExch=.False.
               Else
                  DoCoul=.False.
                  DoExch=.False.
               End If
               DoAOBatch=(DoIntegrals.and.vijkl.gt.CutInt).or.
     &                   (DoFock.and.(DoCoul.or.DoExch)) .or.
     &                   (Batch_On_Disk.and.W2Disc)
!
!--------------Branch out if crude estimate indicates no contributions!
!
               If (.Not.DoAOBatch) Then
                  If (.Not.Batch_On_Disk) Then
!
!                    AO batch is not on the disk
!
                     Go To 300
                  Else If (Batch_On_Disk.and..Not.W2Disc) Then
!
!                    AO batch is on disk! Do a no copy read to
!                    position the next batch on the disc.
!
 1111                Continue
                     Call iRBuf(iWR,2,Copy)
                     Call dRBuf(QInd,2,Copy)
                     Call Store_QLast(QInd)
                     kInts=iWR(1)
                     mInts=iWR(2)
                     If (QInd(1).eq.Quad_ijkl .and.
     &                   QInd(2).eq.RST_triplet) Then
                        If (kInts.ne.nInts) Then
                           Call WarningMessage(2,
     &                                 'Twoel: kInts.ne.nInts!')
                           Write (6,*) 'Twoel: kInts,mInts,nInts=',
     &                                         kInts,mInts,nInts
                           Write (6,*) 'Index,1:',QInd(1),QInd(2),
     &                                  Quad_ijkl,RST_triplet
                           Call Abend()
                        End If
                        If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,
     &                                             NoCopy)
                        Disc = Disc + DBLE(2/RtoI + 2 + mInts)
                        Go To 300
                     Else If (QInd(1).le.Quad_ijkl) Then
                        If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,
     &                                             NoCopy)
                        Disc = Disc + DBLE(2/RtoI + 2 + mInts)
                        Go To 1111
                     Else
                        Call WarningMessage(2,
     &                              'Twoel: batch is lost!')
                        Write (6,*) 'Index,1:',QInd(1),QInd(2),
     &                               Quad_ijkl,RST_triplet
                        Call Abend()
                     End If
                  End If
               End If
!
               If (DoFock) Then
                  If (iShell(1).ge.iShell(2)) Then
                     ij1 = iBasi
                     ij2 = jBasj
                     ij3 = iCmp(1)
                     ij4 = iCmp(2)
                  Else
                     ij1 = jBasj
                     ij2 = iBasi
                     ij3 = iCmp(2)
                     ij4 = iCmp(1)
                  End If
                  If (iShell(3).ge.iShell(4)) Then
                     kl1 = kBask
                     kl2 = lBasl
                     kl3 = iCmp(3)
                     kl4 = iCmp(4)
                  Else
                     kl1 = lBasl
                     kl2 = kBask
                     kl3 = iCmp(4)
                     kl4 = iCmp(3)
                  End If
                  If (iShell(1).ge.iShell(3)) Then
                     ik1 = iBasi
                     ik2 = kBask
                     ik3 = iCmp(1)
                     ik4 = iCmp(3)
                  Else
                     ik1 = kBask
                     ik2 = iBasi
                     ik3 = iCmp(3)
                     ik4 = iCmp(1)
                  End If
                  If (iShell(1).ge.iShell(4)) Then
                     il1 = iBasi
                     il2 = lBasl
                     il3 = iCmp(1)
                     il4 = iCmp(4)
                  Else
                     il1 = lBasl
                     il2 = iBasi
                     il3 = iCmp(4)
                     il4 = iCmp(1)
                  End If
                  If (iShell(2).ge.iShell(3)) Then
                     jk1 = jBasj
                     jk2 = kBask
                     jk3 = iCmp(2)
                     jk4 = iCmp(3)
                  Else
                     jk1 = kBask
                     jk2 = jBasj
                     jk3 = iCmp(3)
                     jk4 = iCmp(2)
                  End If
                  If (iShell(2).ge.iShell(4)) Then
                     jl1 = jBasj
                     jl2 = lBasl
                     jl3 = iCmp(2)
                     jl4 = iCmp(4)
                  Else
                     jl1 = lBasl
                     jl2 = jBasj
                     jl3 = iCmp(4)
                     jl4 = iCmp(2)
                  End If
               End If ! DoFock
!
!--------------Branch point for partial integral storage.
!
               If (Batch_On_Disk.and..Not.W2Disc) Go To 6767
!                                                                      *
!***********************************************************************
!                                                                      *
!              Here if the AO batch will be computed !
!
!--------------Compute actual size of the {a0|c0} block
!
               mabMin=nabSz(Max(la,lb)-1)+1
               If (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nabSz(la+lb-1)+1
               mabMax=nabSz(la+lb)
               mcdMin=nabSz(Max(lc,ld)-1)+1
               If (EQ(CoorM(1,3),CoorM(1,4))) mcdMin = nabSz(lc+ld-1)+1
               mcdMax=nabSz(lc+ld)
               mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
!
!--------------Find the proper centers to start of with the angular
!              momentum on. If la.eq.lb there will exist an
!              ambiguity to which center that angular momentum should
!              be accumulated on. In that case we will use A and C of
!              the order as defined by the basis functions types.
!
               If (iAnga(1).ge.iAnga(2)) Then
                  call dcopy_(3,CoorM(1,1),1,CoorAC(1,1),1)
               Else
                  call dcopy_(3,CoorM(1,2),1,CoorAC(1,1),1)
               End If
               If (iAnga(3).ge.iAnga(4)) Then
                   call dcopy_(3,CoorM(1,3),1,CoorAC(1,2),1)
               Else
                   call dcopy_(3,CoorM(1,4),1,CoorAC(1,2),1)
               End If
!
!--------------Set flags if triangularization will be used
!
               IeqK = EQ(CoorM(1,1),CoorM(1,3))
               JeqL = EQ(CoorM(1,2),CoorM(1,4))
               IJeqKL = IeqK .and. JeqL
!
!--------------Loops to partition the primitives
!
!              Reset pointer ipAOInt if we need to reserve speacial
!              space for the contracted integrals.
!
               IncZet=nAlpha*jPrInc
               IncEta=nGamma*lPrInc
               If (nZeta.ne.IncZet.or.nEta.ne.IncEta) Then
                  mWork2 = nWork2 - nijkl*mabcd
                  ipAOInt=1+nijkl*mabcd
               Else
                  mWork2 = nWork2
                  ipAOInt=1
               End If
!
!
               nZeta_Tot=iGet(Data1(ip_IndZ(1,nZeta),lDCR1),nZeta+1)
               nEta_Tot =iGet(Data2(ip_IndZ(1,nEta ),lDCR2),nEta +1)
!
               iZ13_=ip_HrrMtrx(nZeta)+(NrOpr(lDCRE_)*nHRRAB)/nIrrep
               iE13_=ip_HrrMtrx( nEta)+(NrOpr(lDCRT_)*nHRRCD)/nIrrep
!
               kabcd=0
               Do_TnsCtl=.False.
               NoInts   =.True.
               NoPInts  =.True.
               ipAOInt_=ipAOInt
               iW4_=iW4
!
               Do iZeta = 1, nZeta_Tot, IncZet
                  mZeta=Min(IncZet,nZeta_Tot-iZeta+1)
                  If (lEmpty(Coeff2,nBeta,nBeta,jBasj)) Cycle
!
               Do iEta  = 1, nEta_Tot,  IncEta
                  mEta=Min(IncEta,nEta_Tot-iEta+1)
                  If (lEmpty(Coeff4,nDelta,nDelta,lBasl)) Cycle
!
                  Call C_F_Pointer(C_Loc(Data1(ip_IndZ(1,nZeta),lDCR1)),
     &                             iData1,[nZeta])
                  Call C_F_Pointer(C_Loc(Data2(ip_IndZ(1,nEta ),lDCR2)),
     &                             iData2,[nEta])
                  Call DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,
     &                        nZeta_Tot,nEta_Tot,
     &                        Data1(1,lDCR1),mData1,
     &                        Data2(1,lDCR2),mData2,
     &                        nAlpha,nBeta,nGamma,nDelta,
     &                        iData1,
     &                        Zeta,ZInv,P,KappAB,IndZet,
     &                        iData2,
     &                        Eta,EInv,Q,KappCD,IndEta,
     &                        ix1,iy1,iz1,ix2,iy2,iz2,ThrInt,CutInt,
     &                        vij,vkl,vik,vil,vjk,vjl,
     &                        Prescreen_On_Int_Only,NoInts,iAnga,
     &                        CoorM,CoorAC,
     &                        mabMin,mabMax,mcdMin,mcdMax,nijkl/nComp,
     &                        nabcd,mabcd,Wrk,ipAOInt_,iW4_,
     &                        nWork2,mWork2,
     &                        Data1(iZ13_,lDCR1),Data2(iE13_,lDCR2),
     &                        la,lb,lc,ld,
     &                        iCmp,iShll,NoPInts,
     &                        Dij(1,jOp(1)),mDij,Dkl(1,jOp(2)),mDkl,
     &                        Do_TnsCtl,kabcd,
     &                        Coeff1,iBasi,Coeff2,jBasj,
     &                        Coeff3,kBask,Coeff4,lBasl)
                  Nullify(iData1,iData2)
!
               End Do
               End Do
!
!              Integrals are now returned in Wrk(ipAOInt) or Wrk(iW4)
!
               If (NoPInts) Then
                  If (W2Disc) Then
                     If (Batch_On_Disk) Then
!
!                       If the batch was supposed to be on disk make
!                       a mark.
!
                        mInts=0
!
                        iWR(1)=nInts
                        iWR(2)=mInts
                        Call iWBuf(iWR,RtoI)
                        Call dWBuf(QInd,2)
                        Call Store_QLast(QInd)
!
                        Disc = Disc + DBLE(2/RtoI + 2 + mInts)
                     End If
                  End If
                  Go To 300
               End If
!
!------------- Multiply with factors due to summation over DCR
!
               If (MolWgh.eq.1) Then
                  FactNd = DBLE(nIrrep) / DBLE(LmbdT)
               Else If (MolWgh.eq.0) Then
                  FactNd = u*v*w*x / DBLE(nIrrep**3 * LmbdT)
               Else
                  FactNd = sqrt(u*v*w*x)/DBLE(nirrep*lmbdt)
               End If
               If (FactNd.ne.One) Call DScal_(kabcd*nijkl,FactNd,
     &                                       Wrk(iW4),1)
!
!--------------Apply the transfer equation and transform the spherical
!              harmonic gaussian.
!
               If (Do_TnsCtl) Then
                  Call TnsCtl(Wrk(iW4),nWork2,CoorM,
     &                        nijkl,mabMax,mabMin,mcdMax,mcdMin,
     &                        Data1(iZ13_,lDCR1),Data2(iE13_,lDCR2),
     &                        la,lb,lc,ld,
     &                        iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                        iShll(1),iShll(2),iShll(3),iShll(4),i_Int)
                  ipAOInt=i_Int
                  If (i_Int.eq.1) Then
                     iW3=1+nijkl*nabcd
                  Else
                     iW3=1
                  End If
               Else
!
!---------------- Undo the late Cntrct
!
                  call dcopy_(nabcd*nijkl,Wrk(ipAOInt),1,Wrk(iW3),1)
                  Call DGeTMO(Wrk(iW3),nabcd,nabcd,
     &                        nijkl,Wrk(ipAOInt),nijkl)
!
               End If
!
!--------------Branch point for partial integral storage
!
               If (Batch_On_Disk.and.W2Disc) Then
!
!---------------- Write integrals to current position on disc.
!
                  iOpt=0 ! Always Packing, not run length
                  Call PkR8(iOpt,nInts,nByte,Wrk(ipAOInt),Wrk(iW3))
                  mInts=(nByte+RtoB-1)/RtoB
!
                  iWR(1)=nInts
                  iWR(2)=mInts
                  Call iWBuf(iWR,2)
                  Call dWBuf(QInd,2)
                  Call Store_QLast(Qind)
                  Call dWBuf(Wrk(iW3),mInts)
                  Disc = Disc + DBLE(2/RtoI + 2 + mInts)
!
               End If
!
 6767          Continue
               If (Batch_On_Disk.and..Not.W2Disc) Then
 1112             Continue
                  Call iRBuf(iWR,2,Copy)
                  Call dRBuf(QInd,2,Copy)
                  Call Store_QLast(QInd)
                  kInts=iWR(1)
                  mInts=iWR(2)
                  If (QInd(1).eq.Quad_ijkl .and.
     &                QInd(2).eq.RST_triplet) Then
                     If (kInts.ne.nInts) Then
                        Call WarningMessage(2,
     &                              'Twoel: kInts.ne.nInts!')
                        Write (6,*) 'Twoel: kInts,mInts,nInts=',
     &                                      kInts,mInts,nInts
                        Write (6,*) 'Index,1:',QInd(1),QInd(2),
     &                               Quad_ijkl,RST_triplet
                        Call Abend()
                     End If
                     If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,Copy)
                     Disc = Disc + DBLE(2/RtoI + 2 + mInts)
                     If (mInts.eq.0) Go To 300
                  Else If (QInd(1).lt.Quad_ijkl) Then
                     If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,NoCopy)
                     Disc = Disc + DBLE(2/RtoI + 2 + mInts)
                     Go To 1112
                  Else
                     Call WarningMessage(2,'Twoel: batch is lost!')
                     Write (6,*) 'Index,1:',QInd(1),QInd(2),
     &                            Quad_ijkl,RST_triplet
                     Call Abend()
                  End If
!
                  iOpt=0 ! Always packing, not run length
                  Call UpkR8(iOpt,nInts,nByte,Wrk(iW3),Wrk(ipAOInt))
               End If
!
!--------------Accumulate contributions directly to the symmetry
!              adapted Fock matrix.
!
               mWork3=nWork2-iW3+1
               If (DoFock)
     &         Call FckAcc(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                     Shijij,iShll,iShell,kOp,nijkl/nComp,
     &                     Wrk(ipAOInt),TwoHam,Size(TwoHam),
     &                     Wrk(iW3),mWork3,
     &                     iAO,iAOst,
     &                     iBasi,jBasj,kBask,lBasl,
     &                     Dij(1,jOp(1)),ij1,ij2,ij3,ij4,
     &                     Dkl(1,jOp(2)),kl1,kl2,kl3,kl4,
     &                     Dik(1,jOp(3)),ik1,ik2,ik3,ik4,
     &                     Dil(1,jOp(4)),il1,il2,il3,il4,
     &                     Djk(1,jOp(5)),jk1,jk2,jk3,jk4,
     &                     Djl(1,jOp(6)),jl1,jl2,jl3,jl4,
     &                     FckTmp,nFT,DoCoul,DoExch,ExFac)
!
!              Transform from AO basis to SO basis
!
               If (DoIntegrals)
     &         Call SymAdp(iAnga, iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                     Shijij,iShll,iShell,iAO,kOp,nijkl,
     &                     Aux,nAux,Wrk(ipAOInt),SOInt,nSOInt,NoInts)
!
 300        Continue
 200     Continue
 100  Continue
      End Subroutine TwoEl_Sym_New_Internal
!
      End
!#define _DEBUGPRINT_
      SubRoutine TwoEl_NoSym_New(iS_,jS_,kS_,lS_,
     &           Coor,
     &           iAnga,iCmp,iShell,iShll,iAO,iAOst,
     &           NoInts,iStb,jStb,kStb,lStb,
     &           nAlpha,iPrInc, nBeta,jPrInc,
     &           nGamma,kPrInc,nDelta,lPrInc,
     &           Data1,mData1,nData1,Data2,mData2,nData2,
     &           IJeqKL,kOp,
     &           Dij,mDij,mDCRij,Dkl,mDkl,mDCRkl,Dik,mDik,mDCRik,
     &           Dil,mDil,mDCRil,Djk,mDjk,mDCRjk,Djl,mDjl,mDCRjl,
     &           Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl,
     &           FckTmp,nFT,Zeta,ZInv,IndZet,Kappab,P,nZeta,
     &           Eta,EInv,IndEta,Kappcd,Q,nEta,
     &           SOInt,nSOInt,Wrk,nWork2,
     &           Shijij,nHRRAB,nHRRCD,Aux,nAux)
!***********************************************************************
!                                                                      *
! Object: to generate the SO integrals for four fixed centers and      *
!         fixed basis set types.                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to use Schwartz inequality for pre-  *
!          screening, July 1991.                                       *
!          Modified for direct SCF, January '93                        *
!***********************************************************************
      use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Gateway_Info, only: ThrInt, CutInt
      use Symmetry_Info, only: nIrrep
      use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
      use Int_Options, only: ExFac, Thize, W2Disc, IntOnly=>PreSch
      use Int_Options, only: Disc_Mx, Disc, Quad_ijkl
      use k2_arrays, only: TwoHam=>pFq, Dens=>pDq
      use Breit, only: nComp
      use Constants
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
#include "print.fh"
#include "twoswi.fh"
      Real*8 Coor(3,4), CoorAC(3,2),
     &       Data1(mData1),Data2(mData2),
     &       Zeta(nZeta), ZInv(nZeta), KappAB(nZeta), P(nZeta,3),
     &       Eta(nEta),   EInv(nEta),  KappCD(nEta),  Q(nEta,3),
     &       Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj),
     &       Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl),
     &       SOInt(nSOInt),Wrk(nWork2), QInd(2), Aux(nAux),FckTmp(nFT),
     &       Dij(mDij,mDCRij),Dkl(mDkl,mDCRkl),Dik(mDik,mDCRik),
     &       Dil(mDil,mDCRil),Djk(mDjk,mDCRjk),Djl(mDjl,mDCRjl)
      Integer IndZet(nZeta),IndEta(nEta),iAO(4), kOp(4),
     &        iAnga(4), iCmp(4), iShell(4), iShll(4), iAOst(4), iWR(2)
      Logical NoPInts, Shijij, AeqB, CeqD, AeqC, ABeqCD,
     &        EQ, Copy, NoCopy, Do_TnsCtl, IJeqKL,IeqK,JeqL,
     &        Pij, Pkl, Pijkl, Pik, Pjl,
     &        lEmpty, Prescreen_On_Int_Only, DoCoul, DoExch,
     &        Scrij, Scrkl, Scrik, Scril, Scrjk, Scrjl,
     &        Batch_On_Disk, DoAOBatch, All_Spherical, NoInts
      Data Copy/.True./, NoCopy/.False./
#include "SysDef.fh"
      External EQ, lEmpty
      Interface
      Integer Function iGet(A,n)
      Integer :: n
      Real*8, Target :: A(*)
      End Function iGet
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *
!     Declaration of statement functions to compute canonical index
!
      Integer ixyz, nabSz
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
!
      Call TwoEl_NoSym_New_Internal(Data1,Data2)
!
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iStb)
         Call Unused_integer(jStb)
         Call Unused_integer(kStb)
         Call Unused_integer(lStb)
         Call Unused_integer(iPrInc)
         Call Unused_integer(kPrInc)
         Call Unused_integer(nData1)
         Call Unused_integer(nData2)
         Call Unused_real_array(FckTmp)
         Call Unused_real_array(SoInt)
         Call Unused_integer(nHRRAB)
         Call Unused_integer(nHRRCD)
         Call Unused_real_array(Aux)
      End If
!
!     This is to allow type punning without an explicit interface
      Contains
      Subroutine TwoEl_NoSym_New_Internal(Data1,Data2)
      Real*8, Target :: Data1(*),Data2(*)
      Integer, Pointer :: iData1(:),iData2(:)
!
      All_Spherical=Shells(iShll(1))%Prjct.and.
     &              Shells(iShll(2))%Prjct.and.
     &              Shells(iShll(3))%Prjct.and.
     &              Shells(iShll(4))%Prjct
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Twoel: Data1',' ',
     &             Data1,nAlpha*nBeta,nDArray-1)
      Call RecPrt('Twoel: Data2',' ',
     &             Data2,nGamma*nDelta,nDArray-1)
      Call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
      Call RecPrt('Coeff2',' ',Coeff2,nBeta,jBasj)
      Call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
      Call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
#endif
!
      RST_triplet=One
      QInd(2)=RST_triplet
      kOp(1)=0
      kOp(2)=0
      kOp(3)=0
      kOp(4)=0
!
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      iSmAng=la+lb+lc+ld
!
!switch (to generate better start orbitals...)
      AeqB = EQ(Coor(1,1),Coor(1,2))
      If (NDDO .AND. .NOT.AeqB) Go To 99
      CeqD = EQ(Coor(1,3),Coor(1,4))
      If (NDDO .AND. .NOT.CeqD) Go To 99
!switch
!
      AeqC = EQ(Coor(1,1),Coor(1,3))
      ABeqCD = AeqB .and. CeqD .and. AeqC
      If (ABeqCD .and. Mod(iSmAng,2).eq.1) Go To 99
!     For Spherical Gaussians, batches like
!     (DS|SS), (FP|SS) and (FS|PS) vanish as well
      If (ABeqCD .and. All_Spherical .and.
     &    2*Max(la,lb,lc,ld).gt.iSmAng) Go To 99
!
      nab = iCmp(1)*iCmp(2)
      ncd = iCmp(3)*iCmp(4)
      nijkl = iBasi*jBasj*kBask*lBasl*nComp
      nabcd = nab*ncd
      nInts =nijkl*nabcd
      ipAOInt=1
      iW3=1+nInts
      iW4=1
!
      vijkl = Data1(ip_abMax(nZeta))
     &      * Data2(ip_abMax( nEta))
!
      Batch_On_Disk = (vijkl.gt.Thize) .and.
     &       (Disc+DBLE(nInts+2+2/RtoI).le.Disc_Mx)
!
      Prescreen_On_Int_Only = IntOnly
      If (DoIntegrals) Prescreen_On_Int_Only=.True.
      If (Batch_On_Disk) Prescreen_On_Int_Only = .True.
!
      If (DoFock) Then
         vij = Dij(mDij,1)
         vkl = Dkl(mDkl,1)
!        Coulomb contributions
         Scrkl = vij*vijkl.ge.ThrInt
         Scrij = vkl*vijkl.ge.ThrInt
         DoCoul  = Scrij .or. Scrkl
!
!        Exchange contributions
         vik = Dik(mDik,1)/Four
         vil = Dil(mDil,1)/Four
         vjk = Djk(mDjk,1)/Four
         vjl = Djl(mDjl,1)/Four
         Scrjl = vik*vijkl.ge.ThrInt
         Scrjk = vil*vijkl.ge.ThrInt
         Scril = vjk*vijkl.ge.ThrInt
         Scrik = vjl*vijkl.ge.ThrInt
         DoExch = Scrjl .or. Scrjk .or. Scril .or. Scrik
         If (FckNoClmb) DoCoul=.False.
         If (FckNoExch) DoExch=.False.
      Else
         DoCoul=.False.
         DoExch=.False.
      End If
      DoAOBatch=(DoIntegrals.and.vijkl.gt.CutInt).or.
     &          (DoFock.and.(DoCoul.or.DoExch)) .or.
     &          (Batch_On_Disk.and.W2Disc)
!
!-----Branch out if crude estimate indicates no contributions!
!
      If (.Not.DoAOBatch) Then
         If (.Not.Batch_On_Disk) Then
            Go To 99
         Else If (Batch_On_Disk.and..Not.W2Disc) Then
 1111       Continue
            Call iRBuf(iWR,2,Copy)
            Call dRBuf(QInd,2,Copy)
            Call Store_QLast(QInd)
            kInts=iWR(1)
            mInts=iWR(2)
            If (QInd(1).eq.Quad_ijkl) Then
               If (kInts.ne.nInts) Then
                  Call WarningMessage(2,
     &                        'Twoel: kInts.ne.nInts!')
                  Write (6,*) 'Twoel: kInts,mInts,nInts=',
     &                                kInts,mInts,nInts
                  Write (6,*) 'Index,1:',QInd(1),Quad_ijkl
                  Call Abend()
               End If
               If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,NoCopy)
               Disc = Disc + DBLE(2/RtoI + 2 + mInts)
               Go To 99
            Else If (QInd(1).lt.Quad_ijkl) Then
               If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,NoCopy)
               Disc = Disc + DBLE(2/RtoI + 2 + mInts)
               Go To 1111
            Else
               Call WarningMessage(2,'Twoel: batch is lost!')
               Write (6,*) 'Index,1:',QInd(1),QInd(2),
     &                      Quad_ijkl,RST_triplet
               Call Abend()
            End If
         End If
      End If
!
!
!-----Branch point for partial integral storage
!
      If (Batch_On_Disk.and..Not.W2Disc) Go To 6767
!                                                                      *
!***********************************************************************
!                                                                      *
!     Here if the AO batch will be computed!
!
!-----Compute actual size of the {a0|c0} block
!
      mabMin=nabSz(Max(la,lb)-1)+1
      If (EQ(Coor(1,1),Coor(1,2))) mabMin = nabSz(la+lb-1)+1
      mabMax=nabSz(la+lb)
      mcdMin=nabSz(Max(lc,ld)-1)+1
      If (EQ(Coor(1,3),Coor(1,4))) mcdMin = nabSz(lc+ld-1)+1
      mcdMax=nabSz(lc+ld)
      mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
!
!-----Find the proper centers to start of with the angular
!     momentum on. If la.eq.lb there will exist an
!     ambiguity to which center that angular momentum should
!     be accumulated on. In that case we will use A and C of
!     the order as defined by the basis functions types.
!
      If (iAnga(1).ge.iAnga(2)) Then
         call dcopy_(3,Coor(1,1),1,CoorAC(1,1),1)
      Else
         call dcopy_(3,Coor(1,2),1,CoorAC(1,1),1)
      End If
      If (iAnga(3).ge.iAnga(4)) Then
         call dcopy_(3,Coor(1,3),1,CoorAC(1,2),1)
      Else
         call dcopy_(3,Coor(1,4),1,CoorAC(1,2),1)
      End If
!
!-----Set flags if triangularization will be used
!
      IeqK = EQ(Coor(1,1),Coor(1,3))
      JeqL = EQ(Coor(1,2),Coor(1,4))
      IJeqKL = IeqK .and. JeqL
!
!-----Loops to partion the primitives
!
      IncZet=nAlpha*jPrInc
      IncEta=nGamma*lPrInc
      If (nZeta.ne.IncZet.or.nEta.ne.IncEta) Then
         mWork2 = nWork2 - nijkl*mabcd
         ipAOInt=1+nijkl*mabcd
      Else
         mWork2 = nWork2
         ipAOInt=1
      End If
!
      nZeta_Tot=iGet(Data1(ip_IndZ(1,nZeta)),nZeta+1)
      nEta_Tot =iGet(Data2(ip_IndZ(1,nEta )),nEta +1)
#ifdef _DEBUGPRINT_
      Write (6,*) 'nZeta_Tot, IncZet=',nZeta_Tot, IncZet
      Write (6,*) 'nEta_Tot,  IncEta=',nEta_Tot,  IncEta
#endif
!
      kabcd=0
      Do_TnsCtl=.False.
      NoInts=.True.
      NoPInts = .True.
      ipAOInt_=ipAOInt
      iW4_=iW4
!
      Do iZeta = 1, nZeta_Tot, IncZet
         mZeta=Min(IncZet,nZeta_Tot-iZeta+1)
         If (lEmpty(Coeff2,nBeta,nBeta,jBasj)) Cycle
!
         Do iEta  = 1, nEta_Tot,  IncEta
            mEta=Min(IncEta,nEta_Tot-iEta+1)
            If (lEmpty(Coeff4,nDelta,nDelta,lBasl)) Cycle
!
            Call C_F_Pointer(C_Loc(Data1(ip_IndZ(1,nZeta))),
     &                             iData1,[nZeta])
            Call C_F_Pointer(C_Loc(Data2(ip_IndZ(1,nEta))),
     &                             iData2,[nEta])
            Call DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,
     &                  nZeta_Tot,nEta_Tot,
     &                  Data1,mData1,
     &                  Data2,mData2,
     &                  nAlpha,nBeta,nGamma,nDelta,
     &                  iData1,
     &                  Zeta,ZInv,P,KappAB,IndZet,
     &                  iData2,
     &                  Eta,EInv,Q,KappCD,IndEta,
     &                  1,1,1,1,1,1,ThrInt,CutInt,
     &                  vij,vkl,vik,vil,vjk,vjl,
     &                  Prescreen_On_Int_Only,NoInts,iAnga,
     &                  Coor,CoorAC,
     &                  mabMin,mabMax,mcdMin,mcdMax,nijkl/nComp,
     &                  nabcd,mabcd,Wrk,ipAOInt_,iW4_,
     &                  nWork2,mWork2,
     &                  Data1(ip_HrrMtrx(nZeta)),
     &                  Data2(ip_HrrMtrx(nEta) ),
     &                  la,lb,lc,ld,
     &                  iCmp,iShll,NoPInts,
     &                  Dij(1,1),mDij,Dkl(1,1),mDkl,Do_TnsCtl,kabcd,
     &                  Coeff1,iBasi,Coeff2,jBasj,
     &                  Coeff3,kBask,Coeff4,lBasl)
            Nullify(iData1,iData2)
!
         End Do
      End Do
!
      If (NoPInts) Then
         If (W2Disc) Then
            If (Batch_On_Disk) Then
               iOpt=0
               mInts=0
!
               iWR(1)=nInts
               iWR(2)=mInts
               Call iWBuf(iWR,2)
               QInd(1)=Quad_ijkl
               Call dWBuf(QInd,2)
               Call Store_QLast(QInd)
!
               Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            End If
         End If
         Go To 99
      End If
!
!-----Apply the transfer equation and transform the spherical
!     harmonic gaussian.
!
      If (Do_TnsCtl) Then
         Call TnsCtl(Wrk(iW4),nWork2,Coor,
     &               nijkl,mabMax,mabMin,mcdMax,mcdMin,
     &               Data1(ip_HrrMtrx(nZeta)),
     &               Data2(ip_HrrMtrx(nEta) ),
     &               la,lb,lc,ld,
     &               iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &               iShll(1),iShll(2),iShll(3),iShll(4),i_Int)
         ipAOInt=i_Int
         If (i_Int.eq.1) Then
            iW3=1+nijkl*nabcd
         Else
            iW3=1
         End If
      Else
!
!--------Undo the late Cntrct
!
         call dcopy_(nijkl*nabcd,Wrk(ipAOInt),1,Wrk(iW3),1)
         Call DGeTMO(Wrk(iW3),nabcd,nabcd,nijkl,Wrk(ipAOInt),nijkl)
!
      End If
#ifdef _DEBUGPRINT_
      Call RecPrt('(AB|CD)',' ',Wrk(ipAOInt),nijkl/nComp,
     &            nComp*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4))
#endif
!
!-----Branch point for partial integral storage
!
      If (Batch_On_Disk.and.W2Disc) Then
!
!--------Write integrals to current position on disc.
!
         iOpt=0
         Call PkR8(iOpt,nInts,nByte,Wrk(ipAOInt),Wrk(iW3))
         mInts=(nByte+RtoB-1)/RtoB
!
         iWR(1)=nInts
         iWR(2)=mInts
!        Write (*,*) 'nInts,mInts=',nInts,mInts
         Call iWBuf(iWR,2)
         QInd(1)=Quad_ijkl
         Call dWBuf(QInd,2)
         Call Store_QLast(QInd)
         Call dWBuf(Wrk(iW3),mInts)
!
         Disc = Disc + DBLE(2/RtoI + 2 + mInts)
!
      End If
!
 6767 Continue
      If (Batch_On_Disk.and..Not.W2Disc) Then
 1112    Continue
         Call iRBuf(iWR,2,Copy)
         Call dRBuf(QInd,2,Copy)
         Call Store_QLast(QInd)
         kInts=iWR(1)
         mInts=iWR(2)
         If (QInd(1).eq.Quad_ijkl) Then
            If (kInts.ne.nInts) Then
               Call WarningMessage(2,
     &                     'Twoel: kInts.ne.nInts!')
               Write (6,*) 'Twoel: kInts,mInts,nInts=',
     &                             kInts,mInts,nInts
               Write (6,*) 'Index,1:',QInd(1),Quad_ijkl
               Call Abend()
            End If
            If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,Copy)
            Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            If (mInts.eq.0) Go To 99
         Else If (QInd(1).lt.Quad_ijkl) Then
            If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,NoCopy)
            Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            Go To 1112
         Else
            Call WarningMessage(2,'Twoel: batch is lost!')
            Write (6,*) 'Index,1:',QInd(1),QInd(2),
     &                   Quad_ijkl,RST_triplet
            Call Abend()
         End If
!
         iOpt=0
         Call UpkR8(iOpt,nInts,nByte,Wrk(iW3),Wrk(ipAOInt))
      End If
!
!-----Accumulate contributions directly to the Fock matrix.
!
      If (DoFock)
     &Call FckAcc_NoSymq(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                  Shijij, iShll, iShell, nijkl,
     &                  Wrk(ipAOInt),TwoHam,Dens,Size(TwoHam),
     &                  iAO,iAOst,
     &                  iBasi,jBasj,kBask,lBasl,DoCoul,DoExch,
     &                  vij,vkl,vik,vil,vjk,vjl,ExFac)
!
      If (DoIntegrals) Then
         If (ipAOInt.ne.1) Then
            call dcopy_(nijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),
     &                  Wrk(ipAOInt),1,Wrk(1),1)
            ipAOInt=1
         ENd If
         iPer = 1
         Pij= iS_.eq.jS_
         Pkl= kS_.eq.lS_
         Pik= iS_.eq.kS_
         Pjl= jS_.eq.lS_
         Pijkl= Pij. and. Pkl .and. Pik .and. Pjl
         If (Pij)   iPer = iPer*2
         If (Pkl)   iPer = iPer*2
         If (Pijkl) iPer = iPer*2
         q4 = DBLE(8)/DBLE(iPer)
         If (nIrrep.eq.1) q4 = One
         If (q4.ne.One) Call DScal_(nijkl*iCmp(1)*iCmp(2)
     &                            *iCmp(3)*iCmp(4),q4,Wrk(ipAOInt),1)
      End If
  99  Continue
      End Subroutine TwoEl_NoSym_New_Internal
!
      End
