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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine PLF2(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
!***********************************************************************
!                                                                      *
!  object: to sift and index the petite list format integrals.         *
!                                                                      *
!          the indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          May '90                                                     *
!***********************************************************************
      use SOAO_Info, only: iAOtSO
      use k2_arrays, only: Sew_Scr
      use lw_Info, only: lwSyB, lwInt, lwSqn
      use Gateway_Info, only: ThrInt
      use sort_data, only: DimSyB, lSll
      use Constants, only: One
      Implicit None
!
      Integer ijkl,iCmp,jCmp,kCmp,lCmp,iBas,jBas,kBas,lBas
      Real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4), iSOs(4)

      Integer i, j, iTri, nUt, iAOSti, iAOStj, iAOStk, iAOStl,
     &        iAOi, iAOj, iAOk, iAOl, i1, i2, i3, i4, nij, mij,
     &        iSO, jSO, kSO, lSO, iSOi, jSOj, kSOk, lSOl, nijkl,
     &        iSOij, iSOkl, ijklCmp, iBin
      Real*8 AInt
#ifdef _DEBUGPRINT_
      Real*8 r1, r2
      Real*8, External:: DDot_
#endif
!
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!
#ifdef _DEBUGPRINT_
      r1=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
      r2=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
      Write (6,*) ' Sum=',r1
      Write (6,*) ' Dot=',r2
      Call RecPrt(' In Plf2: AOInt',' ',
     &                              AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
#endif
!
!     Allocate space to store integrals together with their
!     Symmetry batch and sequence number.
!     To avoid conflicts in using memory this is done in the
!     subroutine PSOAO
!
      nUt=-1
      nij=DimSyB(1,1)
      mij=lSll(1)/nij
!     Write (*,*) 'nij,mij=',nij,mij
!
!     quadruple loop over elements of the basis functions angular
!     description. loops are reduced to just produce unique SO integrals
!     observe that we will walk through the memory in AOint in a
!     sequential way.
!
      iAOsti=iAOst(1)
      iAOstj=iAOst(2)
      iAOstk=iAOst(3)
      iAOstl=iAOst(4)
      iAOi=iAO(1)
      iAOj=iAO(2)
      iAOk=iAO(3)
      iAOl=iAO(4)
!
      ijklCmp=iCmp*jCmp*kCmp*lCmp
      Call DCopy_(ijkl*2*ijklCmp,[One],0,Sew_Scr(lwSyB),1)
!
      Do 100 i1 = 1, iCmp
         iSOs(1)=iAOtSO(iAOi+i1,kOp(1))+iAOsti
         Do 200 i2 = 1, jCmp
            iSOs(2)=iAOtSO(iAOj+i2,kOp(2))+iAOstj
            Do 300 i3 = 1, kCmp
               iSOs(3)=iAOtSO(iAOk+i3,kOp(3))+iAOstk
               Do 400 i4 = 1, lCmp
                  iSOs(4)=iAOtSO(iAOl+i4,kOp(4))+iAOstl
!
                iSO =iSOs(1)
                jSO =iSOs(2)
                kSO =iSOs(3)
                lSO =iSOs(4)
!
                nijkl = 0
                Do 120 lSOl = lSO, lSO+lBas-1
                   Do 220 kSOk = kSO, kSO+kBas-1
                      iSOkl = iTri(kSOk,lSOl)
                      Do 320 jSOj = jSO, jSO+jBas-1
                         Do 420 iSOi = iSO, iSO+iBas-1
                            nijkl = nijkl + 1
                            AInt=AOint(nijkl,i1,i2,i3,i4)
                            If (Abs(AInt).lt.ThrInt) Go To 420
                            iSOij = iTri(iSOi,jSOj)
!
!                           Write (*,*) 'iSOij,iSOkl=',iSOij,iSOkl
!
                            nUt=nUt+1
                            Sew_Scr(lwInt+nUt)=Aint
                            iBin=(iSOkl-1)/mij
!                           Write (*,*) 'iBin=',iBin+1
                            Sew_Scr(lwSyB+nUt)=DBLE(iBin+1)
                            Sew_Scr(lwSqN+nUt)=DBLE((iSOkl-1-iBin*mij)
     &                                      *nij+iSOij)
!                           Write (*,*) 'iSq=',Sew_Scr(lwSqN+nUt)
!
                            If (iSOij.ne.iSOkl) Then
                               nUt=nUt+1
                               Sew_Scr(lwInt+nUt)=Aint
                               iBin=(iSOij-1)/mij
!                              Write (*,*) 'iBin=',iBin+1
                               Sew_Scr(lwSyB+nUt)=DBLE(iBin+1)
                               Sew_Scr(lwSqN+nUt)=
     &                               DBLE((iSOij-1-iBin*mij)*nij+iSOkl)
!                              Write (*,*) 'iSq=',Sew_Scr(lwSqN+nUt)
                            End If
!
420                      Continue
320                   Continue
220                Continue
120             Continue
!
400            Continue
300         Continue
200      Continue
100   Continue
!
!     pass the integral to phase 1 of the bin sorting algorithm
!
      Call SORT1A(nUt+1,Sew_Scr(lwInt),Sew_Scr(lwSqN),Sew_Scr(lwSyB))
      nUt=0
!
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
      End If
      End Subroutine PLF2
