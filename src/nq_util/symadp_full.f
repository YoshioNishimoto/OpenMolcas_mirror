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
* Copyright (C) 1991,2021, Roland Lindh                                *
************************************************************************
      Subroutine SymAdp_Full(AOIntegrals, nBfn, PrpInt, nPrp,
     &                       list_s,nlist_s,Fact,ndc)
************************************************************************
*                                                                      *
* Object: to transform the one-electon matrix elements from AO basis   *
*         to SO basis.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1991                                             *
************************************************************************
      use iSD_data
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info,     only: iAOtSO
      use nq_Grid,       only: iBfn_Index
      use Basis_Info,    only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOIntegrals(nBfn,nBfn), PrpInt(nPrp), Fact(ndc,ndc)
      Integer list_s(2,nlist_s)
      Integer nOp(2)
      Integer, Parameter:: iTwoj(0:7)=[1,2,4,8,16,32,64,128]
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*     Call RecPrt('AOIntegrals',' ',AOIntegrals,nBfn,nBfn)
      loper=1
      Do j1 = 0, nIrrep-1
         Do iBfn = 1, nBfn
            ilist_s = iBfn_Index(2,iBfn)
            iCmp    = iBfn_Index(3,iBfn)
            indAO1  = iBfn_Index(6,iBfn)
            iSkal   = list_s(1,ilist_s)
            kDCRE   = list_s(2,ilist_s)
            iAO     = iSD( 7,iSkal)
            mdci    = iSD(10,iSkal)
            iShell  = iSD(11,iSkal)
            nOp(1) = NrOpr(kDCRE)
            xa = DBLE(iChTbl(j1,nOp(1)))
            If (iAOtSO(iAO+iCmp,j1)<0) Cycle
            iSO1=iAOtSO(iAO+iCmp,j1)

            Do j2 = 0, nIrrep-1
               j12 = iEor(j1,j2)
               If (iAnd(lOper,iTwoj(j12)).eq.0) Cycle

               iPnt = iPntSO(j1,j2,lOper,nbas)

               Do jBfn = 1, nBfn
                  jlist_s = iBfn_Index(2,jBfn)
                  If (jlist_s<ilist_s) Cycle
                  jCmp    = iBfn_Index(3,jBfn)
                  indAO2  = iBfn_Index(6,jBfn)
                  jSkal   = list_s(1,jlist_s)
                  kDCRR   = list_s(2,jlist_s)
                  jAO     = iSD( 7,jSkal)
                  mdcj    = iSD(10,jSkal)
                  jShell  = iSD(11,jSkal)
                  nOp(2) = NrOpr(kDCRR)
                  xb = DBLE(iChTbl(j2,nOp(2)))
                  If (iAOtSO(jAO+jCmp,j2)<0) Cycle
                  iSO2=iAOtSO(jAO+jCmp,j2)

                  If (iShell.eq.jShell .and. iCmp<jCmp .and.
     &                nOp(1).eq.nOp(2)) Cycle
*
                  iSO=iSO1+IndAO1-1
                  jSO=iSO2+IndAO2-1

*------------     Diagonal symmetry block
                  If (iShell.ne.jShell .and. iSO1.eq.iSO2 .and.
     &                iSO<jSO) Cycle
                  If (iShell.eq.jShell .and. iSO1.eq.iSO2 .and.
     &                nOp(1).eq.nOp(2) .and. iSO<jSO) Cycle
                  xaxb=xa*xb
                  If (iShell.eq.jShell .and. iSO1.eq.iSO2 .and.
     &                nOp(1).ne.nOp(2) .and. iSO==jSO) xaxb=xaxb*Two
                  Indij=iPnt + iTri(iSO,jSO)

*                 If (Indij==9) Then
*                 Write (*,*) 'Indij=',Indij
*                 Write (*,*) 'iShell,jShell=',iShell,jShell
*                 Write (*,*) 'iSO1,iSO2=',iSO1,iSO2
*                 Write (*,*) 'iSO,jSO=',iSO,jSO
*                 Write (*,*) 'iBfn,jBfn=',iBfn,jBfn
*                 Write (*,*) 'nOp(1),nOp(2)=',nOp(1),nOp(2)
*                 Write (*,*) 'Fact(mdci,mdcj),xaxb=',Fact(mdci,mdcj),
*    &                                                xaxb
*                 Write (*,*) 'AOIntegrals(iBfn,jBfn)=',
*    &                         AOIntegrals(iBfn,jBfn)
*                 End If
                  PrpInt(Indij) = PrpInt(Indij)
     &                          + Fact(mdci,mdcj)*xaxb
     &                          * AOIntegrals(iBfn,jBfn)
*
               End Do ! jBfn
            End Do    ! j2

         End Do       ! iBfn
      End Do          ! j1
*     Call RecPrt('PrpInt',' ',PrpInt,1,nPrp)
*
      Return
      End
