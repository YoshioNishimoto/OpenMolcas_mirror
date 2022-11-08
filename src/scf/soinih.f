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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2017,2022, Roland Lindh                                *
************************************************************************
      SubRoutine SOiniH()
************************************************************************
*                                                                      *
*     purpose: generate initial Hessian (diagonal) from                *
*              orbital energies (for second order update)              *
*                                                                      *
*     called from: WfCtl                                               *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use Orb_Type, only: OrbType
      use InfSCF, only: nSym, nFro, nOrb, nOcc
      use SCF_Arrays, only: EOrb, HDiag, CMO_Ref
      use Constants, only: Zero, Four
      Implicit None
*
*     declaration local variables
      Integer iD, nD
      Integer iSym,ii,ia,ioffs,iHoffs,nOccmF,nOrbmF
      Real*8, Parameter:: Hii_Min=0.05D0
      Real*8, Parameter:: Hii_Max=1.00D0
*     Real*8 :: Hii
*     Integer i
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*     Set the array to silly large values. In the case of UHF these
*     will remain but should not make any difference. They are actully
*     needed to make the rs-rfo code work.
*

*     Compute the diagonal values of the Fock matrix, stored in EOrb.
      Call Mk_EOrb(CMO_Ref,Size(CMO_Ref,1),Size(CMO_Ref,2))

      nD   =Size(EOrb,2)
      HDiag(:)=1.0D+99
*
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'nD=',nD
      Write (6,*) 'kOV(:)=',kOV(:)
      Do iD = 1, nD
         Write (6,*) 'iD=',iD
         Write (6,'(A,8I3)') 'nOcc',(nOcc(iSym,iD),iSym=1,nSym)
      End Do
      Write (6,'(A,8I3)') 'nOrb',(nOrb(iSym),iSym=1,nSym)
#endif
      iHoffs=1
      Do iD = 1, nD
*
         ioffs=1
         Do iSym=1,nSym
*
*            loop over all occ orbitals in sym block
*
             ioffs=ioffs+nFro(iSym)
             nOccmF=nOcc(iSym,iD)-nFro(iSym)
             nOrbmF=nOrb(iSym)-nFro(iSym)
*
#ifdef _DEBUGPRINT_
             iHoffs_ = iHoffs
#endif
             Do ii=ioffs,ioffs+nOccmF-1
*
*               loop over all virt orbitals in sym block
*
                Do ia=ioffs+nOccmF,ioffs+nOrbmF-1
*
                   If (OrbType(ia,iD).eq.OrbType(ii,iD))
     &             HDiag(iHoffs)=Four*(EOrb(ia,iD)-EOrb(ii,iD))
     &                             /DBLE(nD)
                   If (HDiag(iHoffs)<Zero) Then
                       Write (6,*) 'Hii<0.0, Hii=',HDiag(iHoffs)
                       HDiag(iHoffs)=Max(Hii_Max,Abs(HDiag(iHoffs)))
                   Else If (Abs(HDiag(iHoffs)).lt.Hii_Min) Then
                       HDiag(iHoffs)=Hii_Min
                       Write (6,*) 'Abs(Hii)<0.05'
                   End If
*
                   iHoffs=iHoffs+1
*
                End Do  ! ia
*
             End Do     ! ii
*
#ifdef _DEBUGPRINT_
             Write (6,*) 'nOccmF,nOrbmF=',nOccmF,nOrbmF
             If ((nOrbmF-nOccmF)*nOccmF.gt.0)
     &          Call RecPrt('HDiag',' ',HDiag(iHoffs_),
     &                      nOrbmF-nOccmF,nOccmF)
#endif
*
             ioffs=ioffs+nOrbmF
*
          End Do ! iSym
      End Do ! iD

#ifdef _DEBUGPRINT_
      Call RecPrt('HDiag',' ',HDiag(:),1,Size(HDiag))
#endif
      Return
      End
