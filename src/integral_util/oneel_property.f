!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine OneEl_Property(Kernel,KrnlMm,Label,ip,lOper,nComp,
     &                          CCoor,nOrdOp,rNuc,rHrmt,iChO,
     &                          D_tot,nDens,Property,Sig)
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      use Integral_Interfaces, only: int_kernel, int_mem,
     &                               OneEl_Integrals
      use Constants, only: One
      use stdalloc, only: mma_deallocate
      Implicit None
      Procedure(int_kernel) :: Kernel
      Procedure(int_mem) :: KrnlMm
      Character(LEN=8) Label
      Integer nComp, nDens, nOrdOp
      Real*8 CCoor(3,nComp), rNuc(nComp), Property(nComp), D_tot(nDens)
      Integer ip(nComp), lOper(nComp), iChO(nComp)
      Real*8 rHrmt, Sig

      Real*8, Allocatable:: Integrals(:)
      Integer, External:: n2Tri
      Integer LenTot, iComp, iSmLbl, nInt
      Real*8, external:: DDot_
!                                                                      *
!***********************************************************************
!                                                                      *
      If (rHrmt.ne.One) Then
         Call WarningMessage(2,'OneEl_Property: rHrmt.ne.One')
         Call Abend()
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Compute the one-electron integrals
!
      Call OneEl_Integrals(Kernel,KrnlMm,Label,ip,lOper,nComp,
     &                     CCoor,nOrdOp,rHrmt,iChO,Integrals)
!                                                                      *
!***********************************************************************
!                                                                      *
!                    P O S T P R O C E S S I N G                       *
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
      Call PrMtrx(Label,lOper,nComp,ip,Integrals)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Compute properties
!
      LenTot=0
      Do iComp = 1, nComp
         iSmLbl = lOper(iComp)
!                                                                      *
!***********************************************************************
!                                                                      *
!--------Compute properties directly from integrals
!
         nInt=n2Tri(iSmLbl)
         LenTot = LenTot + nInt + 4
         If (nInt.ne.0) Then
            Call CmpInt(Integrals(ip(iComp)),nInt,nBas,nIrrep,iSmLbl)
            If (nInt.ne.nDens) Then
               Call WarningMessage(2,'OneEl_Property: nInt.ne.nDens')
               Write (6,*) 'nInt=',nInt
               Write (6,*) 'nDens',nDens
               Call Abend()
            End If
            Property(iComp)=rNuc(iComp)
     &                     -Sig*DDot_(nDens,D_Tot,1,
     &                                      Integrals(ip(iComp)),1)
         Else
            Property(iComp)=rNuc(iComp)
         End If
!
      End Do  ! iComp
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Deallocate memory for integral
!
      Call mma_deallocate(Integrals)
!
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End SubRoutine OneEl_Property
