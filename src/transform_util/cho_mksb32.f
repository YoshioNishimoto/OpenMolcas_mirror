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
* Copyright (C) 2005, Giovanni Ghigo                                   *
************************************************************************
      Subroutine MkExSB32(iAddSB,LenSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV, iAddSBt,LenSBt)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(3,2) (p,q secondary,active) of  *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical SameLx

*   - SubBlock 3 2
      LenSB = nSsh(iSymA) * nAsh(iSymB)
      Call GetMem('SB','Allo','Real',iAddSB,LenSB)
      If (iSymA.EQ.iSymB .and. iSymI.EQ.iSymJ .and. iI.EQ.iJ) then
*       SB 3,2 = (SB 2,3)+
        Call Trnsps(nSsh(iSymB),nAsh(iSymA),Work(iAddSBt),Work(iAddSB))
        Return
      EndIf

*     Build Lx
      Call GetMem('Lx','Allo','Real',iAddLx0,nSsh(iSymA)*numV)
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL3(iSymA,iSymI,iI,numV, LxType,iIx, iAddLx0,SameLx)

*     Build Ly
      Call GetMem('Ly','Allo','Real',iAddLy0,nAsh(iSymB)*numV)
      Call MkL2(iSymB,iSymJ,iJ,numV, LxType,iIx, iAddLy0,SameLx)

*     Generate the SubBlock
      Call DGEMM_('N','T',nAsh(iSymB),nSsh(iSymA),numV,1.0d0,
     &    Work(iAddLy0),nAsh(iSymB), Work(iAddLx0),nSsh(iSymA),
     &                     0.0d0,Work(iAddSB),nAsh(iSymB) )

      Call GetMem('Ly','Free','Real',iAddLy0,nAsh(iSymB)*numV)
      Call GetMem('Lx','Free','Real',iAddLx0,nSsh(iSymA)*numV)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(LenSBt)
      End

      Subroutine MkCouSB32(AddSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           July 2005                                                  *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(3,2) (p secondary, q active) of *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"

      Real*8, Allocatable:: AddSB(:)

*   - SubBlock 3 2
      LenSB = nSsh(iSymA) * nAsh(iSymB)
      Call mma_allocate(AddSB,LenSB,Label='AddSB')
      Call GetMem('SBt','Allo','Real',iAddSBt,LenSB)

*     Define Lab
      iAddAB = iMemTCVX(5,iSymA,iSymB,1)
      LenAB  = LenSB
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)'     MkCouSB32: TCVE(',iSymA,',',iSymB,')'
c      Write(6,'(8F10.6)')(Work(iAddAB+k),k=0,LenAB*numV-1)
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------

*     Build Lij
      LenLij = numV
      Call GetMem('Lij','Allo','Real',iAddLij,LenLij)
      Call MkLij(iSymI,iSymJ,iI,iJ,numV, iAddLij)

*     Generate the SubBlock
      Call DGEMM_('N','N',LenAB,1,numV,1.0d0,
     &    Work(iAddAB),LenAB, Work(iAddLij),LenLij,
     &                0.0d0,Work(iAddSBt),LenSB )
      Call Trnsps(nSsh(iSymA),nAsh(iSymB),Work(iAddSBt),AddSB)

      Call GetMem('Lij','Free','Real',iAddLij,LenLij)
      Call GetMem('SBt','Free','Real',iAddSBt,LenSB)

      Return
      End
