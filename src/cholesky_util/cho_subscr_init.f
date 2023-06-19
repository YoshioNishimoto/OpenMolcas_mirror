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
      SubRoutine Cho_SubScr_Init()
!
!     Purpose: initialize screening in vector subtraction.
!
      use ChoSubScr, only: DSubScr, DSPNm
      use stdalloc
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"
      Integer l_DSubScr

      l_DSubScr = nnBstR(1,1)
      Do iSym = 2,nSym
         l_DSubScr = max(l_DSubScr,nnBstR(iSym,1))
      End Do
      Call mma_allocate(DSubScr,l_DSubScr,Label='DSubScr')

      Call mma_allocate(DSPNm,nnShl,Label='DSPNm')

      End
