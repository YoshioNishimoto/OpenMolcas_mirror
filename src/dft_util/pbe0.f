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
* Copyright (C) 2005, Per Ake Malmqvist                                *
************************************************************************
      Subroutine PBE0(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object: To compute the PBE0 density functional and its first         *
* derivatives. This functional is defined as a hybrid functional       *
* having 25% HF exchange, 75% PBE exchange, and 100% PBE correlation.  *
*   See   J.P. Perdew, M. Ernzerhof, and K. Burke,                     *
*      J. Chem. Phys. 105, 9982 (1996).                                *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. December 2005                *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *
      CoeffA=1.0D0*CoefR
      Call CPBE(mGrid,CoeffA,iSpin,F_xc)

      CoeffB=0.75D0*CoefX
      Call XPBE(mGrid,CoeffB,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
