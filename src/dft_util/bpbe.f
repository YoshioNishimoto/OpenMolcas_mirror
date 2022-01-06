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
* Copyright (C) 2001, Roland Lindh                                     *
************************************************************************
      Subroutine BPBE(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object: B88 + PBE combination. Resembles popular BP86                *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. March 2001                              *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange
*
      Coeff=One*CoefX
      Call Diracx(mGrid,iSpin,F_xc,Coeff)
*
*---- Becke 88 Exchange
*
      Coeff=One*CoefX
      Call xB88(mGrid,
     &          Coeff,iSpin,F_xc)
*
*---- Lee-Yang-Parr Correlation
*
      Coeff=One*CoefR
      Call CPBE(mGrid,
     &         Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
