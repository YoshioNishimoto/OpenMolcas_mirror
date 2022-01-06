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
*               2010, Grigory A. Shamov                                *
************************************************************************
      Subroutine RGE2(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object: To compute the sum of the RGE2 exchange functional           *
*         and the PBEsol correlation. The combination is RGE2          *
*         as described in Ruzsinsky et al, JCTC 2009, 5, 763           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author: template by Per Ake Malmqvist,                          *
*             University of Lund, SWEDEN. December 2005                *
*              Grigory A Shamov. U of Manitoba, Feb 2010               *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *
      CoeffA=One*CoefR
      Call CPBEsol(mGrid,CoeffA,iSpin,F_xc)

      CoeffB=One*CoefX
      Call XRGE2(mGrid,CoeffB,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
