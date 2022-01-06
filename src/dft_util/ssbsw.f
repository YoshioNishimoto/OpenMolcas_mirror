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
*               2010, Grigory A. Shamov                                *
*               2017, Giovanni Li Manni                                *
*               2017, Aron Cohen                                       *
************************************************************************
      Subroutine SSBSW(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object:     Combination of excnange SSB-sw and PBE correlation terms *
*             ref (secondary): Swart, Sola, Bickelhaupt                *
*             J.Chem.Phys. 131 (2009) 094103. Note that it isnt SSB-D  *
*                                                                      *
*      Author: G. Li Manni & A. Cohen, Max Planck Institute Stuttgart  *
*              Summer 2017, edited in Cambridge (UK) & Palermo (Sicily)*
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *
*---- SSBSW Exchange -- unlike OPTX, SSBSW has its LDA part included !
      Coeff=1.0d0*CoefX
      Call xSSBSW(mGrid,Coeff,iSpin,F_xc)
*
*---- PBE Correlation
      Coeff=1.0d0*CoefR
      Call CPBE(mGrid,Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
