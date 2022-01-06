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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine Overlap(mGrid,iSpin)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      use nq_Grid, only: Rho, vRho
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8, Parameter:: T_x=1.0D-20
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
      Rho_Min=T_X*1.0D-2
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
         d_alpha=Rho(1,iGrid)
         DTot=Two*d_alpha
         If (DTot.lt.T_X) Go To 199
*
*------- Accumulate contributions to the integrated density
*
         F_xc(iGrid)=F_xc(iGrid)+Dtot
*
         vRho(1,iGrid)=One
*
 199     Continue
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/=1
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
         d_alpha=Max(Rho_min,Rho(1,iGrid))
         d_beta =Max(Rho_min,Rho(2,iGrid))
         DTot=d_alpha+d_beta
         If (DTot.lt.T_X) Go To 299
*
*------- Accumulate contributions to the integrated density
*
         F_xc(iGrid)=F_xc(iGrid)+Dtot
*
         vRho(1,iGrid)=One
         vRho(2,iGrid)=One
*
 299     Continue
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      Return
      End
