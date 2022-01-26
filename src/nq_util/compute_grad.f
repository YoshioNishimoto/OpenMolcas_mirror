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
      Real*8 Function Compute_Grad(Weights,mGrid,iSpin,T_X)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use nq_Grid, only: Rho, Sigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Weights(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*
      Compute_Grad=Zero
      Rho_min=T_X*1.0D-2
*
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
         d_alpha=Rho(1,iGrid)
         DTot=d_alpha
         If (DTot.lt.T_X) Go To 199
         Gamma=Sqrt(Sigma(1,iGrid))
*
*------- Accumulate contributions to the integrated Tau
*
         Compute_Grad = Compute_Grad + Two*Gamma*Weights(iGrid)
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
         Gamma=Sqrt(Sigma(1,iGrid)+Two*Sigma(2,iGrid)+Sigma(3,iGrid))
*
*------- Accumulate contributions to the integrated density
*
         Compute_Grad = Compute_Grad + Gamma*Weights(iGrid)
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
