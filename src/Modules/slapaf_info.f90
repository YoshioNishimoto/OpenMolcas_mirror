!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
Module Slapaf_Info
implicit none
Private
Public:: Cx, Gx, Gx0, NAC, Free_Slapaf
Real*8, Allocatable:: Cx(:,:,:)     ! list of Cartesian coordinates
Real*8, Allocatable:: Gx(:,:,:)     ! list of Cartesian Gradients, State 1
Real*8, Allocatable:: Gx0(:,:,:)    ! list of Cartesian Gradients, State 2 for optimization of conical intersections
Real*8, Allocatable:: NAC(:,:)      ! list of Cartesian non-adiabatic coupling vector

Contains
  Subroutine Free_Slapaf()
#include "stdalloc.fh"
  If (Allocated(Cx)) Call mma_deallocate(Cx)
  If (Allocated(Gx)) Call mma_deallocate(Gx)
  If (Allocated(Gx0)) Call mma_deallocate(Gx0)
  If (Allocated(NAC)) Call mma_deallocate(NAC)
  End Subroutine Free_Slapaf
End Module Slapaf_Info
