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
Module Slapaf_Parameters
implicit none
Private
Public:: iRow, iRow_c, iInt, nFix, ddV_Schlegel
Integer:: iRow=0
Integer:: iRow_c=0
Integer:: iInt=0
Integer:: nFix=0

Logical:: ddV_Schlegel=.False.
End Module Slapaf_Parameters
