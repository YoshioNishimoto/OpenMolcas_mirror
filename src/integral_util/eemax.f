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
      Subroutine EEMax(n,A,B)
      Implicit None
      Integer n
      Real*8 A(n), B(n)

      Integer i

      Do i = 1, n
         If (B(i).gt.A(i)) A(i)=B(i)
      End Do

      End Subroutine EEMax
