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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************
      SubRoutine PickUp(Tri,Vec,n)
!***********************************************************************
!                                                                      *
!     purpose: Pick up diagonal elements from triangular matrix and    *
!              put it to a vector                                      *
!                                                                      *
!     input:                                                           *
!       Tri     : triangular matrix                                    *
!       n       : dimension                                            *
!                                                                      *
!     output:                                                          *
!       Vec     : vector                                               *
!                                                                      *
!***********************************************************************
!
      Implicit None
!
      Integer n
      Real*8 Tri(n*(n+1)/2),Vec(n)

      Integer ij, i
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!
      ij = 0
      Do i = 1, n
         ij = ij + i
         Vec(i) = Tri(ij)
      End Do
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      End SubRoutine PickUp
