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
        subroutine MkOE (OE)
!
!        OEo(i) <- OE(i)
!        OEv(a) <- OE(a)
!
        implicit none
#include "chcc1.fh"
        real*8 OE(1:(nv+no))
!
!        help variables
        integer a,i
!
        do i=1,no
          OEo(i)=OE(i)
        end do
!
        do a=1,nv
          OEv(a)=OE(no+a)
        end do
!
        return
        end
