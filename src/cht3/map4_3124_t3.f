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
        subroutine Map4_3124_t3 (A,B,d1,d2,d3,d4)
!
!       this routine do:
!       map B(2314) <- A(1234)
!
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d2,1:d3,1:d1,1:d4)
!
!       help variables
        integer i1,i2,i3,i4
!
        do i4=1,d4
        do i1=1,d1
        do i3=1,d3
        do i2=1,d2
        b(i2,i3,i1,i4)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
!
        return
        end
