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
        subroutine pack32_23 (AA,BB,d1,d2)
!
! this routine do :
!
! A(a,b,c) -> B(a,bc)  bc : b>=c,  dimb must eq dimc
!
        implicit none
        integer d1,d2,a,b,c,bc
        real*8 AA(d1,d2,d2),BB(d1,(d2*(d2+1)/2))
!
        bc=0
        do b=1,d2
        do c=1,b
        do a=1,d1
!
        bc=bc+1
        BB(a,bc)=AA(a,b,c)
        end do
        end do
        end do
!
        return
        end
