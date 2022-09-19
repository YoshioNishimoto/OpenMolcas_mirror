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
        subroutine generate_juzekOE (oe,oeh,oep,no,nv)
!
! this routine do :
!
! modifies standar oe record to one used by DIRCC
!
! oeh = ( oe_occ(alpha) ... oe_occ(beta) )   (alpha=beta for this
!                                             implementation)
!
! oep = ( oe_virt(alpha) ... oe_virt(beta) ) (alpha=beta for this
!                                             implementation)
        implicit none
        integer i,no,nv
!
        real*8 oe(*),oeh(*),oep(*)
!
!
        do i=1,no
        oeh(i)=oe(i)
        oeh(i+no)=oe(i)
        end do
!
        do i=1,nv
        oep(i)=oe(no+i)
        oep(i+nv)=oe(no+i)
        end do
!
        return
        end
