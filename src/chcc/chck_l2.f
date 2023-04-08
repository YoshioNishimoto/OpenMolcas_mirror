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
        subroutine Chck_L2 (L2,dima,dimb,adda,addb)
!
!        this routine test L2
!
        implicit none
#include "chcc1.fh"
        integer dima,dimb,adda,addb
         real*8 L2(1:nc,1:dima,1:dimb)
!
!        help var
        integer m,a,b,bad,i,ntot
        real*8 s
!
        bad=0
        ntot=0
!
        do b=addb+1,addb+dimb
        do a=adda+1,adda+dima
        do m=1,nc
!
          s=L2k(m,a,b)
           do i=1,no
           s=s-L1k(m,i,a)*T1c(b,i)
           end do
!
           if (abs(L2(m,a-adda,b-addb)-s).gt.1.0d-10) then
            bad=bad+1
             L2(m,a-adda,b-addb)=s
          end if
          ntot=ntot+1

        end do
        end do
        end do
!
        write (6,*) ' L2   ',bad,ntot
!
        return
        end
