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
        subroutine Calc_Bc
!
!        this routine calc Bc
!        Bc(a,b,be,ga) =        (a,be|b,ga)
!                     - S(i)   (a,be,|b,i).t1(ga,i)
!                     - S(i)   (a,i,|b,ga).t1(be,i)
!                     + S(i,j) (a,i|b,j).t1(be,i).t1(ga,j)
!
        implicit none
#include "chcc1.fh"
!
!        help var
        integer i,a,b,be,ga
!       integer j
        real*8 s
!
        do ga=1,nv
        do be=1,nv
        do b=1,nv
        do a=1,nv
!
!
!1            (a,be|b,ga)
          s=Q4(be,a,ga,b)
!
!2,3          - S(i)   (a,be,|b,i).t1(ga,i)
!         - S(i)   (a,i,|b,ga).t1(be,i)
          do i=1,no
           s=s-Q3(a,be,b,i)*T1c(ga,i)
           s=s-Q3(b,ga,a,i)*T1c(be,i)
          end do

!4        + S(i,j) (a,i|b,j).t1(be,i).t1(ga,j)
!          do j=1,no
!          do i=1,no
!          s=s+Q21(a,i,b,j)*T1c(be,i)*T1c(ga,j)
!          end do
!          end do
!
!
          Bc(a,b,be,ga)=s
!
        end do
        end do
        end do
        end do
!
!
        return
        end
