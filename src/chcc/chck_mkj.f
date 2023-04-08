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
        subroutine Chck_mkJ
!
!        make J(be,u,i,a)
!
        implicit none
#include "chcc1.fh"
!
        integer be,u,i,a
        integer j,b
        real*8 sj
!
        do a=1,nv
        do i=1,no
        do u=1,no
        do be=1,nv
!
          sj=0.0d0
!
          do j=1,no
           sj=sj-Q1(a,i,j,u)*T1c(be,j)
          end do
!
          do b=1,nv
           sj=sj+Q3(b,be,a,i)*T1c(b,u)
          end do
!
          do j=1,no
          do b=1,nv
           sj=sj+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(be,b,u,j)/2
           sj=sj-Q21(a,i,b,j)*(T2c(b,be,u,j)/2+T1c(b,u)*T1c(be,j))
          end do
          end do
!
          sj=sj+Q21(be,u,a,i)
!
          Jc(be,i,u,a)=sj
!
        end do
        end do
        end do
        end do
!
        write (6,*) ' J done'
!
        return
        end
