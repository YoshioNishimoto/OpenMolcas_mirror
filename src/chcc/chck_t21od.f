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
        subroutine Chck_T21od (T21,beSGrp,gaSGrp)
!
!        test T2n+
!
        implicit none
#include "chcc1.fh"
        integer beSGrp,gaSGrp
        real*8 T21(1:32,1:32,1:no*(no-1)/2)
        integer a,b,u,v,be,ga,bega,uv,bad,gap,bep
        real*8 s
!
        if (beSGrp.eq.2) then
          bep=nv/2
        else
          bep=0
        end if
!
        if (gaSGrp.eq.2) then
          gap=nv/2
        else
          gap=0
        end if

        bad=0
!
        uv=0
        do u=2,no
        do v=1,u-1
        uv=uv+1
!
          bega=0
          do be=1,nv/2
          do ga=1,nv/2
          bega=bega+1
!
            s=0.0d0
            do a=1,nv
            b=a
            s=s+(Q4(b,gap+ga,a,bep+be)+Q4(b,bep+be,a,gap+ga))*          &
     &          (T2c(b,a,v,u)+T2c(b,a,u,v))/4
            end do
!
            s=0.0d0
            do a=2,nv
            do b=1,a-1
            s=s+(Q4(b,gap+ga,a,bep+be)-Q4(b,bep+be,a,gap+ga))*          &
     &          (T2c(b,a,v,u)-T2c(b,a,u,v))/2
            end do
            end do
!
          if (abs(T21(be,ga,uv)-s).gt.1.0d-10) then
            bad=bad+1
!        write (6,99) be,ga,u,v
!99        format (4(i3,1x))
          end if
          T21(be,ga,uv)=s
!
          end do
          end do
!
        end do
        end do
!
        if (bad.eq.0) then
        write (6,*) ' Chck T2 OK ', bad
        else
        write (6,*) ' Chck T2 Bug !!!!!!! ', bad
        end if
!
        return
        end
