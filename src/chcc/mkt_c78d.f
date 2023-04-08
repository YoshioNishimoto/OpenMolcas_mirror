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
        subroutine MkT_C78d (T2,Tp,Tm,dimbe,dimbepp,addbepp,no)
!
!        this routine do:
!       T2n(be',ga',u,v) <<-
!        C7                    T2+(bega",uv)
!        C8                    T2-(bega",uv)
!        for beGrp=gaGrp, beSGrp=gaSGrp
!        N.B. calc only contributions to be',ga' (not ga',be')
!
        implicit none
        integer dimbe,dimbepp,addbepp,no
        real*8 T2(1:dimbe,1:dimbe,1:no,1:no)
        real*8 Tp(1:dimbepp*(dimbepp+1)/2,1:no*(no+1)/2)
        real*8 Tm(1:dimbepp*(dimbepp-1)/2,1:no*(no-1)/2)
!
!        help variables
        integer u,v,be,ga,uv,bega,bep,gap
        real*8 fact
!
!
!1        Distribute symmetric T2+ on proper possitions
!
        uv=0
        do u=1,no
        do v=1,u
        uv=uv+1
        if (u.eq.v) then
          fact=0.5d0
        else
          fact=1.0d0
        end if
!
!          case be".ne.ga"
          bep=addbepp+1
          do be=2,dimbepp
          bega=be*(be-1)/2
          bep=bep+1
          gap=addbepp
          do ga=1,be-1
          bega=bega+1
          gap=gap+1
!
            T2(bep,gap,u,v)=T2(bep,gap,u,v)+Tp(bega,uv)*fact
            T2(bep,gap,v,u)=T2(bep,gap,v,u)+Tp(bega,uv)*fact
!             T2(gap,bep,u,v)=T2(gap,bep,u,v)+Tp(bega,uv)*fact
!            T2(gap,bep,v,u)=T2(gap,bep,v,u)+Tp(bega,uv)*fact
!
          end do
          end do
!
!          case be=ga
          bep=addbepp
          do be=1,dimbepp
          bega=be*(be+1)/2
          bep=bep+1
!
            T2(bep,bep,u,v)=T2(bep,bep,u,v)+Tp(bega,uv)*fact
            T2(bep,bep,v,u)=T2(bep,bep,v,u)+Tp(bega,uv)*fact
!
          end do
!
        end do
        end do
!
!
!2        Distribute anti-symmetric T2- on proper possitions
!
        uv=0
        do u=2,no
        do v=1,u-1
        uv=uv+1
!
          bep=addbepp+1
          bega=0
          do be=2,dimbepp
          bep=bep+1
          gap=addbepp
          do ga=1,be-1
          bega=bega+1
          gap=gap+1
!
            T2(bep,gap,u,v)=T2(bep,gap,u,v)+Tm(bega,uv)
            T2(bep,gap,v,u)=T2(bep,gap,v,u)-Tm(bega,uv)
!            T2(gap,bep,u,v)=T2(gap,bep,u,v)-Tm(bega,uv)
!            T2(gap,bep,v,u)=T2(gap,bep,v,u)+Tm(bega,uv)
!
          end do
          end do
!
        end do
        end do
!
!
        return
        end
