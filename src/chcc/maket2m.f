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
        subroutine MakeT2m (T2m,Tau,aGrp,bGrp,aSGrp,bSGrp,keyT)
!
!       this routine do:
!       Make T2p(i>,(a>b)")  = Tau(i,j,(a>=b)")-Tau(j,i,a>=b)")
!                     from     Tau((a>=b)',i,j)
!        or Transposed (T(ab",ij)
!
!       parameter description:
!       T2m    - T2- array (O)
!       Tau    - Tau array (I)
!       xGrp   - Group of a,b (I)
!       xSGrp  - SubGroup of a,b (I)
!        keyT   - 0 - make T(ij,ab")
!                1 - make T(ab",ij)
!
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
!
        real*8 T2m(1)
        real*8 Tau(1)
        integer aGrp,bGrp,aSGrp,bSGrp,keyT
!
!       help variables
        integer dimi,dimij,dimap,dimbp,dimapp,dimbpp,dimabp,dimabpp
!
        dimi=no
        dimij=no*(no-1)/2
!
        dimap=DimGrpa(aGrp)
        dimbp=DimGrpa(bGrp)
        if (aGrp.eq.bGrp) then
          dimabp=dimap*(dimap+1)/2
        else
          dimabp=dimap*dimbp
        end if
!
        dimapp=DimSGrpa(aSGrp)
        dimbpp=DimSGrpa(bSGrp)
        if (aSGrp.eq.bSGrp) then
          dimabpp=dimapp*(dimapp-1)/2
        else
          dimabpp=dimapp*dimbpp
        end if
!
        if (keyT.eq.0) then
!        T-(ij,ab") case
!
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2pHlp1 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,    &
     &                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2pHlp2 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,    &
     &                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2pHlp3 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,      &
     &                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
!
        else
!        T-(ab",ij) case
!
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2ptHlp1 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,   &
     &                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2ptHlp2 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,   &
     &                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2ptHlp3 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,     &
     &                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
!
        end if
!
        return
        end
