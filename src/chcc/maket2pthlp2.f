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
        subroutine makeT2ptHlp2 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,     &
     &                          dimi,dimij,dimapp,dimbpp,dimap,dimabp)
!
!       this routine do:
!       define T2(+-)((ab)",ij) = Tau((ab)',i,j) +- Tau((ab)',j,i)
!       for the case: aGrp=bGrp but aSGrp>bSGrp
!
!       parameter description:
!       T2p     - array for T2+ (O)
!       Tau     - array for Tau (I)
!       xGrp    - Groups of a',b' (I)
!       xSGrp   - SubGroups of a",b" (I)
!       key     - 0 - T2+; 1 = T2- will be produced
!       dimx    - Dimension of i,(i>=j),a",b",a',(a>=b)' (I)
!
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
!
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimbpp,dimap,dimabp
        real*8 T2p(1:dimapp,1:dimbpp,1:dimij)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
!
!
!       help variables
        integer i,j,ij,app,bpp,ap,abp,appAdd,bppAdd
!
!
!1      def appAdd,bppAdd
!
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
!
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
!
!
!2      define T2(+,-)
!
        if (key.eq.0) then
!
!2.1    define T2+
        ij=0
        do i=1,dimi
        do j=1,i
        ij=ij+1
          do app=1,dimapp
          ap=appAdd+app
          abp=ap*(ap-1)/2+bppAdd
          do bpp=1,dimbpp
          abp=abp+1
            T2p(app,bpp,ij)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
          end do
        end do
        end do
!
        else
!
!2.2    define T2-
        ij=0
        do i=2,dimi
        do j=1,i-1
        ij=ij+1
          do app=1,dimapp
          ap=appAdd+app
          abp=ap*(ap-1)/2+bppAdd
          do bpp=1,dimbpp
          abp=abp+1
            T2p(app,bpp,ij)=Tau(abp,i,j)-Tau(abp,j,i)
          end do
          end do
        end do
        end do
!
        end if
!
        call mv0sv (dimij*dimapp*dimbpp,dimij*dimapp*dimbpp,            &
     &              T2p(1,1,1),0.5d0)
!
!
        return
! Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
