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
        subroutine Ext_W4hlp1 (V2,M1,                                   &
     &                         nc,dima,dimab,dimapp,dimabpp,addapp)
!
!        this routine do:
!        Extract M1(m,a"b") <- V2(m,a'b')
!
        implicit none
        integer nc,dima,dimab,dimapp,dimabpp,addapp
        real*8 V2(1:nc,1:dimab)
        real*8 M1(1:nc,1:dimabpp)
!
!        help variables
!
        integer a,ab,app,bpp,abpp,m
!
        abpp=0
        do app=1,dimapp
        a=addapp+app
        ab=a*(a-1)/2+addapp
        do bpp=1,app
        ab=ab+1
        abpp=abpp+1
!
          do m=1,nc
          M1(m,abpp)=V2(m,ab)
          end do
!
        end do
        end do
!
        return
! Avoid unused argument warnings
        if (.false.) Call Unused_integer(dima)
        end
