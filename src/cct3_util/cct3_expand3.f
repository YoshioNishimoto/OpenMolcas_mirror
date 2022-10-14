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
       subroutine cct3_expand3 (a,b,dimp,dimqr,dimq)
!
!     expand a(p,qr) -> b(p,q,r)
!     assumption: q>r, a(p,q,r)=-a(p,r,q)
!     RISC version
!
       integer dimp,dimq,dimqr
       real*8 a(1:dimp,1:dimqr)
       real*8 b(1:dimp,1:dimq,1:dimq)
!
!     help variables
!
       integer p,q,r,qr
       real*8 scalar
!
       if (dimq.gt.1) then
!
       qr=0
       do 100 q=2,dimq
       do 101 r=1,q-1
       qr=qr+1
!
       do p=1,dimp
       scalar=a(p,qr)
       b(p,q,r)=scalar
       b(p,r,q)=-scalar
       end do
 101    continue
 100    continue
!
       end if
!
       do q=1,dimq
       do p=1,dimp
       b(p,q,q)=0.0d0
       end do
       end do
!
       return
       end
