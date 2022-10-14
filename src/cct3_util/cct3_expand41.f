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
       subroutine cct3_expand41 (a,b,dimpq,dimr,dims,dimp)
!
!     expand a(pq,r,s) -> - b(p,q,s,r)
!     assumption: p>q + antisymetry
!     RISC version
!
       integer dimpq,dimr,dims,dimp
       real*8 a(1:dimpq,1:dimr,1:dims)
       real*8 b(1:dimp,1:dimp,1:dims,1:dimr)
!
!     help variables
!
       integer p,q,r,s,pq
       real*8 scalar
!
       if (dimp.gt.1) then
!
       do 100 s=1,dims
       do 101 r=1,dimr
!
       pq=0
       do 102 p=2,dimp
       do 103 q=1,p-1
       pq=pq+1
!
       scalar=a(pq,r,s)
       b(p,q,s,r)=-scalar
       b(q,p,s,r)=scalar
!
 103    continue
 102    continue
 101    continue
 100    continue
!
       end if
!
       do r=1,dimr
       do s=1,dims
       do p=1,dimp
       b(p,p,s,r)=0.0d0
       end do
       end do
       end do
!
       return
       end
