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
       subroutine t3sglh132 (w,dima,dimb,dimbc,s2,d2,ns)
!
!     this routine add following contribution to W
!     for syma=symb > symc
!
!     W(a,bc) <-  - S2 _i(b) . D2 _jk(a,c)
!     + S2 _i(c) . D2 _jk(a,b)
!
!     w      - W  matrix (I/O)
!     dima   - dimension of a (b) index (I)
!     dimab   - dimension of ab index (I)
!     dimc   - dimension of c index (I)
!     s2     - S2 matrix (I)
!     d2     - D2 matrix (I)
!     ns     - signum of the contribution (+-1) (I)
!
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s2(1:dimb)
       real*8 d2(1:dima,1:dimb)
!
!     help variables
!
       integer a,b,c,bc
       real*8 s
!
       if (ns.eq.1) then
!     phase +1
!
       bc=0
       do 110 b=2,dimb
       s=s2(b)
       do 111 c=1,b-1
       bc=bc+1
       do 112 a=1,dima
       w(a,bc)=w(a,bc)-d2(a,c)*s
 112    continue
 111    continue
 110    continue
!
       bc=0
       do 120 b=2,dimb
       do 121 c=1,b-1
       s=s2(c)
       bc=bc+1
       do 122 a=1,dima
       w(a,bc)=w(a,bc)+d2(a,b)*s
 122    continue
 121    continue
 120    continue
!
       else
!     phase -1
!
       bc=0
       do 210 b=2,dimb
       s=s2(b)
       do 211 c=1,b-1
       bc=bc+1
       do 212 a=1,dima
       w(a,bc)=w(a,bc)+d2(a,c)*s
 212    continue
 211    continue
 210    continue
!
       bc=0
       do 220 b=2,dimb
       do 221 c=1,b-1
       s=s2(c)
       bc=bc+1
       do 222 a=1,dima
       w(a,bc)=w(a,bc)-d2(a,b)*s
 222    continue
 221    continue
 220    continue
!
       end if
!
       return
       end
