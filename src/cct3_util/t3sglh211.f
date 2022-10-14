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
       subroutine t3sglh211 (w,dima,dimab,dimc,s1,d1,ns)
!
!     this routine add following contribution to W
!     for syma=symb;symc
!
!     W(ab,c) <-  + S1 _i(a) . D1 _jk(b,c)
!     - S1 _i(b) . D1 _jk(a,c)
!
!     w      - W  matrix (I/O)
!     dima   - dimension of a (b) index (I)
!     dimab  - dimension of ab  index (I)
!     dimc   - dimension of c index (I)
!     s1     - S1 matrix (I)
!     d1     - D1 matrix (I)
!     ns     - signum of the contribution (+-1) (I)
!
       integer dima,dimab,dimc,ns
       real*8 w(1:dimab,1:dimc)
       real*8 s1(1:dima)
       real*8 d1(1:dima,1:dimc)
!
!     help variables
!
       integer a,b,c,ab
       real*8 s
!
       if (ns.eq.1) then
!     phase + 1
!
       do 100 c=1,dimc
       ab=0
       do 101 a=2,dima
       s=s1(a)
       do 102 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)+d1(b,c)*s
 102    continue
 101    continue
 100    continue
!
       do 110 c=1,dimc
       ab=0
       do 111 a=2,dima
       s=d1(a,c)
       do 112 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)-s1(b)*s
 112    continue
 111    continue
 110    continue
!
       else
!     phase - 1
!
       do 200 c=1,dimc
       ab=0
       do 201 a=2,dima
       s=s1(a)
       do 202 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)-d1(b,c)*s
 202    continue
 201    continue
 200    continue
!
       do 210 c=1,dimc
       ab=0
       do 211 a=2,dima
       s=d1(a,c)
       do 212 b=1,a-1
       ab=ab+1
       w(ab,c)=w(ab,c)+s1(b)*s
 212    continue
 211    continue
 210    continue
!
       end if
!
       return
       end
