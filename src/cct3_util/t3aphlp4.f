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
       subroutine t3aphlp4 (a,b,dimp,dimpq,dimpqr,ns,szkey)
!
!     this routine realize following procedure
!     for symp=symq=symr
!
!     B(pqr)=B(pqr)+ns.(A1(qr,p)-A2(pr,q)+A3(pq,r))
!
!     a      - A  matrix (I)
!     b      - B  matrix (I/O)
!     dimp   - dimension of p (q,r) index (I)
!     dimpq  - dimension of pq index (I)
!     dimpqr - dimension of pqr (I)
!     ns     - singum of the permutation (+-1) (I)
!     szkey  - set zero key (I)
!     = 0 no vanishing
!     = 1 set B=0 at the beggining
!
!     N.B. tie sumacie by sa mozno dali trochu vylepsit
!
       integer dimp,dimpq,dimpqr,ns,szkey
       real*8 a(1:dimpq,1:dimp)
       real*8 b(1:dimpqr)
!
!     help variables
!
       integer p,q,r,pqr,pq,pr,qr,pq0
!
!
       if (szkey.eq.1) then
       call cct3_mv0zero (dimpqr,dimpqr,b)
       end if
!
       if (ns.eq.1) then
!     phase +1
!
       pqr=0
!
       do 100 p=3,dimp
       pq0=(p-1)*(p-2)/2
       qr=0
       do 101 q=2,p-1
       pq=pq0+q
       pr=(p-1)*(p-2)/2
       do 102 r=1,q-1
       pr=pr+1
       qr=qr+1
       pqr=pqr+1
       b(pqr)=b(pqr)+a(qr,p)-a(pr,q)+a(pq,r)
 102    continue
 101    continue
 100    continue
!
       else
!     phase -1
!
       pqr=0
!
       do 200 p=3,dimp
       pq0=(p-1)*(p-2)/2
       qr=0
       do 201 q=2,p-1
       pq=pq0+q
       pr=(p-1)*(p-2)/2
       do 202 r=1,q-1
       pr=pr+1
       qr=qr+1
       pqr=pqr+1
       b(pqr)=b(pqr)-a(qr,p)+a(pr,q)-a(pq,r)
 202    continue
 201    continue
 200    continue
!
!
       end if
!
       return
       end
