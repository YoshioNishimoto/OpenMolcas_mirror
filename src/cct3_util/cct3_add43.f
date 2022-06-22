************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine cct3_add43 (a,b,q,dimp,dimqr,dimr,fact)

c     this routine do:
c     B(p,qr) <-- fact * A(p,r) for given q
c
#include "t31.fh"
       integer dimp,dimqr,dimr,q
       real*8 fact
       real*8 b(1:dimp,1:dimqr)
       real*8 a(1:dimp,1:dimr)
c
c     help variable
c
       integer p,qr,rq,r
c
       if (q.eq.1) goto 101
c
       qr=nshf(q)
       do 100 r=1,q-1
       qr=qr+1
c
       do 50 p=1,dimp
       b(p,qr)=b(p,qr)+fact*a(p,r)
 50     continue
c
 100    continue
c
 101    if (q.eq.dimr) then
       return
       end if
c
c
       do 200 r=q+1,dimr
       rq=nshf(r)+q
       do 150 p=1,dimp
       b(p,rq)=b(p,rq)-fact*a(p,r)
 150    continue
c
 200    continue
c
       return
       end
