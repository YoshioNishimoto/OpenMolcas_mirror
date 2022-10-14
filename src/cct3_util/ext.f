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
       subroutine ext (wrk,wrksize,                                     &
     & nind,exttyp,u,v,x,ssu,ssv,ssx,mapda,mapia,ssa,                   &
     & possb0,mapdb,mapib,ssb,rc)
!
!     thsi routine realize extraction
!
!     A(indA) -> B_u(indB) for given u
!
!     nind   - # of indexex in matrix A  (Input)
!     exttyp - type of extraction :  (Input)
!     1 - A(pqrs) -> B_p (qrs); A(pqr) -> B_p(qr) ; A(pq) -> B_p(q)
!     2 - A(pqrs) -> B_q (prs); A(pqr) -> B_q(p,r); A(pq) -> B_q(p)
!     3 - A(pqrs) -> B_r (pqs); A(pqr) -> B_r(pq)
!     4 - A(pqrs) -> B_s (pqr)
!     5 - A(pqrs) -< A_pq(rs) ; A(pqr) -> B_pq(r)
!     6 - A(pqrs) -> B_qr(p,s); A(pqr) -> B_pq(r)
!     7 - A(pqrs) -> B_rs(pq)
!     8 - A(pqrs) -> B_pqr(s)
!     9 - A(pqrs) -> B_qrs(p)
!
!     u      - value of 1.st fix index (I)
!     v      - value of 2.nd fix index (I)
!     x      - value of 3.rd fix index (I)
!     ssu    - symmetry of 1.st fix index (I)
!     ssv    - symmetry of 2.nd fix index (I)
!     ssx    - symmetry of 3.rd fix index (I)
!     mapda  - direct map matrix corresponding to A  (Input)
!     mapia  - inverse map matrix corresponding to A  (Input)
!     ssa    - overall symetry state  of matrix A  (Input)
!     possb0 - initial possition of matrix B in WRK  (Input)
!     mapdb  - direct map matrix corresponding to B  (Output)
!     mapib  - inverse map matrix corresponding to B  (Output)
!     ssb    - overall symetry state  of matrix B  (Output)
!     rc     - return (error) code  (Output)
!
!     Table of extractions
!
!     nind  exttyp          Operation             Implementation
!     4       1     A(p,q,r,s) -> B _p(q,r,s)          Yes
!     A(pq,r,s)  -> B _p(q,r,s)          Yes
!     A(p,qr,s)  -> B _p(qr,s)           NCI
!     A(p,q,rs)  -> B _p(q,rs)           Yes
!     B(pq,rs)   -> B _p(q,rs)           Yes
!
!     4       2     A(p,q,r,s) -> B _q(p,r,s)          Yes
!     A(pq,r,s)  -> B _q(p,r,s)     No (use 4,1)
!     A(p,qr,s)  -> B _q(p,r,s)          NCI
!     A(p,q,rs)  -> B _q(p,rs)           Yes
!     B(pq,rs)   -> B _q(p,rs)      No (use 4,1)
!
!     4       3     A(p,q,r,s) -> B _r(p,q,s)          Yes
!     A(pq,r,s)  -> B _r(pq,s)           Yes
!     A(p,qr,s)  -> B _r(p,q,s)     No (use 4.2 - NCI)
!     A(p,q,rs)  -> B _r(p,q,s)          Yes
!     B(pq,rs)   -> B _r(pq,s)           Yes
!
!     4       4     A(p,q,r,s) -> B _s(p,q,r)          Yes
!     A(pq,r,s)  -> B _s(pq,r)           Yes
!     A(p,qr,s)  -> B _s(p,qr)           NCI
!     A(p,q,rs)  -> B _s(p,q,r)     No (use 4,3)
!     B(pq,rs)   -> B _s(pq,r)      No (use 4,3)
!
!     4       5     A(p,q,r,s) -> B _p_q (r,s)         Yes
!     A(pq,r,s)  -> B _pq (r,s)          NCI
!     A(p,qr,s)  -> B _p_q (r,s)         NCI
!     A(p,q,rs)  -> B _p_q (rs)          NCI
!     A(pq,rs)   -> B _pq (rs)           Yes
!
!     4       6     any case                           NCI
!
!     4       7     A(p,q,r,s) -> B _r_s (p,q)         Yes
!     A(pq,r,s)  -> B _r_s (pq)          NCI
!     A(p,qr,s)  -> B _r_s (p,q)         NCI
!     A(p,q,rs)  -> B _r_s (rs)          Yes
!     A(pq,rs)   -> B _pq (rs)           Yes
!
!     4      8,9    any case                           NCI
!
!     3   any case  any case                           NCI
!
!     2       1     A(p,q)     -> B _p (q)             Yes
!     A(pq)      -> B _p (q)             NCI
!
!     2       2     A(p,q)     -> B _q (p)             Yes
!     A(pq)      -> B _q (p)             NCI
!
!
#include "t31.fh"
#include "wrk.fh"
!
       integer nind,exttyp,ssa,u,v,x,ssu,ssv,ssx,possb0,ssb,rc
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
!
!     help variables
       integer typa,typb,possa,possb,symp,symq,symr,syms
       integer ia,ib,key,posst
       integer nhelp1,nhelp2,jjind,signum
       integer dimp,dimq,dimr,dims
!
!     To fix some warnings
      symr=0
!0.*  some general tests
!
       if (mapda(0,6).eq.2) then
!     RC=2  : nind=4, typA=2 (NCI)
       rc=2
       return
       end if

!
!0.*  def typa and ssB
!
       typa=mapda(0,6)
       if (exttyp.le.4) then
!     one extetnal index
       ssb=mmul(ssa,ssu)
       else if (exttyp.le.7) then
!     two extetnal indexes
       nhelp1=mmul(ssu,ssv)
       ssb=mmul(ssa,nhelp1)
       else
!     3 external indeses
       nhelp1=mmul(ssu,ssv)
       nhelp2=mmul(nhelp1,ssx)
       ssb=mmul(nhelp2,ssa)
       end if
!
!
!
       if (nind.eq.4) then
!
!4    A(pqrs) -> B(stv)
!
       if (exttyp.eq.1) then
!
!4.1  A(pqrs) -> B_p(qrs)
!
!
       if (typa.eq.0) then
!
!4.1.0**** case A(p,q,r,s) -> B_p(q,r,s) ****
!
!4.1.0.*    get mapdb,mapib
       typb=0
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 410 ib=1,mapdb(0,5)
!
!4.1.0.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.1.0.*      find porper A
       ia=mapia(ssu,symq,symr)
!
!4.1.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.1.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimq*dimr*dims
!
!4.1.0.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
!
 410    continue
!
!
       else if (typa.eq.1) then
!
!4.1.1**** case A(pq,r,s) -> B_p(q,r,s) ****
!
       typb=0
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 411 ib=1,mapdb(0,5)
!
!4.1.1.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.1.1.*      def key su>sq - 1 ; su=sq - 2 ; su<sq - 3
       if (ssu.gt.symq) then
       key=1
       else if (ssu.eq.symq) then
       key=2
       else
       key=3
       end if
!
!4.1.1.*      find porper A
       if (key.lt.3) then
       ia=mapia(ssu,symq,symr)
       else
       ia=mapia(symq,ssu,symr)
       end if
!
!4.1.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.1.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
!4.1.1.*      realize extraction
       if (key.eq.1) then
       nhelp1=dimq*dimr*dims
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
       else if (key.eq.2) then
       nhelp1=dimp*(dimp-1)/2
       nhelp2=dimr*dims
       call exth4 (wrk(possa),wrk(possb),dimp,nhelp1,nhelp2,u)
       else
!     key=3
       nhelp1=dimr*dims
       call exth3 (wrk(possa),wrk(possb),dimq,dimp,nhelp1,u,-1)
       end if
!
 411    continue
!
!
       else if (typa.eq.2) then
!
!4.1.2**** case A(p,qr,s) -> B_p(qr,s) ****
!     NCI
!
!
       else if (typa.eq.3) then
!
!4.1.3**** case A(p,q,rs) -> B_p(q,rs) ****
!
       typb=2
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 413 ib=1,mapdb(0,5)
!
!4.1.3.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.1.3.*      find porper A
       ia=mapia(ssu,symq,symr)
!
!4.1.3.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.1.3.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       if (symr.eq.syms) then
       nhelp1=dimq*dimr*(dimr-1)/2
       else
       nhelp1=dimq*dimr*dims
       end if
!
!4.1.3.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
!
 413    continue
!
!
       else if (typa.eq.4) then
!
!4.1.4**** case A(pq,rs) -> B_p(q,rs) ****
!
       typb=2
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 414 ib=1,mapdb(0,5)
!
!4.1.4.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.1.4.*      def key su>sq - 1 ; su=sq - 2 ; su<sq - 3
       if (ssu.gt.symq) then
       key=1
       else if (ssu.eq.symq) then
       key=2
       else
       key=3
       end if
!
!4.1.4.*      find porper A
       if (key.lt.3) then
       ia=mapia(ssu,symq,symr)
       else
       ia=mapia(symq,ssu,symr)
       end if
!
!4.1.4.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.1.4.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
!4.1.4.*      realize extraction
!
       if (key.eq.1) then
       if (symr.eq.syms) then
       nhelp1=dimq*dimr*(dimr-1)/2
       else
       nhelp1=dimq*dimr*dims
       end if
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
!
       else if (key.eq.2) then
       nhelp1=dimp*(dimp-1)/2
       if (symr.eq.syms) then
       nhelp2=dimr*(dimr-1)/2
       else
       nhelp2=dimr*dims
       end if
       call exth4 (wrk(possa),wrk(possb),dimp,nhelp1,nhelp2,u)
!
       else
!     key=3
       if (symr.eq.syms) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
       call exth3 (wrk(possa),wrk(possb),dimq,dimp,nhelp1,u,-1)
       end if
!
 414    continue
!
!
       end if
!
!
!
       else if (exttyp.eq.2) then
!
!4.2  A(pqrs) -> B_q(prs)
!
!
       if ((typa.eq.0).or.(typa.eq.3)) then
!
!4.2.0case A(p,q,r,s) -> B_q(p,r,s)
!
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,3),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 420 ib=1,mapdb(0,5)
!
!4.2.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.2.0.*      find porper A
       ia=mapia(symp,ssu,symr)
!
!4.2.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.2.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimr*dims
!
!4.2.0.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),dimp,dimq,nhelp1,u,1)
!
 420    continue
!
!
       else if (typa.eq.1) then
!
!4.2.1case A(pq,r,s) -> B_q(p,r,s)
!
!     RC=3 : nind=4, typA=1, exttyp=2 (NI - Use exttyp 1)
       rc=3
       return
!
!
       else if (typa.eq.2) then
!
!4.2.2case A(p,qr,s) -> B_q(p,r,s)
!
!     RC=4 : nind=4, typA=2 (NCI)
       rc=4
       return
!
!
       else if (typa.eq.3) then
!
!4.2.3case A(p,q,rs) -> B_q(p,rs)
!
       typb=2
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,3),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 423 ib=1,mapdb(0,5)
!
!4.2.3.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.2.3.*      find porper A
       ia=mapia(symp,ssu,symr)
!
!4.2.3.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.2.3.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       if (symr.eq.syms) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
!
!4.2.3.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),dimp,dimq,nhelp1,u,1)
!
 423    continue
!
!
       else if (typa.eq.4) then
!
!4.2.4case A(pq,r,s) -> B_q(p,r,s)
!
!     RC=5 : nind=4, typA=4, exttyp=2 (NI - Use exttyp 1)
       rc=5
       return
!
!
       end if
!
!
!
       else if (exttyp.eq.3) then
!
!4.3  A(pqrs) -> B_r(pqs)
!
!
       if (typa.eq.0) then
!
!4.3.0case A(p,q,r,s) -> B_r(p,q,s)
!
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 430 ib=1,mapdb(0,5)
!
!4.3.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.3.0.*      find porper A
       ia=mapia(symp,symq,ssu)
!
!4.3.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.3.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq
!
!4.3.0.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
!
 430    continue
!
!
       else if (typa.eq.1) then
!
!4.3.1case A(pq,r,s) -> B_r(pq,s)
!
       typb=1
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 431 ib=1,mapdb(0,5)
!
!4.3.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.3.1.*      find porper A
       ia=mapia(symp,symq,ssu)
!
!4.3.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.3.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       if (symp.eq.symq) then
       nhelp1=dimp*(dimp-1)/2
       else
       nhelp1=dimp*dimq
       end if
!
!4.3.1.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
!
 431    continue
!
!
       else if (typa.eq.2) then
!
!4.3.2case A(p,qr,s) -> B_r(p,q,s)
!
!     RC=6 : nind=4, typA=2 (NCI)
       rc=6
       return
!
!
       else if (typa.eq.3) then
!
!4.3.3case A(p,q,rs) -> B_r(p,q,s)
!
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 433 ib=1,mapdb(0,5)
!
!4.3.3.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.3.3.*      def key su>ss - 1 ; su=ss - 2 ; su<ss - 3
       if (ssu.gt.syms) then
       key=1
       else if (ssu.eq.syms) then
       key=2
       else
       key=3
       end if
!
!4.3.3.*      find porper A
       if (key.lt.3) then
       ia=mapia(symp,symq,ssu)
       else
       ia=mapia(symp,symq,symr)
       end if
!
!4.3.3.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.3.3.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
!4.3.3.*      realize extraction
       if (key.eq.1) then
       nhelp1=dimp*dimq
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
       else if (key.eq.2) then
       nhelp1=dimp*dimq
       nhelp2=dimr*(dimr-1)/2
       call exth5 (wrk(possa),wrk(possb),nhelp1,dimr,nhelp2,u)
       else
!     key=3
       nhelp1=dimp*dimq
       call exth3 (wrk(possa),wrk(possb),nhelp1,dims,dimr,u,-1)
       end if
!
 433    continue
!
!
       else if (typa.eq.4) then
!
!4.3.4case A(pq,rs) -> B_r(pq,s)
!
       typb=1
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 434 ib=1,mapdb(0,5)
!
!4.3.4.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
!
!4.3.4.*      def key su>ss - 1 ; su=ss - 2 ; su<ss - 3
       if (ssu.gt.syms) then
       key=1
       else if (ssu.eq.syms) then
       key=2
       else
       key=3
       end if
!
!4.3.4.*      find porper A
       if (key.lt.3) then
       ia=mapia(symp,symq,ssu)
       else
       ia=mapia(symp,symq,syms)
       end if
!
!4.3.4.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.3.4.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
!4.3.4.*      realize extraction
       if (key.eq.1) then
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)/2
       else
       nhelp1=dimp*dimq
       end if
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
       else if (key.eq.2) then
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)/2
       else
       nhelp1=dimp*dimq
       end if
       nhelp2=dimr*(dimr-1)/2
       call exth5 (wrk(possa),wrk(possb),nhelp1,dimr,nhelp2,u)
       else
!     key=3
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)/2
       else
       nhelp1=dimp*dimq
       end if
       nhelp1=nhelp1*dimr
       call exth2 (wrk(possa),wrk(possb),nhelp1,dims,u,-1)
       end if
!
 434    continue
!
!
       end if
!
!
!
       else if (exttyp.eq.4) then
!
!4.4  case A(pqrs) -> B_s(pqr)
!
!
       if (typa.eq.0) then
!
!4.4.0case A(p,q,r,s) -> B_s(p,q,r)
!
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,3),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 440 ib=1,mapdb(0,5)
!
!4.4.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       symr=mapdb(ib,5)
!
!4.4.0.*      find porper A
       ia=mapia(symp,symq,symr)
!
!4.4.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.4.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq*dimr
!
!4.4.0.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp1,dims,u,1)
!
 440    continue
!
!
       else if (typa.eq.1) then
!
!4.4.1case A(pq,r,s) -> B_s(pq,r)
!
       typb=1
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,3),0,ssb,   &
     & possb0,posst,mapdb,mapib)
!
       do 441 ib=1,mapdb(0,5)
!
!4.4.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       symr=mapdb(ib,5)
!
!4.4.1.*      find porper A
       ia=mapia(symp,symq,symr)
!
!4.4.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.4.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
!4.4.1.*      realize extraction
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)*dimr/2
       else
       nhelp1=dimp*dimq*dimr
       end if
       call exth2 (wrk(possa),wrk(possb),nhelp1,dims,u,1)
!
 441    continue
!
!
       else if (typa.eq.2) then
!
!4.4.2case A(p,qr,s) -> B_s(p,qr)
!
!     RC=7 : nind=4, typA=2  (NCI)
       rc=7
       return
!
!
       else if (typa.eq.3) then
!
!4.4.3case A(p,q,rs) -> B_s(p,q,r)
!
!     RC=8 : nind=4, typA=3, exttyp=4 (NI - Use exttyp 3)
       rc=8
       return
!
!
       else if (typa.eq.4) then
!
!4.4.4case A(pq,rs) -> B_s(pq,r)
!
!     RC=9 : nind=4, typA=4, exttyp=4 (NI - Use exttyp 3)
       rc=9
       return
!
!
       end if
!
       else if (exttyp.eq.5) then
!
!4.5  case A(pqrs) -> B_pq(rs)
!
       if (typa.eq.0) then
!
!4.5.0case A(p,q,r,s) -> B_p_q(r,s)
!
       typb=0
       call cct3_grc0 (2,typb,mapda(0,3),mapda(0,4),0,0,ssb,            &
     & possb0,posst,mapdb,mapib)
!
       do 450 ib=1,mapdb(0,5)
!
!4.5.0.*      def symmetry of all indices
       symr=mapdb(ib,3)
       syms=mapdb(ib,4)
!
!4.5.0.*      find porper A
       ia=mapia(ssu,ssv,symr)
!
!4.5.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.5.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq
       nhelp2=dimr*dims
       jjind=(v-1)*dimp+u
!
!4.5.0.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,1)
!
 450    continue
!
       else if (typa.eq.4) then
!
!4.5.4case A(pq,rs) -> B_pq(rs)
!
       typb=1
       call cct3_grc0 (2,typb,mapda(0,3),mapda(0,4),0,0,ssb,            &
     & possb0,posst,mapdb,mapib)
!
       if (ssu.ge.ssv) then
!4.5.4.1      case symp >= symq
!
       do 4541 ib=1,mapdb(0,5)
!
!4.5.4.1.*      def symmetry of all indices
       symr=mapdb(ib,3)
       syms=mapdb(ib,4)
!
!4.5.4.1.*      find porper A
       ia=mapia(ssu,ssv,symr)
!
!4.5.4.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.5.4.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
       if (ssu.eq.ssv) then
       nhelp1=dimp*(dimp-1)/2
       else
       nhelp1=dimp*dimq
       end if
!
       if (symr.eq.syms) then
       nhelp2=dimr*(dimr-1)/2
       else
       nhelp2=dimr*dims
       end if
!
!4.5.4.1.*        def joined index and signum
       signum=1
       if (ssu.eq.ssv)then
       if (u.gt.v) then
       jjind=nshf(u)+v
       else if (u.eq.v) then
       jjind=0
       signum=0
       else
       jjind=nshf(v)+u
       signum=-1
       end if
       else
       jjind=(v-1)*dimp+u
       end if
!
!4.5.4.1.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,           &
     & signum)
!
 4541   continue
!
       else
!4.5.4.2      case symp < symq
!
       do 4542 ib=1,mapdb(0,5)
!
!4.5.4.2.*      def symmetry of all indices
       symr=mapdb(ib,3)
       syms=mapdb(ib,4)
!
!4.5.4.2.*      find porper A
       ia=mapia(ssv,ssu,symr)
!
!4.5.4.2.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.5.4.2.*      def dimensions
       dimq=dimm(mapda(0,1),mapda(ia,3))
       dimp=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
       nhelp1=dimp*dimq
!
       if (symr.eq.syms) then
       nhelp2=dimr*(dimr-1)/2
       else
       nhelp2=dimr*dims
       end if
!
!4.5.4.2.*        def joined index and signum
       signum=-1
       jjind=(u-1)*dimp+v
!
!4.5.4.2.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,           &
     & signum)
!
 4542   continue
!
       end if
!
       else
!     RC=10 , nindA=4, exttyp=5, typa is nou 0,4 (NCI)
       rc=10
       return
       end if
!
       else if (exttyp.eq.7) then
!
!4.7  case         A(pqrs) -> B_rs(pq)
!
       if (typa.eq.0) then
!
!4.7.0case A(p,q,r,s) -> B_r_s(p,q)
!
       typb=0
       call cct3_grc0 (2,typb,mapda(0,1),mapda(0,2),0,0,ssb,            &
     & possb0,posst,mapdb,mapib)
!
       do 470 ib=1,mapdb(0,5)
!
!4.7.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
!
!4.7.0.*      find porper A
       ia=mapia(symp,symq,ssu)
!
!4.7.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.7.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq
       nhelp2=dimr*dims
       jjind=(v-1)*dimr+u
!
!4.7.0.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,1)
!
 470    continue
!
       else if (typa.eq.3) then
!
!4.7.3case A(p,q,rs) -> B_rs(p,q)
!
       typb=2
       call cct3_grc0 (2,typb,mapda(0,1),mapda(0,2),0,0,ssb,            &
     & possb0,posst,mapdb,mapib)
!
       if (ssu.ge.ssv) then
!4.7.3.1      case symr >= syms
!
       do 4731 ib=1,mapdb(0,5)
!
!4.7.3.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
!
!4.7.3.1.*      find porper A
       ia=mapia(symp,symq,ssu)
!
!4.7.3.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.7.3.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
       if (ssu.eq.ssv) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
!
       nhelp2=dimp*dimq
!
!4.7.3.1.*        def joined index and signum
       signum=1
       if (ssu.eq.ssv)then
       if (u.gt.v) then
       jjind=nshf(u)+v
       else if (u.eq.v) then
       jjind=0
       signum=0
       else
       jjind=nshf(v)+u
       signum=-1
       end if
       else
       jjind=(v-1)*dimr+u
       end if
!
!4.7.3.1.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,           &
     & signum)
!
 4731   continue
!
       else
!4.7.3.2      case symp < symq
!
       do 4732 ib=1,mapdb(0,5)
!
!4.7.3.2.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
!
!4.7.3.2.*      find porper A
       ia=mapia(symp,symq,ssv)
!
!4.7.3.2.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.7.3.2.*      def dimensions
       dimq=dimm(mapda(0,1),mapda(ia,3))
       dimp=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
       nhelp1=dimr*dims
!
       nhelp2=dimp*dimq
!
!4.7.3.2.*        def joined index and signum
       signum=-1
       jjind=(u-1)*dimr+v
!
!4.7.3.2.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,           &
     & signum)
!
 4732   continue
!
       end if
!
!
       else if (typa.eq.4) then
!
!4.7.4case A(pq,rs) -> B_rs(pq)
!
       typb=1
       call cct3_grc0 (2,typb,mapda(0,1),mapda(0,2),0,0,ssb,            &
     & possb0,posst,mapdb,mapib)
!
       if (ssu.ge.ssv) then
!4.7.4.1      case symr >= syms
!
       do 4741 ib=1,mapdb(0,5)
!
!4.7.4.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
!
!4.7.4.1.*      find porper A
       ia=mapia(symp,symq,ssu)
!
!4.7.4.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.7.4.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
       if (ssu.eq.ssv) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
!
       if (symp.eq.symq) then
       nhelp2=dimp*(dimp-1)/2
       else
       nhelp2=dimp*dimq
       end if
!
!4.7.4.1.*        def joined index and signum
       signum=1
       if (ssu.eq.ssv)then
       if (u.gt.v) then
       jjind=nshf(u)+v
       else if (u.eq.v) then
       jjind=0
       signum=0
       else
       jjind=nshf(v)+u
       signum=-1
       end if
       else
       jjind=(v-1)*dimr+u
       end if
!
!4.7.4.1.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,           &
     & signum)
!
 4741   continue
!
       else
!4.7.4.2      case symp < symq
!
       do 4742 ib=1,mapdb(0,5)
!
!4.7.4.2.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
!
!4.7.4.2.*      find porper A
       ia=mapia(symp,symq,ssv)
!
!4.7.4.2.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!4.7.4.2.*      def dimensions
       dimq=dimm(mapda(0,1),mapda(ia,3))
       dimp=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
!
       nhelp1=dimr*dims
!
       if (symp.eq.symq) then
       nhelp2=dimp*(dimp-1)/2
       else
       nhelp2=dimp*dimq
       end if
!
!4.7.4.2.*        def joined index and signum
       signum=-1
       jjind=(u-1)*dimr+v
!
!4.7.4.2.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,           &
     & signum)
!
 4742   continue
!
       end if
!
       else
!     RC=11 , nindA=4, exttyp=7, typa is nou 0,4 (NCI)
       rc=11
       return
       end if
!
       else
!     RC=12  nindA=4, exttyp 1,2,3,4,5 or 7 (Stup/NCI)
       rc=12
       return
!
       end if
!
       else if (nind.eq.2) then
!
       if (exttyp.eq.1) then
!
!2.1  case A(p q) -> A _p(q)
!
       if (typa.eq.0) then
!
!2.1.0case A(p,q) -> B _p(q)
!
       typb=0
       call cct3_grc0 (1,typb,mapda(0,2),0,0,0,ssb,                     &
     & possb0,posst,mapdb,mapib)
!
       do 210 ib=1,mapdb(0,5)
!
!2.1.0.*      def symmetry of all indices
       symq=mapdb(ib,3)
!
!2.1.0.*      find porper A
       ia=mapia(ssu,1,1)
!
!2.1.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!2.1.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
!
!2.1.0.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),dimp,dimq,u,1)
!
 210    continue
!
       else
!     RC=13  : nindA=2, exttyp=1, typa is not 0 (NCI)
       rc=13
       return
       end if
!
       else if (exttyp.eq.2) then
!
!2.2  case A(p q) -> A _q(p)
!
       if (typa.eq.0) then
!
!2.2.0case A(p,q) -> B _q(p)
!
       typb=0
       call cct3_grc0 (1,typb,mapda(0,1),0,0,0,ssb,                     &
     & possb0,posst,mapdb,mapib)
!
       do 220 ib=1,mapdb(0,5)
!
!2.2.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
!
!2.2.0.*      find porper A
       ia=mapia(symp,1,1)
!
!2.2.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
!
!2.2.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
!
!2.2.0.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),dimp,dimq,u,1)
!
 220    continue
!
!
       else
!     RC=14  : nindA=2, exttyp=2, typa is not 0 (NCI)
       rc=14
       return
       end if
!
       else
!     RC=15  : nindA=2, exttyp is not 1,2 (Stup)
       rc=15
       return
!
       end if
!
       else
!     RC=16  nindA is not 2 or 4 (NCI)
       rc=16
       return
!
       end if
!
       return
! Avoid unused argument warnings
      if (.false.) call Unused_integer(x)
       end
