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
       subroutine cct3_map (wrk,wrksize,                                &
     & nind,p,q,r,s,                                                    &
     & mapda,mapia,ssa,mapdb,mapib,possb0,posst,rc)
!
!     this routine realize mappings
!
!     B(indb) <-- A(inda)
!     where inda are order of indexes in mtx A and indb = Perm(inda)
!
!     nind   - number of indexes in matrix A (and B)  (Input)
!     p      - possition of 1-st index od mtx A in mtx B  (Input)
!     q      - possition of 2-nd index od mtx A in mtx B  (Input)
!     r      - possition of 3-rd index od mtx A in mtx B  (Input)
!     s      - possition of 4-th index od mtx A in mtx B  (Input)
!     mapda  - direct map matrix corresponding to A  (Input)
!     mapia  - inverse map matrix corresponding to A  (Input)
!     ssa    - overall symetry state  of matrix A  (Input)
!     mapdb  - direct map matrix corresponding to B  (Output)
!     mapib  - inverse map matrix corresponding to B  (Output)
!     possb0 - initial possition of matrix B in WRK  (Input)
!     posst  - final possition of matrix B in WRK (Output, not used yet)
!     rc     - return (error) code  (Output)
!
!
!     The table of implemented permutations
!
!     ninda  typA   p  q  r  s            Operation                 Implemented
!
!     4     0     all 24 comb.    A(1,2,3,4) -> B(p,q,r,s)             Yes
!
!     4     1     1  2  3  4      A(12,3,4)  -> B(12,3,4)              Yes
!     4     1     1  2  4  3      A(12,3,4)  -> B(12,4,3)              Yes
!     4     1     2  3  1  4      A(12,3,4)  -> B(3,12,4)              Yes
!     4     1           4  1      A(12,3,4)  -> B(4,12,3)              Yes
!     4     1     3  4  1  2      A(12,3,4)  -> B(3,4,12)              Yes
!     4     1           2  1      A(12,3,4)  -> B(4,3,12)              Yes
!     4     1     other comb.                                          No
!
!     4     2     1  2  3  4      A(1,23,4)  -> B(1,23,4)              Yes
!     4     2     4  2  3  1      A(1,23,4)  -> B(4,23,1)              Yes
!     4     2     3  1  2  4      A(1,23,4)  -> B(23,1,4)              Yes
!     4     2     4        3      A(1,23,4)  -> B(23,4,1)              Yes
!     4     2     1  3  4  2      A(1,23,4)  -> B(1,4,23)              Yes
!     4     2     2        1      A(1,23,4)  -> B(4,1,23)              Yes
!     4     2     other comb.                                          No
!
!     4     3     1  2  3  4      A(1,2,34)  -> B(1,2,34)              Yes
!     4     3     2  1  3  4      A(1,2,34)  -> B(2,1,34)              Yes
!     4     3     3  4  1  2      A(1,2,34)  -> B(34,1,2)              Yes
!     4     3     4  4            A(1,2,34)  -> B(34,2,1)              Yes
!     4     3     1  4  2  3      A(1,2,34)  -> B(1,34,2)              Yes
!     4     3     4  1            A(1,2,34)  -> B(2,34,1)              Yes
!     4     3     other comb.                                          No
!
!     4     4     1  2  3  4      A(12,34)   -> B(12,34)               Yes
!     4     4     3  4  1  2      A(12,34)   -> B(32,12)               Yes
!     4     4     other comb.                                          No
!
!     3     0     all 6 comb.     A(1,2,3)   -> B(p,q,r)               Yes
!
!     3     1     1  2  3  -      A(12,3)    -> B(12,3)                Yes
!     3     1     2  3  1  -      A(12,3)    -> B(3,12)                Yes
!     3     1     other comb.                                          No
!
!     3     2     1  2  3         A(1,23)    -> B(1,23)                Yes
!     3     2     3  1  2         A(1,23)    -> B(23,1)                Yes
!     3     2     other comb.                                          No
!
!     2     0     all 2 comb.     A(1,2)     -> B(p,q)                 Yes
!
!     2     1     1  2  -  -      A(12)      -> B(12)                  Yes
!     2     1     other comb.                                          No
!
!
#include "t31.fh"
#include "wrk.fh"

!
       integer nind,p,q,r,s,ssa,possb0,posst,rc
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
!
       integer type(1:4),dl(1:4),sa(1:4)
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7
       integer ia,ib,typ,newtyp
!
       rc=0
!
!     some test of correctness
!
       nhelp1=p+q+r+s
       if ((nind.eq.1).and.(nhelp1.eq.1)) then
       else if ((nind.eq.2).and.(nhelp1.eq.3)) then
       else if ((nind.eq.3).and.(nhelp1.eq.6)) then
       else if ((nind.eq.4).and.(nhelp1.eq.10)) then
       else
!     RC=1 - uncorrect indexes
       rc=1
       return
       end if
!
!     ******** No permutation *******
!
       if (nind.eq.1) then
       call cct3_noperm (wrk,wrksize,                                   &
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
!
       if ((nind.eq.2).and.(p.eq.1).and.(q.eq.2)) then
       call cct3_noperm (wrk,wrksize,                                   &
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
!
       if ((nind.eq.3).and.(p.eq.1).and.(q.eq.2).and.(r.eq.3)) then
       call cct3_noperm (wrk,wrksize,                                   &
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
!
       if ((nind.eq.4).and.(p.eq.1).and.(q.eq.2).and.(r.eq.3)           &
     &     .and.(s.eq.4)) then
       call cct3_noperm (wrk,wrksize,                                   &
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
!
       if (nind.eq.2) then
!
!     *********** 2 index ***********
!
       typ=mapda(0,6)
!
       if (typ.eq.0) then
!
!2.1  map A(p,q) -> B(q,p)
!
!     get mapdb,mapib
!
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       call cct3_grc0 (nind,0,type(1),type(2),0,0,ssa,                  &
     & possb0,posst,mapdb,mapib)
!
       do 210 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 210
!
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       ib=mapib(sa(1),1,1)
!
!     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
!
!     dimp,dimq
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,p,q,1)
!
 210    continue
!
       else
!     RC=2 - bad typ for: nind=2
       rc=2
       end if
!
       else if (nind.eq.3) then
!
!     *********** 3 index ***********
!
       typ=mapda(0,6)
!
       if (typ.eq.0) then
!
!3.1  map A(p,q,r) -> B(p1,q1,r1)
!
!     get mapdb,mapib
!
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call cct3_grc0 (nind,0,type(1),type(2),type(3),0,ssa,            &
     & possb0,posst,mapdb,mapib)
!
       do 310 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 310
!
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
!
!     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
!
!     dimp,dimq,r
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,          &
     &  nhelp5,p,q,r,1)
!
 310    continue
!
       else if (typ.eq.1) then
!
!3.2  map A(pq,r)
!     => only sophystical order of p,q,r is 2,3,1 i.e. B(r,pq)
!     => newtyp=2
!
       if ((p.eq.2).and.(q.eq.3).and.(r.eq.1)) then
!
!     get mapdb,mapib
!
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call cct3_grc0 (nind,2,type(1),type(2),type(3),0,ssa,            &
     & possb0,posst,mapdb,mapib)
!
       do 320 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 320
!
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
!
!     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
!
!     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
!
!     dimpq
       if (sa(2).eq.sa(3)) then
       nhelp7=nhelp3*(nhelp3-1)/2
       else
       nhelp7=nhelp3*nhelp4
       end if
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp7,nhelp5,2,1,1)
!
 320    continue
!
       else
!     RC=3 - bad order of p,q,r for: nind=3,typ=1
       rc=3
       end if
!
       else if (typ.eq.2) then
!
!3.3  map A(p,qr)
!     => only sophystical order of p,q,r is 3,1,2 i.e. B(qr,p)
!     => typb=1
!
       if ((p.eq.3).and.(q.eq.1).and.(r.eq.2)) then
!
!     get mapdb,mapib
!
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call cct3_grc0 (nind,1,type(1),type(2),type(3),0,ssa,            &
     & possb0,posst,mapdb,mapib)
!
       do 330 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 330
!
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
!
!     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
!
!     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
!
!     dimpq
       if (sa(1).eq.sa(2)) then
       nhelp7=nhelp4*(nhelp4-1)/2
       else
       nhelp7=nhelp4*nhelp5
       end if
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp7,2,1,1)
!
 330    continue
!
       else
!     RC=4 - bad order of p,q,r for: nind=3,typ=2
       rc=4
       end if
!
       else
!     RC=5 - bad type for: nind=3
       rc=5
       end if
!
       else if (nind.eq.4) then
!
!     *********** 4 index ***********
!
       typ=mapda(0,6)
!
       if (typ.eq.0) then
!
!4.1  map A(p,q,r,s)
!
!     get mapdb,mapib
!
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call cct3_grc0 (nind,0,type(1),type(2),type(3),type(4),ssa,      &
     & possb0,posst,mapdb,mapib)
!
       do 410 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 410
!
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
!
!     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
!
!     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       nhelp6=dimm(mapda(0,4),sa(s))
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,          &
     & nhelp5,nhelp6,p,q,r,s,1)
!
 410    continue
!
       else if ((typ.gt.0).and.(typ.lt.4)) then
!
!4.2  map A(p,q,r,s) with one limitation
!
!     def newtype and do corresponding tests
!
       if (typ.eq.1) then
       if ((p.gt.3).or.((q-p).ne.1)) then
!     RC=6 bad order of p,q,r,s, for: nind=4,typ=1
       rc=6
       return
       else
       nhelp1=p
       end if
       else if (typ.eq.2) then
       if ((q.gt.3).or.((r-q).ne.1)) then
!     RC=7 bad order of p,q,r,s, for: nind=4,typ=2
       rc=7
       return
       else
       nhelp1=q
       end if
       else if (typ.eq.3) then
       if ((r.gt.3).or.((s-r).ne.1)) then
!     RC=8 bad order of p,q,r,s, for: nind=4,typ=3
       rc=8
       return
       else
       nhelp1=r
       end if
       end if
       newtyp=nhelp1
!
!     get mapdb,mapib
!
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call cct3_grc0 (nind,newtyp,type(1),type(2),type(3),type(4),ssa, &
     & possb0,posst,mapdb,mapib)

       do 420 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 420
!
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
!
!     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
!
!     dimp,dimq,dimr,dims
       dl(p)=dimm(mapda(0,1),sa(p))
       dl(q)=dimm(mapda(0,2),sa(q))
       dl(r)=dimm(mapda(0,3),sa(r))
       dl(s)=dimm(mapda(0,4),sa(s))
!
       if ((newtyp.eq.1).and.(sa(1).eq.sa(2))) then
!     B(p1q1,r1,s1) case => oldtyp can be 2,3
       nhelp7=dl(1)*(dl(1)-1)/2
!
       if (typ.eq.1) then
!     A(pq,r,s) case
!     A(pq,r,s) -> A(pq,s,r)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),           &
     &  dl(s),1,3,2,1)
       else if (typ.eq.2) then
!     A(p,qr,s) case
       if (p.eq.3) then
!     A (p,qr,s) -> B(qr,p,s)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),     &
     &  2,1,3,1)
       else
!     A (p,qr,s) -> B(qr,s,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),     &
     &  3,1,2,1)
       end if
       else
!     A(p,q,rs) case
       if (p.eq.3) then
!     A (p,q,rs) -> B(rs,p,q)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),            &
     &  nhelp7,2,3,1,1)
       else
!     A (p,q,rs) -> B(rs,q,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),            &
     &  nhelp7,3,2,1,1)
       end if
       end if
!
       else if ((newtyp.eq.2).and.(sa(2).eq.sa(3))) then
!     B(p1,q1r1,s1) case => oldtyp can be 1,3
       nhelp7=dl(2)*(dl(2)-1)/2
       if (typ.eq.1) then
!     A(pq,r,s) case
       if (r.eq.1) then
!     A (pq,r,s) -> B(r,pq,s)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),     &
     &  2,1,3,1)
       else
!     A (pq,r,s) -> B(s,pq,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),     &
     &  2,3,1,1)
       end if
       else if (typ.eq.2) then
!     A(p,qr,s) case
!     A(p,qr,s) -> B(s,qr,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),     &
     & 3,2,1,1)
       else
!     A(p,q,rs) case
       if (p.eq.1) then
!     A (p,q,rs) -> B(p,rs,q)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,     &
     & 1,3,2,1)
       else
!     A (p,q,rs) -> B(q,rs,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,     &
     & 3,1,2,1)
       end if
       end if
!
       else if ((newtyp.eq.3).and.(sa(3).eq.sa(4))) then
!     B(p1,q1,r1s1) case => oldtyp can be 1,2
       nhelp7=dl(3)*(dl(3)-1)/2
       if (typ.eq.1) then
!     A(pq,r,s) case
       if (r.eq.2) then
!     A (pq,r,s) -> B(r,s,pq)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),     &
     & 3,1,2,1)
       else
!     A (pq,r,s) -> B(s,r,pq)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),     &
     & 3,2,1,1)
       end if
       else if (typ.eq.2) then
!     A(p,qr,s) case
       if (p.eq.1) then
!     A (p,qr,s) -> B(p,s,qr)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),     &
     & 1,3,2,1)
       else
!     A (p,qr,s) -> B(s,p,qr)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),     &
     & 2,3,1,1)
       end if
       else
!     A(p,q,rs) case
!     A(p,q,rs) -> A(q,p,rs)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,     &
     &  2,1,3,1)
       end if
!
       else
!     B(p1,q1,r1,s1) case
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),dl(r),      &
     & dl(s),p,q,r,s,1)
       end if
!
 420    continue
!
!
       else if (typ.eq.4) then
!
!4.2  map A(pq,rs) -> B(rs,pq)
!
!     tests
!
       if ((p.ne.3).and.(q.ne.4).and.(r.ne.1).and.(s.ne.2)) then
!     RC=9 : nind=4, typA=4 (bad order of p,q,r,s, Stup)
       rc=9
       return
       end if
!
!     get mapdb,mapib
!
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call cct3_grc0 (nind,4,type(1),type(2),type(3),type(4),ssa,      &
     & possb0,posst,mapdb,mapib)
!
       do 430 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 430
!
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
!
!     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
!
!     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       nhelp6=dimm(mapda(0,4),sa(s))
!
       if ((sa(p).eq.sa(q)).and.(sa(r).eq.sa(s))) then
!     A(pq,rs) -> B(rs,pq)
       nhelp3=nhelp3*(nhelp3-1)/2
       nhelp5=nhelp5*(nhelp5-1)/2
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,2,1,1)
!
       else if (sa(p).eq.sa(q)) then
!     A(pq,r,s) -> B(r,s,pq)
       nhelp3=nhelp3*(nhelp3-1)/2
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,nhelp6,   &
     & 3,1,2,1)
!
       else if (sa(r).eq.sa(s)) then
!     A(p,q,rs) -> B(rs,p,q)
       nhelp5=nhelp5*(nhelp5-1)/2
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,   &
     & 2,3,1,1)
!
       else
!     A(p,q,r,s) -> B(r,s,p,q)
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,   &
     & nhelp6,3,4,1,2,1)
       end if
!
 430    continue
!
       else
!     RC=10 bad type for: nind=4
       rc=10
       end if
!
       else
!
!     *********  more than 4 indexes
!
!     RC=11 bad number of indexes
       rc=11
       end if
!
!
       return
       end
