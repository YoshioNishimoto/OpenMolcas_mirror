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
       subroutine cct3_add (wrk,wrksize,                                &
     & ninda,nindb,nindext,typext,u,v,ssu,ssv,factor,                   &
     & mapda,ssa,mapdb,mapib,ssb,rc)
!
!     this routine do:
!     B(indb) = B(indb) + factor * A(inda)
!
!     ninda   - # of indexes in matrix A (1-4)
!     nindb   - # of indexes in matrix B (1-4)
!     nindext - # of external (frozen, fixed) indexes (now:0-2)
!     typext  - characterize external indexes as follows:
!     0 - no frozen index
!     1 - frozen index p
!     2 - frozen index q
!     3 - frozen index r
!     4 - frozen index s
!     5 - frozen indexes p,q
!     6 - frozen indexes r,s
!     u       - value of first external index (if any, else 0)
!     v       - value of second external index (if any, else 0)
!     ssu     - symmetry of u (if any, else 1)
!     ssv     - symmetry of v (if any, else 1)
!     factor  - multiplicative factor (see def)
!     mapda   - direct map matrix corresponding to A (see docc.txt)
!     ssa     - overall spin state of matrix A
!     mapdb   - direct map matrix corresponding to B (see docc.txt)
!     mapib   - inverse map matrix corresponding to B (see docc.txt)
!     ssb     - overall spin state of matrix B
!     rc      - return (error) code
!
!     Table of present implementations:
!
!     nindB  nindxet  typext  typB  =>   Implemented
!     >4                                   No
!
!     4       0        0    0-4            Yes
!     4       1       1-4   0,4            Yes
!     4       1       1-4   2,3            No
!     4       2        5    0,4            Yes
!     4       2        5    2,3            No
!     4       2       6-n                  No
!     4       3                            No
!
!     3       0       0     0-2            Yes
!     3       1       1-3    0             Yes
!     3       1       1-3   1,2            No
!     3       2                            No
!
!     2       0       0     0,1            Yes
!     2       1       1-2    0             Yes
!     2       1       1-2    1             No
!     2       2                            No
!
!     1                                    No
!
!
!     !N.B. oprav co je oznacene c@!
!
#include "t31.fh"
#include "wrk.fh"
!
       integer ninda,nindb,nindext,typext,u,v,ssu,ssv,ssa,ssb,rc
       real*8 factor
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
!
       integer mapib(1:8,1:8,1:8)
!
!     help variables
!
       integer sa1,sa2,sa3,ssp,ssq,pq
       integer :: nhelp1=0,nhelp2=0,nhelp3=0,nhelp4=0,nhelp5=0
       integer :: nhelp6=0,nhelp7=0,nhelp8=0,nhelp9=0,nhelp10=0
       integer ia,ib,ibm
       integer typa,typb,p,q
       real*8 fact
!      To fix some 'uninitialized' warnings
       p=0
       q=0
!     general tests
!
       nhelp1=nindA+nindext
!
       if (nhelp1.ne.nindb) then
!     RC=1  : incompatible (nindA, nindB and nindext, Stup)
       rc=1
       return
       end if
!
       nhelp1=mmul(ssu,ssv)
       nhelp1=mmul(ssa,nhelp1)
       if (nhelp1.ne.ssb) then
!     RC=2  : incompatible (ssa, ssb ,ssu and ssv, Stup)
       rc=2
       return
       end if
!
       typa=mapda(0,6)
       typb=mapdb(0,6)
       fact=factor
!
       if (nindext.gt.0) then
       if ((typb.ge.1).and.(typb.le.3)) then
!     RC=3 : nindext>0, typB is 1,2 or 3 (NCI)
       rc=3
       return
       end if
       end if

!
       if ((nindext.eq.2).and.(typb.eq.4)) then
!
!     def p,q,ssp,ssq,fact(new)
!
       if (ssu.gt.ssv) then
!     ssu>ssv
       p=u
       q=v
       ssp=ssu
       ssq=ssv
       fact=factor
       else if (ssu.eq.ssv) then
!     ssu=ssv
       if (u.ge.v) then
       p=u
       q=v
       ssp=ssu
       ssq=ssv
       fact=factor
       else
       p=v
       q=u
       ssp=ssv
       ssq=ssu
       fact=-factor
       end if
       else
!     ssu<ssv
       p=v
       q=u
       ssp=ssv
       ssq=ssu
       fact=-factor
       end if
!
       end if
!
       if (nindb.eq.4) then
!
!     **********  -> B(pqrs) **********
!
       if (nindext.eq.0) then
!
!400  case B(pqrs) <-- A(pqrs) or B(p,q,r,s) <-- A(p,q,r,s)
!
!     tests
!
       if (typa.ne.typb) then
!     RC=4 : nindB=4, nindext=0 (TypA incompatible with TypB ,Stup)
       rc=4
       return
       end if
!
       do 400 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,sa2,sa3)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 400
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
       call cct3_add10 (wrk(nhelp2),wrk(nhelp3),nhelp1,fact)
!
 400    continue
!
       else if (nindext.eq.1) then
!
       if (typext.eq.1) then
!
       if (typb.eq.0) then
!
!4110 case B(p,q,r,s) <-- A(q,r,s)
!
!     tsets
!
       if (typa.ne.0) then
!     RC=5 : nindB=4, nindeext=1, typext=1, typB=0 (typA is not 0, Stup)
       rc=5
       return
       end if
!
       do 4110 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(ssu,sa1,sa2)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4110
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       nhelp8=nhelp5*nhelp6*nhelp7
!
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp8,fact)
!
 4110   continue
!
       else if (typb.eq.4) then
!
!4114 case B(pq,rs) <-- A(q,rs)
!
!     tsets
!
       if (typa.ne.2) then
!     RC=6  : nindB=4, nindeext=1, typext=1, typB=4 (typA is not 2, Stup)
       rc=6
       return
       end if
!
       do 4114 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(ssu,sa1,sa2)
       ibm=mapib(sa1,ssu,sa2)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4114
!
!     def possA
       nhelp2=mapda(ia,1)
!
!
       if (ssu.gt.sa1) then
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
!
!     def possB
       nhelp3=mapdb(ib,1)
!     def fictive dimensions
       nhelp9=nhelp5*nhelp8
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,fact)
!
       else if (ssu.eq.sa1) then
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
!
!     def possB
       nhelp3=mapdb(ib,1)
!     def fictive dimensions
       nhelp9=nhelp4*(nhelp4-1)/2
       call cct3_add41 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp4,nhelp9,nhelp8,fact)
!
       else
!     ssu<sa1  B(qp,rs) <-- -A_p (q,rs)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ibm,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ibm,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ibm,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ibm,6))
!
!     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
!
!     def possB-
       nhelp3=mapdb(ibm,1)
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp4,nhelp5,nhelp8,-fact)
!
       end if
!
 4114   continue
!
       else
!     RC=7 : nindB=4, nindext=1, typext=1, (typA is not 0 or 4, (NCI))
       rc=7
       return
       end if
!
       else if (typext.eq.2) then
!
       if (typb.eq.0) then
!
!4120 case B(p,q,r,s) <-- A(p,r,s)
!
!     tsets
!
       if (typa.ne.0) then
!     RC=8 : nindB=4, nindeext=1, typext=2, typB=0 (typA is not 0, Stup)
       rc=8
       return
       end if
!
       do 4120 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,ssu,sa2)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4120
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       nhelp8=nhelp6*nhelp7
!
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp4,nhelp5,nhelp8,fact)
!
 4120   continue
!
       else if (typb.eq.4) then
!
!4124 case B(pq,rs) <-- A(p,rs)
!
!     tsets
!
       if (typa.ne.2) then
!     RC=9 : nindB=4, nindeext=1, typext=2, typB=4 (typA is not 2, Stup)
       rc=9
       return
       end if
!
       do 4124 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,ssu,sa2)
       ibm=mapib(ssu,sa1,sa2)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4124
!
!     def possA
       nhelp2=mapda(ia,1)
!
!
       if (sa1.gt.ssu) then
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
!
!     def possB
       nhelp3=mapdb(ib,1)
       call cct3_add32 (wrk(nhelp2),                                    &
     &  wrk(nhelp3),u,nhelp4,nhelp5,nhelp8,fact)
!
!
       else if (sa1.eq.ssu) then
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
!
!     def possB
       nhelp3=mapdb(ib,1)
!     def fictive dimensions
       nhelp9=nhelp4*(nhelp4-1)/2
       call cct3_add42 (wrk(nhelp2),wrk(nhelp3),                        &
     & u,nhelp4,nhelp9,nhelp8,fact)
!
!
       else
!     sa1<ssu  B(qp,rs) <-- -A_q (p,rs)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ibm,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ibm,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ibm,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ibm,6))
!
!     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
!
!     def possB-
       nhelp3=mapdb(ibm,1)
!     def fictive index
       nhelp9=nhelp8*nhelp5
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,-fact)
!
       end if
!
 4124   continue
!
       else
!     RC=10: nindB=4, nindext=1, typext=2, (typA is not 0 or 4, (NCI))
       rc=10
       return
       end if
!
       else if (typext.eq.3) then
!
       if (typb.eq.0) then
!
!4130 case B(p,q,r,s) <-- A(p,q,s)
!
!     tsets
!
       if (typa.ne.0) then
!     RC=11: nindB=4, nindeext=1, typext=3, typB=0 (typA is not 0, Stup)
       rc=11
       return
       end if
!
       do 4130 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,sa2,ssu)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4130
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       nhelp8=nhelp4*nhelp5
!
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp8,nhelp6,nhelp7,fact)
!
 4130   continue
!
       else if (typb.eq.4) then
!
!4134 case B(pq,rs) <-- A(pq,s)
!@!   oprav to tak ako v typext 1 a 2
!
!     tsets
!
       if (typa.ne.1) then
!     RC=12: nindB=4, nindeext=1, typext=3, typB=4 (typA is not 1, Stup)
       rc=12
       return
       end if
!
       do 4134 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,sa2,ssu)
       ibm=mapib(sa1,sa2,sa3)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4134
!
!     def possA
       nhelp2=mapda(ia,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       if (sa1.eq.sa2) then
       nhelp8=nhelp4*(nhelp4-1)/2
       else
       nhelp8=nhelp4*nhelp5
       end if
!
       if (ssu.gt.sa3) then
!
!     def possB
       nhelp3=mapdb(ib,1)
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp8,nhelp6,nhelp7,fact)
!
       else if (ssu.eq.sa3) then
!
!     def possB
       nhelp3=mapdb(ib,1)
!     def fictive dimensions
       nhelp9=nhelp6*(nhelp6-1)/2
       call cct3_add43 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp8,nhelp9,nhelp6,fact)
!
       else
!     ssu<sa3  B(pq,sr) <-- -A_r (pq,s)
!     def possB-
       nhelp3=mapdb(ibm,1)
!     def fictive dimension
       nhelp9=nhelp8*nhelp7
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp9,nhelp6,-fact)
!
       end if
!
 4134   continue
!
       else
!     RC=13: nindB=4, nindext=1, typext=3, (typA is not 0 or 4, (NCI))
       rc=13
       return
       end if
!
       else if (typext.eq.4) then
!
       if (typb.eq.0) then
!
!4140 case B(p,q,r,s) <-- A(p,q,r)
!
!     tsets
!
       if (typa.ne.0) then
!     RC=14: nindB=4, nindeext=1, typext=4, typB=0 (typA is not 0, Stup)
       rc=14
       return
       end if
!
       do 4140 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,sa2,sa3)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4140
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       nhelp8=nhelp4*nhelp5*nhelp6
!
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp7,fact)
!
 4140   continue
!

       else if (typb.eq.4) then
!
!4144 case B(pq,rs) <-- A(pq,r)
!@!   oprav to tak ako v typext 1 a 2
!
!     tsets
!
       if (typa.ne.1) then
!     RC=15: nindB=4, nindeext=1, typext=4, typB=4 (typA is not 1, Stup)
       rc=15
       return
       end if
!
       do 4144 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,sa2,sa3)
       ibm=mapib(sa1,sa2,ssu)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4144
!
!     def possA
       nhelp2=mapda(ia,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     def fictive dimensions
       if (sa1.eq.sa2) then
       nhelp8=nhelp4*(nhelp4-1)/2
       else
       nhelp8=nhelp4*nhelp5
       end if
!
       if (sa3.gt.ssu) then
!
!     def possB
       nhelp3=mapdb(ib,1)
!     def fictive dimension
       nhelp9=nhelp8*nhelp6
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp7,fact)
!
       else if (sa3.eq.ssu) then
!
!     def possB
       nhelp3=mapdb(ib,1)
!     def fictive dimensions
       nhelp9=nhelp6*(nhelp6-1)/2
       call cct3_add44 (wrk(nhelp2),                                    &
     &  wrk(nhelp3),u,nhelp8,nhelp9,nhelp6,fact)
!
       else
!     sa3<ssu  B(pq,sr) <-- -A_s (pq,r)
!     def possB-
       nhelp3=mapdb(ibm,1)
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp8,nhelp7,nhelp6,-fact)
!
       end if
!
 4144   continue
!
       else
!     RC=16: nindB=4, nindext=1, typext=4, (typA is not 0 or 4, (NCI))
       rc=16
       return
       end if
!
       else
!     RC=17: nindB=4, nindext=1, typext=@  (Stup)
       rc=17
       return
       end if
!
       else if (nindext.eq.2) then
!
       if (typext.eq.5) then
!
       if (typb.eq.0) then
!
!4250 case B(p,q,r,s) <-- A(r,s)
!
!     tests
!
       if ((typb.eq.0).and.(typa.ne.0)) then
!     RC=18: nindB=4, nindext=2, typext=5, typB=0 (typA is not 0, Stup)
       rc=18
       return
       end if
!
       do 4250 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
!@
       nhelp1=mmul(ssu,ssv)
       nhelp1=mmul(nhelp1,sa1)
       nhelp1=mmul(nhelp1,ssb)
       if (nhelp1.ne.sa2) then
       write(6,*) ' Add Bpqrs <- Ars incorrect',ssp,ssq,sa1,nhelp1,sa1, &
     & sa2
       goto 4250
       end if
!@@
!
       ib=mapib(ssu,ssv,sa1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4250
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     calc joined pq index
       pq=(v-1)*nhelp4+u
!
!     calc fictive lengths
       nhelp9=nhelp6*nhelp7
       nhelp10=nhelp4*nhelp5
!
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),pq,nhelp10,nhelp9,fact)
!
 4250   continue
!
       else if (typb.eq.4) then
!
!4254 case B(pq,rs) <-- A(rs)
!
!     tests
!
       if ((typb.eq.4).and.(typa.ne.1)) then
!     RC=19: nindB=4, nindext=2, typext=5, typB=4 (typA is not 1, Stup)
       rc=19
       return
       end if
!
       do 4254 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
!
!@
       nhelp1=mmul(ssp,ssq)
       nhelp1=mmul(nhelp1,sa1)
       nhelp1=mmul(nhelp1,ssb)
       if (nhelp1.ne.sa2) then
       write(6,*) ' Add Bpqrs <- Ars incorrect'
       goto 4254
       end if
!@@
       ib=mapib(ssp,ssq,sa1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4254
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
!
!     calc joined pq index and fictive length of pq pair
       if (ssp.eq.ssq) then
       pq=(p-1)*(p-2)/2+q
       nhelp10=nhelp4*(nhelp4-1)/2
       else
       pq=(q-1)*nhelp4+p
       nhelp10=nhelp4*nhelp5
       end if
!
!     calc fictive lengths
       if (sa1.eq.sa2) then
       nhelp9=nhelp6*(nhelp6-1)/2
       else
       nhelp9=nhelp6*nhelp7
       end if
!
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),pq,nhelp10,nhelp9,fact)
!
 4254   continue
!
       else
!     RC=20: nindB=4, nindext=2, typext=5 (typB is not 0 or 4, NCI)
       rc=20
       return
       end if
!
       else if (typext.eq.6) then
!
!426  case B(p,q,r,s) <-- A(p,q) and B(pq,rs) <-- A(pq)
!
!     RC=21: nindB=4, nindext=2, typext=6, NCI)
       rc=21
       return
!
       else
!     RC=22: nindB=4, nindext=2, (typext is not 5 or 6, NCI)
       rc=22
       return
       end if
!
       else
!     RC=23: nindB=4, nindext>2 (NCI)
       rc=23
       return
       end if
!
       else if (nindb.eq.3) then
!
!     **********  -> B(pqr) **********
!
       if (nindext.eq.0) then
!
!300  case B(pqr) <-- A(pqr)
!
!     tests
!
       if (typa.ne.typb) then
!     RC=24: nindB=3, nindext=0 (TypA incompatible with TypB ,Stup)
       rc=24
       return
       end if
!
       do 300 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
!
       ib=mapib(sa1,sa2,1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 300
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
       call cct3_add10 (wrk(nhelp2),wrk(nhelp3),nhelp1,fact)
!
 300    continue
!
       else if (nindext.eq.1) then
!
       if (typext.eq.1) then
!
!311  case B(p,q,r) <-- A(q,r)
!
       if ((typa.eq.0).and.(typb.eq.0)) then
!
!311  case B(p,q,r) <-- A(q,r)
!
       do 311 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
!
       ib=mapib(ssu,sa1,1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 311
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
!
!     def fictive dimensions
       nhelp7=nhelp5*nhelp6
!
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp7,fact)
!
 311    continue
!
       else
!     RC=25: nindB=3, nindext=1, typext=1 (tybA,B is not 0, NCI)
       rc=25
       return
       end if
!
       else if (typext.eq.2) then
!
       if ((typa.eq.0).and.(typb.eq.0)) then
!
!312  case B(p,q,r) <-- A(p,r)
!
       do 312 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
!
       ib=mapib(sa1,ssu,1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 312
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
!
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),                        &
     &  u,nhelp4,nhelp6,nhelp7,fact)
!
 312    continue
!
       else
!     RC=26: nindB=3, nindext=1, typext=2 (tybA,B is not 0, NCI)
       rc=26
       return
       end if
!
       else if (typext.eq.3) then
!
       if ((typa.eq.0).and.(typb.eq.0)) then
!
!313  case B(p,q,r) <-- A(p,q)
!
       do 313 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
!
       ib=mapib(sa1,sa2,1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 313
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq,dimr
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
!
!     def fictive dimensions
       nhelp7=nhelp4*nhelp5
!
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp7,nhelp6,fact)
!
 313    continue
!
       else
!     RC=27: nindB=3, nindext=1, typext=3 (tybA,B is not 0, NCI)
       rc=27
       return
       end if
!
       else
!     RC=28: nindB=3 , typext=@ (Stup)
       rc=28
       return
       end if
!
       else
!     RC=29: nindB=3, nindext>1 (NCI)
       rc=29
       return
       end if
!
       else if (nindb.eq.2) then
!
!     **********  -> B(pq) **********
!
       if (nindext.eq.0) then
!
!200  case B(p,q) <-- A(p,q)
!
       do 200 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
!
       ib=mapib(sa1,1,1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 200
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
       call cct3_add10 (wrk(nhelp2),wrk(nhelp3),nhelp1,fact)
!
 200    continue
!
       else if (nindext.eq.1) then
!
       if (typext.eq.1) then
!
!211  case B(p,q) <-- A(q)
!
       if ((typa.eq.0).and.(typb.eq.0)) then
!
       do 211 ia=1,mapda(0,5)
!
       ib=mapib(ssu,1,1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 211
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
!
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,fact)
!
 211    continue
!
       else
!     RC=30: nindB=2, nindext=1, typext=1, (typA,B is not 0, NCI)
       rc=30
       return
       end if
!
       else if (typext.eq.2) then
!
!212  case B(p,q) <-- A(p)
!
       if ((typa.eq.0).and.(typb.eq.0)) then
!
       do 212 ia=1,mapda(0,5)
!
       sa1=mapda(ia,3)
       ib=mapib(sa1,1,1)
!
!     def length
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 212
!
!     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
!
!     def dimp,dimq
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
!
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,fact)
!
 212    continue
!
       else
!     RC=31: nindB=2, nindext=1, typext=2, (typA,B is not 0, NCI)
       rc=31
       return
       end if
!
       else
!     RC=32: nindB=2, nindext=1, typext=@ (Stup)
       rc=32
       return
       end if
!
       else
!     RC=33: nindB=2, ininext>1 (NCI)
       rc=33
       return
       end if
!
       else
!     RC=34: nindb less then 2 (NCI/Stup)
       rc=34
       return
       end if
!
       return
       end
