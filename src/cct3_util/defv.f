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
       subroutine defv (wrk,wrksize,                                    &
     & deftyp,possv0,mapdv,mapiv,ssv,                                   &
     & mapdr,mapir,ssr,rc)
!
!     this routine define v mediate as:
!     V_i(abc) = <ab||ic>  from integrals, stored in R_i(abc)
!     R _i(a,bc)bbb is a matrix, where integrals <ab|ic> are stored
!     for b>=c
!
!     deftyp - type of definition (see table) (I)
!     possv0 - initial address of V (I)
!     mapdv  - direct map of V (O)
!     mapiv  - inverse map of V (O)
!     ssv    - overall spin of V (O)
!     mapdr  - direct map of R (I)
!     mapir  - inverse map of R (I)
!     ssr    - overall spin of R (I)
!     rc     - return (error) code (O)
!
!     Table of definitions
!
!     deftyp         Operation                 Implement.
!     1      V(ab,c)aaa  = R(abc)-R(bac)         Yes
!     2      V(ab,c)bbb  = R(abc)-R(bac)         Yes
!     3      V(a,b,c)abb = R(abc)                Yes
!     4      V(a,b,c)aba = -R(bac)                Yes
!
#include "t31.fh"
#include "wrk.fh"
!
       integer deftyp,possv0,ssv,ssr,rc
       integer mapdv(0:512,1:6),mapdr(0:512,1:6)
       integer mapiv(1:8,1:8,1:8),mapir(1:8,1:8,1:8)
!
!     help variables
!
       integer posst,possv,possr1,possr2
       integer iv,ir1,ir2,syma,symb,symc
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5
       integer nhelp6,nhelp7,nhelp8,nhelp9,nhelp10
!
!
!0.*  def mapdv,mapiv
       if (deftyp.eq.1) then
!     case V(ab,c)aaa
       call cct3_grc0 (3,1,3,3,3,0,ssr,possv0,posst,mapdv,mapiv)
       else if (deftyp.eq.2) then
!     case V(ab,c)bbb
       call cct3_grc0 (3,1,4,4,4,0,ssr,possv0,posst,mapdv,mapiv)
       else if (deftyp.eq.3) then
!     case V(a,b,c)abb
       call cct3_grc0 (3,0,3,4,4,0,ssr,possv0,posst,mapdv,mapiv)
       else if (deftyp.eq.4) then
!     case V(a,b,c)aba
       call cct3_grc0 (3,0,3,4,3,0,ssr,possv0,posst,mapdv,mapiv)
       else
!     RC=1 , deftyp out of range (1-4) (Stup)
       rc=1
       return
       end if
!
!0.*  some tests
!

!0.*  define spin
       ssv=ssr
!
       if ((deftyp.eq.1).or.(deftyp.eq.2)) then
!
!12   case V(ab,c)aaa,bbb
!
       do 12 iv=1,mapdv(0,5)
!
!12.* def symmetries
       syma=mapdv(iv,3)
       symb=mapdv(iv,4)
       symc=mapdv(iv,5)
!
       if (syma.eq.symb) then
!12.1 syma=symb
!
       if (symb.eq.symc) then
!12.1.1            case syma=symb=symc
!
!12.1.1.*     def possitions of V, R1
       possv=mapdv(iv,1)
       ir1=mapir(syma,symb,1)
       possr1=mapdr(ir1,1)
!
!12.1.1.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nhelp1*(nhelp1+1)/2
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp3=nhelp4*(nhelp4-1)/2
       nhelp5=nvb(syma)-nhelp4
!
!12.1.1.*     realize definition of V
       call defvhlp1 (wrk(possr1),wrk(possv),                           &
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5)
!
       else
!12.1.2            case syma=symb.ne.symc
!
!12.1.2.*     def possitions of V, R1
       possv=mapdv(iv,1)
       if (syma.ge.symc) then
       ir1=mapir(syma,symb,1)
       else
       ir1=mapir(syma,symc,1)
       end if
       possr1=mapdr(ir1,1)
!
!12.1.2.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp3=nhelp4*(nhelp4-1)/2
       nhelp5=dimm(mapdv(0,3),symc)
       nhelp6=nvb(syma)-nhelp4
       nhelp7=nvb(symc)-nhelp5
!
!12.1.2.*     realize definition of V
       if (syma.ge.symc) then
       call defvhlp21 (wrk(possr1),wrk(possv),nhelp1,nhelp2,            &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       else
       call defvhlp22 (wrk(possr1),wrk(possv),nhelp1,nhelp2,            &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       end if
!
       end if
!
       else
!12.2 syma > symb
!
       if (syma.eq.symc) then
!12.2.1     case syma>symb , syma=symc
!
!12.2.1.*     def possitions of V, R1, R2
       possv=mapdv(iv,1)
!     R1 is permuted, since (a=c)>b
       ir1=mapir(syma,symc,1)
       possr1=mapdr(ir1,1)
       ir2=mapir(symb,syma,1)
       possr2=mapdr(ir2,1)
!
!12.2.1.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nvb(symc)
       nhelp4=nhelp1*(nhelp1+1)/2
       nhelp5=dimm(mapdv(0,1),syma)
       nhelp6=dimm(mapdv(0,2),symb)
       nhelp7=dimm(mapdv(0,3),symc)
       nhelp8=nvb(syma)-nhelp5
       nhelp9=nvb(symb)-nhelp6
       nhelp10=nvb(symc)-nhelp7
!
!12.2.1.*     realize definition of V
       call defvhlp3 (wrk(possr1),wrk(possr2),wrk(possv),               &
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,                              &
     & nhelp6,nhelp7,nhelp8,nhelp9,nhelp10)
!
       else if (symb.eq.symc) then
!12.2.2     case syma>symb , symb=symc
!
!12.2.2.*     def possitions of V, R1, R2
       possv=mapdv(iv,1)
       ir1=mapir(syma,symb,1)
       possr1=mapdr(ir1,1)
       ir2=mapir(symb,syma,1)
       possr2=mapdr(ir2,1)
!
!12.2.2.*     def dimensions
       nhelp1=nvb(syma)
       nhelp3=nvb(symb)
       nhelp4=nvb(symc)
       nhelp2=nhelp3*(nhelp3+1)/2
       nhelp5=dimm(mapdv(0,1),syma)
       nhelp6=dimm(mapdv(0,2),symb)
       nhelp7=dimm(mapdv(0,3),symc)
       nhelp8=nvb(syma)-nhelp5
       nhelp9=nvb(symb)-nhelp6
       nhelp10=nvb(symc)-nhelp7
!
!12.2.2.*     realize definition of V
       call defvhlp4 (wrk(possr1),wrk(possr2),wrk(possv),               &
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,                              &
     & nhelp6,nhelp7,nhelp8,nhelp9,nhelp10)
!
       else
!12.2.3     case syma>symb , symc.ne.syma,symb
!
!
!12.2.3.*     def possitions of V, R1, R2
       possv=mapdv(iv,1)
!
       if (symb.ge.symc) then
       ir1=mapir(syma,symb,1)
       else
       ir1=mapir(syma,symc,1)
       end if
       possr1=mapdr(ir1,1)
!
       if (syma.ge.symc) then
       ir2=mapir(symb,syma,1)
       else
       ir2=mapir(symb,symc,1)
       end if
       possr2=mapdr(ir2,1)
!
!12.2.3.*     def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
       nhelp8=nvb(symb)-nhelp5
       nhelp9=nvb(symc)-nhelp6
!
!12.2.3.*     realize definition of V
       if ((syma.gt.symc).and.(symb.gt.symc)) then
       call defvhlp51 (wrk(possr1),wrk(possr2),wrk(possv),              &
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,                              &
     & nhelp6,nhelp7,nhelp8,nhelp9)
       else if ((syma.gt.symc).and.(symb.lt.symc)) then
       call defvhlp52 (wrk(possr1),wrk(possr2),wrk(possv),              &
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,                              &
     & nhelp6,nhelp7,nhelp8,nhelp9)
       else if ((syma.lt.symc).and.(symb.gt.symc)) then
       call defvhlp53 (wrk(possr1),wrk(possr2),wrk(possv),              &
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,                              &
     & nhelp6,nhelp7,nhelp8,nhelp9)
       else if ((syma.lt.symc).and.(symb.lt.symc)) then
       call defvhlp54 (wrk(possr1),wrk(possr2),wrk(possv),              &
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,                              &
     & nhelp6,nhelp7,nhelp8,nhelp9)
       end if
!
       end if
!
       end if
!
 12     continue
!
!
       else if (deftyp.eq.3) then
!
!3    case V(a,b,c)abb
!
       do 3 iv=1,mapdv(0,5)
!
!3.*  def symmetries
       syma=mapdv(iv,3)
       symb=mapdv(iv,4)
       symc=mapdv(iv,5)
!
       if (symb.eq.symc) then
!3.1  symb=symc
!
!3.1.*def possitions of V, R1
       possv=mapdv(iv,1)
       ir1=mapir(syma,symb,1)
       possr1=mapdr(ir1,1)
!
!3.1.*def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nhelp2*(nhelp2+1)/2
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4

!3.1.*realize definition of V
       call defvhlp7 (wrk(possr1),wrk(possv),nhelp1,nhelp2,             &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
!
       else
!     case symb.ne.symc
!
!3.1.*def possitions of V, R1
       possv=mapdv(iv,1)
       if (symb.ge.symc) then
       ir1=mapir(syma,symb,1)
       else
       ir1=mapir(syma,symc,1)
       end if
       possr1=mapdr(ir1,1)
!
!3.1.*def dimensions
       nhelp1=nvb(syma)
       nhelp2=nvb(symb)
       nhelp3=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
!
!3.1.*realize definition of V
       if (symb.gt.symc) then
       call defvhlp61 (wrk(possr1),wrk(possv),nhelp1,nhelp2,            &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       else
       call defvhlp62 (wrk(possr1),wrk(possv),nhelp1,nhelp2,            &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
       end if
!
       end if
!
 3      continue
!
       else if (deftyp.eq.4) then
!
!4    case V(a,b,c)aba
!
       do 4 iv=1,mapdv(0,5)
!
!4.*  def symmetries
       syma=mapdv(iv,3)
       symb=mapdv(iv,4)
       symc=mapdv(iv,5)
!
       if (syma.eq.symc) then
!4.1  syma=symc
!
!4.1.*def possitions of V, R2
       possv=mapdv(iv,1)
       ir2=mapir(symb,syma,1)
       possr2=mapdr(ir2,1)
!
!4.1.*def dimensions
       nhelp1=nvb(symb)
       nhelp2=nvb(syma)
       nhelp3=nhelp2*(nhelp2+1)/2
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
       nhelp8=nvb(symc)-nhelp6
!
!4.1.*realize definition of V
       call defvhlp9 (wrk(possr2),wrk(possv),nhelp1,nhelp2,             &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
!
       else
!     case symb.ne.symc
!
!4.1.*def possitions of V, R2
       possv=mapdv(iv,1)
       if (syma.ge.symc) then
       ir2=mapir(symb,syma,1)
       else
       ir2=mapir(symb,symc,1)
       end if
       possr2=mapdr(ir2,1)
!
!4.1.*def dimensions
       nhelp1=nvb(symb)
       nhelp2=nvb(syma)
       nhelp3=nvb(symc)
       nhelp4=dimm(mapdv(0,1),syma)
       nhelp5=dimm(mapdv(0,2),symb)
       nhelp6=dimm(mapdv(0,3),symc)
       nhelp7=nvb(syma)-nhelp4
       nhelp8=nvb(symc)-nhelp6
!
!4.1.*realize definition of V
       if (syma.gt.symc) then
       call defvhlp81 (wrk(possr2),wrk(possv),nhelp1,nhelp2,            &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
       else
       call defvhlp82 (wrk(possr2),wrk(possv),nhelp1,nhelp2,            &
     & nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
       end if
!
       end if
!
 4      continue
!
       else
!     not implemented deftyp
!
       end if
!
       return
       end
