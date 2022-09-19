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
        subroutine gen_vvvo(occ_ind,w3,l1_1,l2_1,tmp)
!
! this routine do
!
! regenerate VVVo integrals from cholesky vectors
!
! -------------------
!
! structure of the cholesky vector files :
!
!       L1(m,I ,A') L1vcxx xx - Group of A'
!
!       L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
!                          yy - Group of B'
!
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
!
        integer a,b,c,dima,dimb,dimc,occ_ind
        integer length
        integer lasta,lastb,lastc
!
        real*8 w3(1:(nv*(nv+1))/2,1:nv)
        real*8 tmp(*),l1_1(*),l2_1(*)
!
        integer a_tmp,b_tmp,c_tmp
!
! algoritmus je dobry ak maxdim > no
! inak treba vymenit citanie L1 za L2
!
! dalo by sa to urobit podstatne lepsie, kedby
! dircc nevyzadoval VVV ako (ab,c) ale ako (a,b,c)
!
! mozno urob sort L1i (m,c') <- L1(m,i,c')
! ---
!
!1         loop over a'
!
        do a=1,NvGrp
!
!2        loop over b'
!
        do b=1,a
!
!2.1        read L2(m,a',b')
!
         if (a.eq.b) then  ! a=b
! open the pertinent file
!
!mp@        dima=nv/NvGrp
!mp@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
        dima=DimGrpaR(a)
!
!mp!         write (6,*) 'dima = ',dima
         dimb=dima
         length=(dima*(dima+1)*nc)/2
!mp!         write (6,*) 'length L2Name(a,b) = ',L2Name(a,b),length
!mp!     write (6,*) 'file size (g77) = ',16+length*8
!mp!         write (6,*) 'file size L2Name(a,b) (ifort) = ',
!mp!     & L2Name(a,b),8+length*8
!
        call GetX_t3 (tmp,length,LunAux,L2Name(a,b),1,1)
!
         else ! a>b
! open the pertinent file
!
!mp@@        dima=nv/NvGrp
!mp@@        dimb=nv/NvGrp
!mp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
!mp@@         if (b.eq.NvGrp) dimb=nv-((NvGrp-1)*dimb)
        dima=DimGrpaR(a)
        dimb=DimGrpaR(b)
!
!mp!         write (6,*) 'dima, dimb = ',dima,dimb
         length=dima*dimb*nc
!mp!         write (6,*) 'length L2Name(a,b) = ',L2Name(a,b),length
!mp!     write (6,*) 'file size (g77) = ',16+length*8
!mp!         write (6,*) 'file size L2Name(a,b) (ifort) = ',
!mp!     & L2Name(a,b),8+length*8
!
        call GetX_t3 (tmp,length,LunAux,L2Name(a,b),1,1)
        end if
!
!2.2        map  L2_1(a',b',m) <- tmp(m,a',b')
!
        if (a.eq.b) then ! expand and map
! expand and map l2_1 (a',b',m) <- tmp (m,ab')
        call exMap3_231 (tmp,l2_1,nc,dima)
        else
        call Map3_231_t3 (tmp,l2_1,nc,dima,dimb)
        end if
!
!3         loop over c'
!
        do c=1,NvGrp
!
!3.1        read L1(m,i,c')
!
!mp@@        dimc=nv/NvGrp
!mp@@         if (c.eq.NvGrp) dimc=nv-((NvGrp-1)*dimc)
        dimc=DimGrpaR(c)
!
!mp!         write (6,*) 'dimc = ',dimc
         length=nc*no*dimc
!mp!         write (6,*) 'length L1Name(c) = ',L1Name(c),length
!mp!         write (6,*) 'file size L1Name(c) (ifort) = ',
!mp!     & L1Name(c),8+8*length
!
        call GetX_t3 (tmp,length,LunAux,L1Name(c),1,1)
!
!3.2        extract l1_1 (m,c')_i <- tmp (m,i,c')
! toto by sa dalo nahradit mapovanim
!
        call ext_o_32 (tmp,l1_1,nc,no,dimc,occ_ind)
!
!3.2.1        zero tmp
!
        call zeroma(tmp,1,dima*dimb*dimc)
!
!3.3         mult tmp (a',b',c') <- L2_1 (a',b',m) l1_1 (m,c')
!
               call mc0c1a3b                                            &
     & (dima*dimb,nc,nc,dimc,dima*dimb,dimc,                            &
     & dima*dimb,nc,dimc,l2_1,l1_1,tmp)
!
!3.4        add W(ab,c) <- tmp (a'b'c')
!
!mp@@        lasta=(a-1)*(nv/NvGrp)
!mp@@        lastb=(b-1)*(nv/NvGrp)
!mp@@        lastc=(c-1)*(nv/NvGrp)
!
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
!
          lastb=0
        if (b.gt.1) then
          do b_tmp=1,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
!
          lastc=0
        if (c.gt.1) then
          do c_tmp=1,c-1
          lastc=lastc+DimGrpaR(c_tmp)
          end do
        end if
!
!mp!        write (6,'(A,3(i4),2x,3(i4))') 'lasta, lastb, lastc = ',
!mp!     & lasta,lastb,lastc,a,b,c

! sme v gen_vvvo
        call grow_w3 (w3,tmp,                                           &
     & nv,nv,dima,dimb,dimc,lasta,lastb,lastc)
!
!3.5        end loop over c'
        end do
!4        end loop over b'
        end do
!5        end loop over a'
        end do
!
        return
        end
