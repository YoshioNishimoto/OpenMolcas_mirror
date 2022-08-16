************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
*    Copyright (C) 2020, Bruno Tenorio                                 *
************************************************************************
      SUBROUTINE MKRTDM2(IFSBTAB1,IFSBTAB2,ISSTAB,
     &                  MAPORB,DET1,DET2,
     &                  IF21,IF12,NRT2M,RT2M,
     &                  ISTATE,JSTATE,SPIN)
! The spin coupling matrix elements have the following index-code:
             !SPIN=1 means  K2V (AAB+BBB)
             !SPIN=-1 means SDA (AAA+BBA)
             !SPIN=2 means: bbb
             !SPIN=3 means: aaa
             !SPIN=4 means: aab
             !SPIN=5 means: bba
             !SPIN=6 means: aba
             !SPIN=7 means: bab
! Notice, SPIN here has nothing to do with the spin quantum number. It
! is just a printing code.

      IMPLICIT NONE
      INTEGER SPIN
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER ISSTAB(*),MAPORB(*),NRT2M
      REAL*8 DET1(*),DET2(*)
      REAL*8 RT2M(NRT2M)
      INTEGER NASHT,NASORB,LORBTB,IP,I
      REAL*8 Srt12, Srt13
      REAL*8 GVAL,GAAA,GAAB,GABA,GABB,GBAA,GBAB,GBBA,GBBB
      REAL*8 GVALS, GVALT
!      !INTEGER LSYM1,MSPOJ1,LSYM2,MSPROJ2,ISYOP,MS2OP
!      !INTEGER MPLET1,MPLET2, MPLETD
      INTEGER IAJBLA,IAJBLB,IBJALA,IBJALB,IBJAJB
      INTEGER IBJBJB,IAJAJB
      INTEGER IAJALA,IAJALB,IBJBLA,IBJBLB,IAJBJA,IBJBJA
      INTEGER LORB,JORB,IORB 
      INTEGER JORBA,JORBB,LORBA,LORBB,IORBA,IORBB
      INTEGER ITABS,JTABS,LTABS,JLTABS,IJLTABS
!      INTEGER IBKA,IBKB,IJ,IJIJ,IORBA,IORBB,ITU,ITUVX
!      INTEGER IVABS,IVX,IXABS,JALA,JALB,JBJA,JBLA,JBLB
!      INTEGER JORBA,JORBB,KORB,KORBA,KORBB,LORB,LORBA,LORBB
      INTEGER NASGEM,NSRT2M,ISTATE,JSTATE,SIGNLJ
      LOGICAL IF21,IF12
      !INTEGER job1,job2,ist,jst
#include "symmul.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
      Real*8, Allocatable:: SRT2M(:)
C Given two CI expansions, using a biorthonormal set of SD''s,
C calculate the 2-particle transition density matrix
C in the biorthonormal active orbital basis.
C It will build the contribution from high spin (J->beta,L->beta)
C and low spin (J->beta,L->alpha).
      LORBTB=ISSTAB(3)
C Pick out nr of active orbitals from orbital table:
      NASORB=IWORK(LORBTB+3)
      NASHT=NASORB/2
      NASGEM=(NASORB*(NASORB-1))/2
      !NSRT2M=NASORB*(NASORB*(NASORB-1))/2
      NSRT2M=NASORB**3
      Call mma_allocate(SRT2M,nSRT2M,Label='SRT2M')
      SRT2M(:)=0.0D0
      !ISYOP   = MUL(LSYM1,LSYM2)
      !MS2OP   = MSPROJ1-MSPROJ2
      !MPLETD =  MPLET1 - MPLET2
        CALL SRTDM2(IWORK(LORBTB),ISSTAB,
     &              IFSBTAB1,IFSBTAB2,DET1,DET2,
     &              IF21,IF12,NSRT2M,SRT2M)
      !write(6,*)"BRN SRT2M print",nSRT2M
      !DO I=1,NSRT2M
      !  write(6,*) "SRT2",I,SRT2M(I)
      !END DO   
C Mapping from active spin-orbital to active orbital in external order.
C Note that these differ, not just because of the existence of two
C spin-orbitals for each orbital, but also because the active orbitals
C (external order) are grouped by symmetry and then RAS space, but the
C spin orbitals are grouped by subpartition.
      IP=0
      Srt12=1.0D0/Sqrt(2.0D0)
      Srt13=1.0D0/Sqrt(3.0D0)
      IAJALA=0      ! dummy initialize
      IAJALB=0      ! dummy initialize
      IAJBLA=0      ! dummy initialize
      IAJBLB=0      ! dummy initialize
      IBJALA=0      ! dummy initialize
      IBJALB=0      ! dummy initialize
      IBJBLA=0      ! dummy initialize
      IBJBLB=0      ! dummy initialize
      IAJBJA=0      ! dummy initialize for J=L
      IBJBJA=0      ! dummy initialize for J=L
      IBJAJB=0      ! dummy initialize for J=L
C For high spin density it will keep only beta,beta,beta.
C Notice that beta,beta,beta when J=L is zero.
C For low spin density it'll keep only alpha,beta,alpha.
      DO IORB=1,NASHT
       IORBA=2*IORB-1
       IORBB=2*IORB
       ITABS=MAPORB(IORBA)
       DO JORB=1,NASHT
        JORBA=2*JORB-1
        JORBB=2*JORB
        JTABS=MAPORB(JORBA)
        DO LORB=1,NASHT
         LORBA=2*LORB-1
         LORBB=2*LORB
         LTABS=MAPORB(LORBA)
         JLTABS=LTABS+NASHT*(JTABS-1)
         IF(JORB.GT.LORB) THEN ! When J>L
          SIGNLJ=1.0D0
          IAJALB=IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
          IAJBLA=IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
          IBJBLB=IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
          IBJALB=IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
          IAJALA=IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
          IBJBLA=IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
          !The spin function
          IF(SPIN.EQ.1) THEN ! K^1/2,1/2= (AAB+BBB)
           GAAB=SRT2M(IAJALB)
           GBBB=SRT2M(IBJBLB)
           GVAL=(GAAB+GBBB)
          ELSE IF(SPIN.EQ.-1) THEN ! SDA. K^1/2,-1/2
           GBBA=SRT2M(IBJBLA)
           GAAA=SRT2M(IAJALA)
           GVAL=(GAAA+GBBA)
          ELSE IF(SPIN.EQ.2) THEN 
           GBBB=SRT2M(IBJBLB)
           GVAL=GBBB
          ELSE IF(SPIN.EQ.3) THEN
           GAAA=SRT2M(IAJALA)
           GVAL=GAAA
          ELSE IF(SPIN.EQ.4) THEN
           GAAB=SRT2M(IAJALB)
           GVAL=GAAB
          ELSE IF(SPIN.EQ.5) THEN
           GBBA=SRT2M(IBJBLA)
           GVAL=GBBA
          ELSE IF(SPIN.EQ.6) THEN
           GABA=SRT2M(IAJBLA)
           GVAL=GABA
          ELSE IF(SPIN.EQ.7) THEN
           GBAB=SRT2M(IBJALB)
           GVAL=GBAB
          END IF
         ELSE IF(JORB.EQ.LORB) THEN
         !ELSE ! When J=L
          IAJALB=IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
          IAJBLA=IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
          IBJBLB=IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
          IBJALB=IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
          IAJALA=IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
          IBJBLA=IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
          !The spin function
          IF(SPIN.EQ.1) THEN ! K^1/2,1/2
           GAAB=SRT2M(IAJALB)
           GBBB=SRT2M(IBJBLB)
           GVAL=(GAAB+GBBB)
          ELSE IF(SPIN.EQ.-1) THEN ! SDA. K^1/2,-1/2
           GBBA=SRT2M(IBJBLA)
           GAAA=SRT2M(IAJALA)
           GVAL=(GAAA+GBBA)
          ELSE IF(SPIN.EQ.2) THEN 
           GBBB=SRT2M(IBJBLB)
           GVAL=GBBB
          ELSE IF(SPIN.EQ.3) THEN
           GAAA=SRT2M(IAJALA)
           GVAL=GAAA
          ELSE IF(SPIN.EQ.4) THEN
           GAAB=SRT2M(IAJALB)
           GVAL=GAAB
          ELSE IF(SPIN.EQ.5) THEN
           GBBA=SRT2M(IBJBLA)
           GVAL=GBBA
          ELSE IF(SPIN.EQ.6) THEN
           GABA=SRT2M(IAJBLA)
           GVAL=GABA
          ELSE IF(SPIN.EQ.7) THEN
           GBAB=SRT2M(IBJALB)
           GVAL=GBAB
          END IF
         ELSE IF(JORB.LT.LORB) THEN ! When J<L
          IAJALB=IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
          IAJBLA=IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
          IBJBLB=IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
          IBJALB=IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
          IAJALA=IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
          IBJBLA=IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
          !The spin function
          IF(SPIN.EQ.1) THEN ! K^1/2,1/2
           GAAB=SRT2M(IAJALB)
           GBBB=SRT2M(IBJBLB)
           GVAL=(GAAB+GBBB)
          ELSE IF(SPIN.EQ.-1) THEN ! SDA. K^1/2,-1/2
           GBBA=SRT2M(IBJBLA)
           GAAA=SRT2M(IAJALA)
           GVAL=(GAAA+GBBA)
          ELSE IF(SPIN.EQ.2) THEN 
           GBBB=SRT2M(IBJBLB)
           GVAL=GBBB
          ELSE IF(SPIN.EQ.3) THEN
           GAAA=SRT2M(IAJALA)
           GVAL=GAAA
          ELSE IF(SPIN.EQ.4) THEN
           GAAB=SRT2M(IAJALB)
           GVAL=GAAB
          ELSE IF(SPIN.EQ.5) THEN
           GBBA=SRT2M(IBJBLA)
           GVAL=GBBA
          ELSE IF(SPIN.EQ.6) THEN
           GABA=SRT2M(IAJBLA)
           GVAL=GABA
          ELSE IF(SPIN.EQ.7) THEN
           GBAB=SRT2M(IBJALB)
           GVAL=GBAB
          END IF         
         END IF
         IJLTABS=ITABS+NASHT*(JLTABS-1)
         RT2M(IJLTABS)=GVAL
        END DO
       END DO
      END DO
      CALL mma_deallocate(SRT2M)
      RETURN
      END
