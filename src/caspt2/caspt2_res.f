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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine CASPT2_Res(VECROT)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "caspt2_grad.fh"
C
C#include "SysDef.fh"
C
      DIMENSION VECROT(*)
C
      !! 1) Calculate the derivative of the CASPT2 energy with respect
      !!    to the amplitude.
      !! 2) In the standard CASPT2, solve the CASPT2 equation. In the
      !!    diagonal CASPT2, compute the lambda directly.
      !!
      !! L_S = U_{TS}*H_{TU}*U_{US}
      !!     + U_{SS}*(H_{SS} + <\Psi_S|H0-E0|\Psi_S>)*U_{SS}
      !!     + <lambda|H|\Psi0> + <lambda|H0-E0+Eshift|\Psi_S>
C
C     write(6,*) "in CASPT2_res"
      IRHS2  = 7
      CALL PSCAVEC(1.0D+00,IRHS,IRHS2)
C
      !! Construct the partial derivative of the target state
      !! The derivative is constructed in IRHS2
      !! The shift parameters are set to zero, because the actual energy
      !! is computed without them. The reference state has to be
      !! multiplied by two, from the above equation for L_S.
      !! For MS-CASPT2, the rotation is mutiplied later.
      SAV=SHIFT
      SAVI=SHIFTI
      SHIFT=0.0d0
      SHIFTI=0.0d0
      CALL SIGMA_CASPT2(2.0D+00,2.0D+00,IVECX,IRHS2)
      SHIFT=SAV
      SHIFTI=SAVI
C
      !! Add the partial derivative contribution for MS-CASPT2
      !! (off-diagonal elements). The derivative is taken with IVECW
      !! and put in IVECC.
C     write (*,*) "Ifmscoup = ", ifmscoup, nstlag
      IF (IFMSCOUP) Then
        Call RHS_ZERO(IVECC)
        Call PSCAVEC(VECROT(jStLag),IRHS2,IRHS2)
        Do iStLag = 1, nStLag
          Scal = VECROT(iStLag)
          If (iStLag.eq.jStLag) Scal = 0.0d+00
          If (ABS(VECROT(iStLag)).le.1.0d-12) Cycle
          Call MS_Res(1,iStLag,jStLag,Scal)
        End Do
        !! Transform to SR representatin (IRHS).
        CALL PTRTOSR(0,IVECC,IRHS)
        !! Add to IRHS2
        Call PLCVEC(1.0D+00,1.0D+00,IRHS,IRHS2)
      End If
C
      !! Finally, solve the lambda equation.
      !! The following is just a copy-and-paste of eqctl2.f and pcg.f,
      !! but some unnecessary lines (comuptation of energy etc.) are
      !! omitted.
      !! The lambda equation is solved with the shift parameters,
      !! as is the case for the T-amplitude.
C
      iVecXbk = iVecX
      iVecRbk = iVecR
      iRHSbk  = iRHS
      iVecX   = iVecR
      iRHS    = 7
      iVecR   = 8
C
      Call PCG_RES(ICONV)
      IF (ICONV .NE. 0) THEN
        WRITE (6,'(" Lambda equation did not converge...")')
        WRITE (6,'(" Continue anyway?")')
      END IF
C
      iVecX   = iVecXbk
      iVecR   = iVecRbk
      iRHS    = iRHSbk
C
      !! For implicit derivative of S
      IF (IFMSCOUP) THEN
        CALL PTRTOSR(1,IVECW,IRHS)
        Call RHS_ZERO(IVECC)
        Do iStLag = 1, nStLag
          Scal = VECROT(iStLag)*0.5d+00
          If (iStLag.eq.jStLag) Scal = Scal*2.0d+00
          If (ABS(VECROT(iStLag)).le.1.0d-12) Cycle
          Call MS_Res(1,iStLag,jStLag,Scal)
        End Do
        CALL PTRTOSR(0,IVECC,IRHS2)
      END IF
C
      !! Restore contravariant and covariant representations of the
      !! non-variational T-amplitude
      CALL PTRTOC(0,IVECX,IVECC)
      CALL PTRTOC(1,IVECX,IVECC2)
C     Do iCase = 1, 13
C       write(6,*) "icase=",icase
C       Do iSym = 1, nSym
C         nIN = nINDEP(iSym,iCase)
C         IF(NIN.EQ.0) Cycle
C         nAS = nASUP(iSym,iCase)
C         nIS = nISUP(iSym,iCase)

C         Call RHS_ALLO(nIN,nIS,lg_V1)
C         !! Read the solution vector
C         Call RHS_Read(nIN,nIS,lg_V1,iCase,iSym,iVecR)
C         do i = 1, nin*nis
C           write(6,'(i3,f20.10)') i,work(lg_v1+i-1)
C         end do
C         Call RHS_Free(nIN,nIS,lg_V1)
C       End Do
C     End Do
C
C
C
      RETURN
C
      END
C
C-----------------------------------------------------------------------
C
      !! RHS_SGMDIA
      SUBROUTINE CASPT2_ResD(Mode,NIN,NIS,lg_W,DIN,DIS)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION DIN(*),DIS(*)

C Apply the resolvent of the diagonal part of H0 to an RHS array

#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
C-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,mW,LDW)
          CALL CASPT2_ResD2(MODE,NROW,NCOL,DBL_MB(mW),LDW,DIN(iLo),
     &                DIS(jLo),SHIFT,SHIFTI)
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
        CALL GA_Sync()
C       CALL GAdSUM_SCAL(DOVL)
      ELSE
        CALL CASPT2_ResD2(MODE,NIN,NIS,WORK(lg_W),NIN,DIN,DIS,
     &                   SHIFT,SHIFTI)
      END IF
#else
      CALL CASPT2_ResD2(MODE,NIN,NIS,WORK(lg_W),NIN,DIN,DIS,
     &                 SHIFT,SHIFTI)
#endif

      END
C
C-----------------------------------------------------------------------
C
      !! RESDIA
      SUBROUTINE CASPT2_ResD2(Mode,NROW,NCOL,W,LDW,DIN,DIS,
     &                  SHIFT,SHIFTI)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION W(LDW,*),DIN(*),DIS(*)

      DO J=1,NCOL
        DO I=1,NROW
          SCAL = 0.0D+00
          If (Mode.eq.1) Then
            DELTA  = SHIFT+DIN(I)+DIS(J)
            DELINV = DELTA/(DELTA**2+SHIFTI**2)
            !! The following SCAL is the actual residual
            SCAL   = 1.0D+00 - (DIN(I)+DIS(J))*DELINV
C           write(6,*) "residue = ", scal
C           if (abs(residue).ge.1.0d-08) write(6,*) "residue = ", scal
            !! Another scaling is required for lambda
            SCAL   =-SCAL*DELINV
          ELse If (Mode.eq.2) Then
            SCAL   =-SHIFTI/(DIN(I)+DIS(J))
          End If
          W(I,J) = SCAL*W(I,J)
        END DO
      END DO
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PCG_RES(ICONV)
      USE INPUTDATA
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

      INTEGER ICONV

      INTEGER I,ITER
      INTEGER IVECP,IVECT,IVECU
      INTEGER Cho_X_GetTol
      EXTERNAL Cho_X_GetTol
      REAL*8 ALPHA,BETA,PR,PT,UR
      REAL*8 ECORR(0:8,0:MXCASE)
      REAL*8 EAIVX,EATVX,EBJAI,EBJAT,EBVAT,EVJAI,EVJTI,EVJTU
      REAL*8 E2NONV
      REAL*8 OVLAPS(0:8,0:MXCASE)
      REAL*8 DSCALE

C Flag to tell wether convergence was obtained
      ICONV = 0

C Lists of coupling coefficients, used for sigma vector
C generation from non-diagonal blocks of H0.
C     CALL GETMEM('LISTS','ALLO','INTE',LLISTS,NLSTOT)
C     CALL MKLIST(iWORK(LLISTS))


C Mnemonic names for vectors stored on LUSOLV, see EQCTL.
C Here, we use the local names IVECP, IVECT, IVECU which are thus
C to be seen as overlayed areas. The true vectors IVECC and IVECC2
C are computed on return from this routine, so for a while we use them
C for temporaries.
      IVECP=IVECC
      IVECT=IVECC2
      IVECU=IVECC2


      ITER=0
      RNORM=0.0d0

C Solve equations for the diagonal case, in eigenbasis:
C Current solution vector X, Current residual vector R
      CALL PSCAVEC(-1.0D00,IRHS,IVECR)
      CALL PRESDIA(IVECR,IVECX,OVLAPS)
      IF(MAXIT.EQ.0) THEN
       IF(IPRGLB.GE.TERSE) THEN
        WRITE(6,*)
        WRITE(6,'(23A5)')('-----',i=1,23)
        WRITE(6,*)' DIAGONAL CASPT2 APPROXIMATION:'
        GOTO 900
       END IF
      END IF

C Pre-conditioned conjugate gradient:
C R <- R - (H0-E0)*X
      CALL SIGMA_CASPT2(-1.0D00,1.0D00,IVECX,IVECR)
      CALL POVLVEC(IVECR,IVECR,OVLAPS)
      RNORM=SQRT(OVLAPS(0,0))
      IF(RNORM.LT.THRCONV) GOTO 900
      IF(IPRGLB.GE.USUAL) THEN
       WRITE(6,*)
       WRITE(6,*) "CASPT2/RASPT2 with level-shift and ",
     *            "(X)MS-CASPT2/RASPT2 are non-variational,"
       WRITE(6,*) "so the Lambda equation has to be solved"//
     *            " for analytic gradients"
       WRITE(6,*) "Following values are nonsense (or I just don't"//
     *            " know the meaning)"
       WRITE(6,*)
C      WRITE(6,*)'The contributions to the second order'//
C    &     ' correlation energy in atomic units.'
       WRITE(6,'(25A5)')('-----',I=1,25)
       WRITE(6,'(2X,A,A)')
     & 'IT.      VJTU        VJTI        ATVX        AIVX        VJAI ',
     & '       BVAT        BJAT        BJAI        TOTAL       RNORM  '
       WRITE(6,'(25A5)')('-----',I=1,25)
      END IF
      CALL PRESDIA(IVECR,IVECP,OVLAPS)
C PCG iteration loops:
C---------------------
 100  CONTINUE
      CALL POVLVEC(IVECP,IVECP,OVLAPS)
      DSCALE=1.0D00/SQRT(OVLAPS(0,0))
      CALL PSCAVEC(DSCALE,IVECP,IVECP)
      CALL POVLVEC(IVECP,IVECR,OVLAPS)
      PR=OVLAPS(0,0)
      CALL SIGMA_CASPT2(1.0D00,0.0D0,IVECP,IVECT)
      CALL POVLVEC(IVECP,IVECT,OVLAPS)
      PT=OVLAPS(0,0)
      ALPHA=PR/PT
      CALL PLCVEC(ALPHA,1.0D00,IVECP,IVECX)
      CALL PLCVEC(-ALPHA,1.0D00,IVECT,IVECR)
      CALL POVLVEC(IVECR,IVECR,OVLAPS)
      RNORM=SQRT(OVLAPS(0,0))
      IF(RNORM.LT.THRCONV) GOTO 900
      ITER=ITER+1
      CALL POVLVEC(IRHS,IVECX,ECORR)
      EVJTU=ECORR(0,1)
      EVJTI=ECORR(0,2)+ECORR(0,3)
      EATVX=ECORR(0,4)
      EAIVX=ECORR(0,5)
      EVJAI=ECORR(0,6)+ECORR(0,7)
      EBVAT=ECORR(0,8)+ECORR(0,9)
      EBJAT=ECORR(0,10)+ECORR(0,11)
      EBJAI=ECORR(0,12)+ECORR(0,13)
      E2NONV=ECORR(0,0)
      IF(IPRGLB.GE.USUAL) THEN
       WRITE(6,'(1X,I3,1X,10F12.6)') ITER,EVJTU,EVJTI,EATVX,EAIVX,
     &                     EVJAI,EBVAT,EBJAT,EBJAI,E2NONV,RNORM
       CALL XFLUSH(6)
      END IF
      IF(ITER.GE.MAXIT) GOTO 800
      CALL PRESDIA(IVECR,IVECU,OVLAPS)
      UR=OVLAPS(0,0)
      BETA=PR/UR
      CALL PLCVEC(BETA,1.0D00,IVECU,IVECP)
      GOTO 100
C---------------------

 800  CONTINUE
      IF(IPRGLB.GE.TERSE) THEN
       WRITE(6,*)
       WRITE(6,*)' NOT CONVERGED AFTER MAX ITERATIONS.'
      END IF
      ICONV = 16
 900  CONTINUE
      IF(IPRGLB.GE.TERSE) THEN
       WRITE(6,'(25A5)')('-----',I=1,25)
       WRITE(6,*)
C      WRITE(6,*)' FINAL CASPT2 RESULT:'
      END IF
C     CALL POVLVEC(IRHS,IVECX,ECORR)
C     EVJTU=ECORR(0,1)
C     EVJTI=ECORR(0,2)+ECORR(0,3)
C     EATVX=ECORR(0,4)
C     EAIVX=ECORR(0,5)
C     EVJAI=ECORR(0,6)+ECORR(0,7)
C     EBVAT=ECORR(0,8)+ECORR(0,9)
C     EBJAT=ECORR(0,10)+ECORR(0,11)
C     EBJAI=ECORR(0,12)+ECORR(0,13)
C     E2NONV=ECORR(0,0)
C     CALL POVLVEC(IVECX,IVECX,OVLAPS)
C     DENORM=1.0D0+OVLAPS(0,0)
C     REFWGT=1.0D00/DENORM
CPAM Insert: Compute the variational second-order energy.
CPAM Use unshifted H0. Save any shifts, then restore them.
C     SAV=SHIFT
C     SAVI=SHIFTI
C     SHIFT=0.0d0
C     SHIFTI=0.0d0
C     CALL SIGMA_CASPT2(1.0d0,0.0d0,IVECX,IVECT)
C     SHIFT=SAV
C     SHIFTI=SAVI
C     CALL POVLVEC(IVECX,IVECT,OVLAPS)
C     E2CORR=2.0D0*E2NONV+OVLAPS(0,0)
CPAM End of insert.
C     ESHIFT=E2CORR-E2NONV
C     E2TOT=EREF+E2CORR

C     IF(IPRGLB.GT.USUAL.or.iprglb.ne.silent) THEN
C       WRITE(6,*)
C       WRITE(6,*)' Correlation energy /Case, /Symm, and sums:'
C       DO IC=1,13
C        WRITE(6,'(1X,A8,9F12.8)')
C    &      CASES(IC),(ECORR(IS,IC),IS=1,NSYM),ECORR(0,IC)
C       END DO
C       WRITE(6,'(1X,A8,9F12.8)')
C    &    'Summed: ', (ECORR(IS,0),IS=1,NSYM),ECORR(0,0)
C     ENDIF

      IF (IPRGLB.GE.TERSE) THEN
      !  WRITE(6,*)

      !  If (.not.Input % LovCASPT2) Then
C     !     WRITE(6,'(6x,a,f18.10)')'Reference energy:     ',EREF
      !     WRITE(6,'(6x,a,f30.20)')'Reference energy:     ',EREF
C     !     WRITE(6,'(6x,a,f18.10)')'E2 (Non-variational): ',E2NONV
      !     WRITE(6,'(6x,a,f30.20)')'E2 (Non-variational): ',E2NONV
      !     IF(SHIFT.NE.0.0d0.or.SHIFTI.ne.0.0d0) THEN
C     !       WRITE(6,'(6x,a,f18.10)')'Shift correction:     ',ESHIFT
      !       WRITE(6,'(6x,a,f30.20)')'Shift correction:     ',ESHIFT
      !     END IF
C     !     WRITE(6,'(6x,a,f18.10)')'E2 (Variational):     ',E2CORR
      !     WRITE(6,'(6x,a,f30.20)')'E2 (Variational):     ',E2CORR
      !     If (.not.Input % FnoCASPT2) Then
C     !        WRITE(6,'(6x,a,f18.10)')'Total energy:         ',E2TOT
      !        write(6,'(6x,a,f30.20)')'Total energy:         ',E2TOT
      !     Else
      !        WRITE(6,'(6x,a,f18.10,a)')'FNO correction:       ',EMP2,
     &!             '   (estimate)   '
      !        WRITE(6,'(6x,a,f13.5)')
      !        E2TOT=E2TOT+EMP2
      !        WRITE(6,'(6x,a,f18.10,a)')'Total energy:         ',E2TOT,
     &!             '   (FNO-CASPT2) '
      !     EndIf
      !     WRITE(6,'(6x,a,f18.10)')'Residual norm:        ',RNORM
C     !     WRITE(6,'(6x,a,f13.5)') 'Reference weight:     ',REFWGT
      !     WRITE(6,'(6x,a,f30.20)') 'Reference weight:     ',REFWGT
      !  Else
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Reference energy:                 ',EREF
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Active-Site E2 (Non-variational): ',E2NONV
      !     IF(SHIFT.NE.0.0d0.or.SHIFTI.ne.0.0d0) THEN
      !       WRITE(6,'(6x,a,f18.10)')
     &!             'Shift correction:                 ',ESHIFT
      !     END IF
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Active-Site E2 (Variational):     ',E2CORR
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Frozen region E2 :                ',EMP2
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Residual norm:                    ',RNORM
      !     WRITE(6,'(6x,a,f13.5)')
     &!             'Reference weight:                 ',REFWGT
      !     WRITE(6,'(6x,a,f13.5)')
      !     E2TOT=E2TOT+EMP2
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Total energy (LovCASPT2):         ',E2TOT
      !  EndIf
      END IF

C
      RETURN
C
      END SUBROUTINE PCG_RES
