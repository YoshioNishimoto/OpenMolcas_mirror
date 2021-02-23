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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
*  Cho_X_getVtra
*
*> @brief
*>    This routine performs a half-MO-transformation of Cholesky vectors stored in reduced storage
*> @author F. Aquilante
*>
*> @details
*> This routine performs a half-MO-transformation of Cholesky vectors stored in reduced
*> storage. For \p DoRead = ``.true.`` the vectors are read from
*> disk using array \p RedVec as scratch space, whereas for
*> \p DoRead = ``.false.`` the reduced vectors must be supplied in
*> array \p RedVec.
*>
*> Given a set of pointers (\p ipChoT) to target
*> arrays, the routine performs a half-MO-transformation of \p NUMV Cholesky
*> vectors of compound symmetry \p ISYM starting with
*> vector \p IVEC1 and returns them in the target arrays.
*> Each pointer should thereby point to a
*> location where the corresponding Cholesky
*> vector of a given unique symmetry pair
*> of indices has to be stored
*>
*> - \p iSwap = ``0``: \f$ L(k,b,J) \f$ is returned
*> - \p iSwap = ``1``: \f$ L(a,k,J) \f$ is returned
*> - \p iSwap = ``2``: \f$ L(k,J,b) \f$ is returned
*> - \p iSwap = ``3``: \f$ L(a,J,k) \f$ is returned
*>
*> - \p iSkip(syma) = ``0``: skip the symmetry block \f$ a \f$. Any vector \f$ L_{ab} \f$ or \f$ L_{ba} \f$
*>                           with \c syma &times; \c symb = \p ISYM won't be returned in the target array
*>
*> - \p IREDC: reduced set in core at the moment of the call to the routine.
*>             Can be set to ``-1`` (= unknown or undefined) by the calling routine.
*>
*> @param[out] irc     return code
*> @param[in]  RedVec  Vectors stored in reduced set(s) [\p DoRead option off] or scratch space for reading reduced vectors [\p DoRead option on]
*> @param[in]  lRedVec size of the \p RedVec
*> @param[in]  IVEC1   first vector to read
*> @param[in]  NUMV    number of vectors to transform starting from \p IVEC1
*> @param[in]  ISYM    compound symmetry of the Cholesky vectors
*> @param[in]  iSwap   type of the full storage for the half transformed Cholesky vectors
*> @param[in]  IREDC   reduced set in core
*> @param[in]  nDen    total number of densities to which MOs refer
*> @param[in]  kDen    first density for which the MO transformation has to be performed
*> @param[in]  ipMOs   matrix (8 &times; \p nDen) of pointers to the MOs coefficients
*> @param[in]  ipChoT  pointers to the half transformed vectors
*> @param[in]  iSkip   skipping parameters for each symmetry block \f$ (ab) \f$ of compound symmetry \p ISYM
*> @param[in]  DoRead  flag for reading the reduced vectors
************************************************************************
      Subroutine Cho_X_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,
     &                         iSwap,IREDC,nDen,kDen,MOs,ipChoT,
     &                         iSkip,DoRead)
      use Data_Structures, only: CMO_Type
      Implicit Real*8 (a-h,o-z)

      Type (CMO_Type) MOs(nDen)

      Dimension RedVec(lRedVec)
      Integer   nDen,kDen
      Integer   ipChoT(8,nDen)
      Integer   iSkip(*)
      Logical   DoRead
      Character(LEN=13), Parameter:: SECNAM = 'Cho_X_GetVtra'

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Integer, Allocatable:: ipVec(:,:)

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
**************************************************

      MXUSD = 0
      MUSED = 0

C zeroing the target arrays
C--------------------------
      Do iSymp=1,nSym
         IF (iSkip(iSymp).ne.0) THEN
            iSymb = muld2h(ISYM,iSymp)
            Do jDen=kDen,nDen
               If (iSwap.eq.0 .or. iSwap.eq.2) Then ! Lpb,J or LpJ,b
                  Call FZero(Work(ipChoT(iSymp,jDen)),
     &                       SIZE(MOs(jDen)%pA(iSymp)%A,1)*
     &                       nBas(iSymb)*
     &                       NUMV)
               ElseIf (iSwap.eq.1 .or. iSwap.eq.3) Then !Laq,J or LaJ,q
                  Call FZero(Work(ipChoT(iSymp,jDen)),
     &                       nBas(iSymp)*
     &                       SIZE(MOs(jDen)%pA(iSymp)%A,1)*
     &                       NUMV)
               EndIf
            End Do
         ENDIF
      End Do

C --- define local pointers to the target arrays
C ----------------------------------------------
      Call mma_allocate(ipVec,8,nDen,Label='ipVec')

      Do jDen=kDen,nDen
         Do i=1,nSym
            ipVec(i,jDen) = ipChoT(i,jDen)
         End Do
      End do
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      IF (DoRead) THEN
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
       JVEC1 = IVEC1
       IVEC2 = JVEC1 + NUMV - 1

       Do While (jVec1.le.iVec2)

        Call CHO_VECRD(RedVec,lRedVec,JVEC1,IVEC2,ISYM,JNUM,IREDC,MUSED)

        MXUSD = MAX(MXUSD,MUSED)

        If (JNUM.le.0 .or. JNUM.gt.(IVEC2-JVEC1+1)) then
           irc=77
           RETURN
        End If

        jVref = JVEC1 - IVEC1 + 1

        Call cho_vTra(irc,RedVec,lRedVec,jVref,JVEC1,JNUM,NUMV,ISYM,
     &            IREDC,iSwap,nDen,kDen,MOs,ipVec,iSkip)

        if (irc.ne.0) then
           return
        endif

        jVec1 = jVec1 + JNUM
*                                                                     *
***********************************************************************
C --- Updating the local pointers to the target arrays (iff iSwap=0,1)
C --------------------------------------------------------------------
        IF(iSwap.eq.0)THEN    ! Lpb,J

         Do jDen=kDen,nDen
          do iSymp=1,nSym
             n1 = SIZE(MOs(jDen)%pA(iSymp)%A,1)
             iSymb=mulD2h(iSymp,ISYM)
             if(iSkip(iSymp).ne.0)then
               ipVec(iSymp,jDen) = ipVec(iSymp,jDen)
     &                           + n1*nBas(iSymb)*JNUM
             endif
          end do
         End Do

        ELSEIF(iSwap.eq.1)THEN    ! Laq,J

         Do jDen=kDen,nDen
          do iSyma=1,nSym
             iSymq=mulD2h(iSyma,ISYM)
             n2 = SIZE(MOs(jDen)%pA(iSymq)%A,1)
             if(iSkip(iSyma).ne.0)then
               ipVec(iSyma,jDen) = ipVec(iSyma,jDen)
     &                           + nBas(iSyma)*n2*JNUM
             endif
          end do
         End Do

        ELSEIF(iSwap.ne.2 .or. iSwap.ne.3)THEN

         write(6,*)SECNAM//': invalid argument. Iswap= ',Iswap
         irc=66
         return

        ENDIF
***********************************************************************
*                                                                     *
       End Do  ! end the while loop
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      ELSE ! only MO transformation
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
       JNUM = NUMV

       Call cho_vTra(irc,RedVec,lRedVec,1,IVEC1,JNUM,NUMV,ISYM,IREDC,
     &              iSwap,nDen,kDen,MOs,ipVec,iSkip)

       if (irc.ne.0) then
          return
       endif
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      END IF
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      Call mma_deallocate(ipVec)

      irc=0

      RETURN
      END
