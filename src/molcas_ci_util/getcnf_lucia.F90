!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2002, Jeppe Olsen                                      *
!***********************************************************************

subroutine GETCNF_LUCIA(KCNF,KTYP,K,ICONF,IREFSM,NEL)
! Obtain configuration number K .
! Occupation in KCNF in form of old RASSCF ( doubly occ orbs first)
! Type in KTYP
!
! Adapted for LUCIA Jeppe Olsen, summer of 02

implicit real*8(A-H,O-Z)
#include "spinfo.fh"
#include "ciinfo.fh"
dimension KCNF(*), ICONF(*)
! Configuration list is assumed to be in the form used
! in LUCIA, i.e. doubly occupied orbitals are flagged by
! a minus

ICNFB1 = 1
ICNFB2 = 1
KTYP = 0
do JTYP=1,NTYP
  JOP = JTYP-1+MINOP
  JCL = (NEL-JOP)/2
  JOCC = JOP+JCL

  NJCNF = NCNFTP(JTYP,IREFSM)
  if ((K >= ICNFB1) .and. (K <= ICNFB1+NJCNF-1)) then
    KREL = K-ICNFB1+1
    KADD = (KREL-1)*JOCC
    KTYP = JTYP
    ! Outdated ...
    !call ICOPY(JOCC,ICONF(ICNFB2+KADD),1,KCNF,1)

    ! Obtain configuration in standard RASSCF form
    IIBOP = 1
    IIBCL = 1
    do KOCC=1,JOCC
      KORB = ICONF(ICNFB2+KADD-1+KOCC)
      if (KORB < 0) then
        ! Doubly occupied orbital
        KCNF(IIBCL) = abs(KORB)
        IIBCL = IIBCL+1
      else
        ! Singly occupied orbital
        KCNF(JCL+IIBOP) = KORB
        IIBOP = IIBOP+1
      end if
    end do
  end if

  ICNFB1 = ICNFB1+NJCNF
  ICNFB2 = ICNFB2+NJCNF*JOCC
end do

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' Output from GETCNF '
  write(6,*) ' ================== '
  write(6,*) ' Input configuration number : ',K
  write(6,*) ' Corresponding type : ',KTYP
  write(6,*) ' Occupation : '
  NOCC = (NEL+KTYP-1+MINOP)/2
  call IWRTMA(KCNF,1,NOCC,1,NOCC)
end if

return

end subroutine GETCNF_LUCIA