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

subroutine MKCLIST(ISM,IDOWN,NOW,IOW,ICASE,ISCR)
! PURPOSE: CONSTRUCT THE COMPRESSED CASE-LIST, I.E.,
!          STORE THE STEP VECTOR FOR ALL POSSIBLE WALKS
!          IN THE ARRAY ICASE. GROUPS OF 15 CASES ARE PACKED
!          INTO ONE INTEGER WORD.

implicit real*8(A-H,O-Z)
#include "rasdim.fh"
#include "general.fh"
#include "gugx.fh"
dimension ISM(NLEV), IDOWN(NVERT,0:3)
dimension NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV)
dimension ICASE(NICASE)
dimension ISCR(3,0:NLEV)
parameter(IVERT=1,ISYM=2,ISTEP=3)

! CLEAR ARRAY NOW. IT WILL BE RESTORED FINALLY

do IHALF=1,2
  do MV=1,NMIDV
    do IS=1,NSYM
      NOW(IHALF,IS,MV) = 0
    end do
  end do
end do

! START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.

do IHALF=1,2
  if (IHALF == 1) then
    IVTSTA = 1
    IVTEND = 1
    LEV1 = NLEV
    LEV2 = MIDLEV
  else
    IVTSTA = MIDV1
    IVTEND = MIDV2
    LEV1 = MIDLEV
    LEV2 = 0
  end if

  ! LOOP OVER VERTICES STARTING AT TOP OF SUBGRAPH

  do IVTOP=IVTSTA,IVTEND
    ! SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH:
    LEV = LEV1
    ISCR(IVERT,LEV) = IVTOP
    ISCR(ISYM,LEV) = 1
    ISCR(ISTEP,LEV) = -1
100 continue
    if (LEV > LEV1) goto 400
    ! FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX:
    IVT = ISCR(IVERT,LEV)
    do ISTP=ISCR(ISTEP,LEV)+1,3
      IVB = IDOWN(IVT,ISTP)
      if (IVB == 0) goto 110
      goto 200
110   continue
    end do
    ! ALT A -- NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
    ISCR(ISTEP,LEV) = -1
    LEV = LEV+1
    goto 100
    ! ALT B -- SUCH AN ARC WAS FOUND. WALK DOWN:
200 continue
    ISCR(ISTEP,LEV) = ISTP
    ISML = 1
    if ((ISTP == 1) .or. (ISTP == 2)) ISML = ISM(LEV)
    LEV = LEV-1
    ISCR(ISYM,LEV) = MUL(ISML,ISCR(ISYM,LEV+1))
    ISCR(IVERT,LEV) = IVB
    ISCR(ISTEP,LEV) = -1
    if (LEV > LEV2) goto 100
    ! WE HAVE REACHED THE LOWER LEVEL. THE WALK IS COMPLETE.
    ! MIDVERTEX NUMBER:
    MV = ISCR(IVERT,MIDLEV)+1-MIDV1
    ! SYMMETRY LABEL OF THIS WALK:
    IWSYM = ISCR(ISYM,LEV2)
    ! ITS ORDERING NUMBER WITHIN THE SAME BATCH OF (IHALF,IWSYM,MV):
    ILND = 1+NOW(IHALF,IWSYM,MV)
    NOW(IHALF,IWSYM,MV) = ILND
    ! CONSEQUENTLY, THE POSITION IMMEDIATELY BEFORE THIS COMPRESSED WALK:
    IPOS = IOW(IHALF,IWSYM,MV)+(ILND-1)*NIPWLK
    ! PACK THE STEPS IN GROUPS OF 15 LEVELS PER INTEGER:
    do LL=LEV2+1,LEV1,15
      IC = 0
      do L=min(LL+14,LEV1),LL,-1
        IC = 4*IC+ISCR(ISTEP,L)
      end do
      IPOS = IPOS+1
      ICASE(IPOS) = IC
    end do
    ! FINISHED WITH THIS WALK. BACK UP ONE LEVEL AND TRY AGAIN:
    LEV = LEV+1
    goto 100
400 continue
  end do
end do

! EXIT

return

end subroutine MKCLIST