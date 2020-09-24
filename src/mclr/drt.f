************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE DRT_MCLR(NVERT0,NVERT,IDRT0,IDOWN0,IVER,IDRT,IDOWN)
C
C     PURPOSE: USING THE UNRESTRICTED DRT TABLE GENERATED BY DRT0 AND
C              THE MASKING ARRAY PRODUCED BY RESTR COPY ALL VALID
C              VERTICES FROM THE OLD TO THE NEW DRT TABLE
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDRT0(NVERT0,5),IDOWN0(NVERT0,0:3)
      DIMENSION IVER(NVERT0)
      DIMENSION IDRT(NVERT,5),IDOWN(NVERT,0:3)
C
C
      DO 30 IV=1,NVERT0
        IVNEW=IVER(IV)
        IF(IVNEW.EQ.0) GOTO 30
        DO 10 ITAB=1,5
          IDRT(IVNEW,ITAB)=IDRT0(IV,ITAB)
10      CONTINUE
        DO 20 IC=0,3
          ID=IDOWN0(IV,IC)
          IDNEW=0
          IF(ID.NE.0) IDNEW=IVER(ID)
          IDOWN(IVNEW,IC)=IDNEW
20      CONTINUE
30    CONTINUE
C
C
C     EXIT
C
      RETURN
      END
