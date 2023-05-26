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
      SUBROUTINE CHO_RDBUF(LENGTH,BUF,IBUF,LENBUF,IUNIT)
!
!     Purpose: read buffer from disk.
!
      Implicit Real*8 (a-h,o-z)
      REAL*8 BUF(LENBUF)
      INTEGER   IBUF(4,LENBUF)

      READ(IUNIT) LENGTH,BUF,IBUF

      END
