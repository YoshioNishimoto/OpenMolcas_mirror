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
       subroutine GetIntPoss
!
!        this routine read T3IntPos array from the first record
!        of t3nam file
!
        implicit none
#include "t31.fh"
!
        integer lun
!
       lun=1
       call daname (lun,t3nam)
       daddr(lun)=0
       call idafile (lun,2,T3IntPos,maxorb,daddr(lun))
       call daclos (lun)
!
       return
       end
