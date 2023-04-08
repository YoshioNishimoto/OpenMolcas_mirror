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
        subroutine MakeChckData (wrk,wrksize,LunAux)
!
!        this routine generate checkeroo data
!        Use id possible only when NaGrp=NaSGRp=1
!
!        assumption nc>=nv>=no (inac preverit dimenzovanie)
!        dimensions: V1    - no*nv*nv2, nv2*nv2
!                    V2    - nc*nv2
!                    V3    - nc*nv*no
!
        implicit none
#include "wrk.fh"
#include "chcc1.fh"
#include "chcc_files.fh"
!
        integer LunAux
!
!        help variables
        integer dim1
        character*6 LunName
        integer PossV1,PossV2,PossV3,PossT
!
        call DistMemChck (PossV1,PossV2,PossV3,PossT)
!
!1        make Q0
        LunName=I0name
        dim1=no*(no+1)*no*(no+1)/4
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
        call MkQ0 (wrk(PossV1))
!
!2        make Q1
        LunName=I1name(1)
        dim1=no*(no+1)*no*nv/2
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
        call MkQ1 (wrk(PossV1))
!
!3        Read Q21
        LunName=I2name(1,1)
        dim1=no*no*nv*nv
        call GetX (Q21(1,1,1,1),dim1,LunAux,LunName,1,1)
!
!
!4        make Q22
        LunName=I3name(1,1)
        dim1=no*(no+1)*nv*(nv+1)/4
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
        call MkQ22 (wrk(PossV1))
!
!5        make Q4, L2k
        LunName=L2name(1,1)
        dim1=nc*nv*(nv+1)/2
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
         call MkL2_chcc (wrk(PossV2))
        dim1=nv*(nv+1)*nv*(nv+1)/4
        call mv0zero (dim1,dim1,wrk(PossV1))
        dim1=nv*(nv+1)/2
        call mc0c1at3b (nc,dim1,nc,dim1,dim1,dim1,                      &
     &                 dim1,nc,dim1,                                    &
     &                 wrk(PossV2),wrk(PossV2),wrk(PossV1))
        call MkQ4 (wrk(PossV1))

!
!        make Q3, L1k
        LunName=L1name(1)
        dim1=nc*no*nv
        call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
         call mv0u (dim1,wrk(PossV3),1,L1k(1,1,1),1)
        dim1=no*nv*nv*(nv+1)/2
        call mv0zero (dim1,dim1,wrk(PossV1))
        dim1=nv*(nv+1)/2
        call mc0c1at3b (nc,dim1,nc,no*nv,dim1,no*nv,dim1,nc,no*nv,      &
     &                 wrk(PossV2),wrk(PossV3),wrk(PossV1))
        call MkQ3 (wrk(PossV1))
!
!
!        make L0
         LunName=L0name
         dim1=nc*no*(no+1)/2
         call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
         call MkL0 (wrk(PossV3))
!
!
!        make OEo,OEv
         call MkOE (wrk(PossOE))
!
        call MkT1T2
!
        return
        end
