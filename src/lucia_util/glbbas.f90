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
Module GLBBAS
Private
#include "mxpdim.fh"
      Integer ::    KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,        &
     &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR,             &
     &              KSBIA,KSBIB,KPNIJ,KIJKK,KFI,KINH1,                  &
     &              KMOAOIN,KMOAOUT,KLOCCLS,NOCCLS_G,KPGINT1(MXPOBS),   &
     &              KINT1O,KFIO,KPGINT1A(MXPOBS),                       &
     &              KMOMO,KSRHO1,KFIZ,                                  &
     &              KINT1_SIMTRH,KINT2_SIMTRH,                          &
     &              KPINT1_SIMTRH,KPINT2_SIMTRH,KINH1_NOCCSYM,          &
     &              KICONF_REO(8),                                      &
     &              KSDREO_I(8),KZ_PTDT(MXPORB+1),                      &
     &              KREO_PTDT(MXPORB+1)
      Real*8, Allocatable:: VEC3(:)

! DETERMINE BASE ADRESSES
!             DFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
!             CFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
!             DTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
!             CONF_OCC(I) : SPACE FOR STORING  NCNSM
!                            CONFIGURATION EXPANSIONS
      Integer, Allocatable:: DFTP(:)
      Integer, Allocatable:: CFTP(:)
      Real*8,  Allocatable:: DTOC(:)
      Type Array
         Integer, Allocatable:: I(:)
      End Type  Array
      Type (Array) CONF_OCC(8)

      Public        KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,        &
     &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR,             &
     &              KSBIA,KSBIB,VEC3,KPNIJ,KIJKK,KFI,KINH1,             &
     &              KMOAOIN,KMOAOUT,KLOCCLS,NOCCLS_G,KPGINT1,           &
     &              KINT1O,KFIO,KPGINT1A,                               &
     &              KMOMO,KSRHO1,KFIZ,                                  &
     &              KINT1_SIMTRH,KINT2_SIMTRH,                          &
     &              KPINT1_SIMTRH,KPINT2_SIMTRH,KINH1_NOCCSYM,          &
     &              CONF_OCC,KICONF_REO,                              &
     &              DFTP,CFTP,DTOC,                                     &
     &              KSDREO_I,KZ_PTDT,                                   &
     &              KREO_PTDT
End Module GLBBAS
