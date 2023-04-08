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
        subroutine UrobI1 (I1,NaGrp,LunAux)
!
!       vyraba fily so simulovanymi I1       vektormi
!       so spravnou strukturou
!
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
!
        integer NaGrp,LunAux
        real*8 I1(1)
!
!       help variables
        integer aGrp,len
        real*8 schem
!
!1      cycle over a,be Groups
!
        do aGrp=1,NaGrp
!
!1.1      def length
          len=no*DimGrpv(aGrp)*no*(no+1)/2
!
!1.2      full I1 with random numbers
          schem=1.0d-2
          call RNFill (len,I1(1),schem)
!
!1.3      open proper file
!         open (unit=LunAux,file=I1Name(aGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,I1Name(aGrp))
!
!1.4      write I1 into proper file
          write (6,*) aGrp,len
          call wri_chcc (LunAux,len,I1(1))
!
        close (LunAux)
!
        end do
!
!
        return
        end
