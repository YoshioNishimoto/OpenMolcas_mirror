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
      program main
#ifdef _FPE_TRAP_
      Use, Intrinsic :: IEEE_Exceptions
#endif
      implicit None
      Character(LEN=20), Parameter:: Module_Name='scf'
      Integer iReturn
#ifdef _FPE_TRAP_
      Call IEEE_Set_Halting_Mode(IEEE_Usual,.True._4)
#endif

      Call Start(Module_Name)
!      Call xml_Open('module',' ',' ','scf')
      Call scf(ireturn)
!      Call xml_Close('module')
      Call Finish(ireturn)

      end Program Main
