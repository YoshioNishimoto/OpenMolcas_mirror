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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine dumps characters in xml format.                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

subroutine xml_sDump(TagName,opt)

use Definitions, only: iwp

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: TagName
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer(kind=iwp) :: opt
interface
  subroutine xml_cDumpc(name_,nx_name) bind(C,name='xml_cdumpc_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    character(kind=c_char) :: name_(*)
    integer(kind=MOLCAS_C_INT) :: nx_name
  end subroutine xml_cDumpc
  subroutine xml_cDumps(name_,nx_name) bind(C,name='xml_cdumps_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    character(kind=c_char) :: name_(*)
    integer(kind=MOLCAS_C_INT) :: nx_name
  end subroutine xml_cDumps
end interface
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if (opt == 0) then
  call xml_cDumps(TagName,len(TagName))
else
  call xml_cDumpc(TagName,len(TagName))
end if
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine xml_sDump
