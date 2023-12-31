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

subroutine Cho_x_Quit(SecNam,Str1,Str2)
! Purpose: stop execution using and print the trace stack.

implicit none
character(len=*), intent(in) :: SecNam, Str1, Str2

call SysAbendMsg(SecNam,Str1,Str2)

end subroutine Cho_x_Quit
