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
! Copyright (C) 2019, Oskar Weser                                      *
!***********************************************************************

module fortran_strings

use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
use Definitions, only: wp, iwp, MOLCAS_C_INT

implicit none
private

public :: str, to_lower, to_upper, operator(.in.), split, count_char, StringWrapper_t, Cptr_to_str, char_array

! This type exists to have an array of string pointers
! and to allow unequally sized strings.
! NOTE: Due to old compilers this had to be
! character(len=1), dimension(:), allocatable
! if possible change it to
! character(len=:), allocatable
type :: StringWrapper_t
  character(len=1), allocatable :: str(:)
end type

!>  @brief
!>    Convert to Fortran string
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is a generic procedure that accepts Fortran integer,
!>  Fortran real, and Fortran arrays with single character elements.
!>
!>  @param[in] A Fortran integer or real, or character array.
interface str
  module procedure :: I_to_str, R_to_str, character_array_to_str
end interface

interface
  pure function strlen_c(c_string) bind(C,name='strlen_wrapper')
    import :: c_ptr, MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: strlen_c
    type(c_ptr), intent(in) :: c_string
  end function strlen_c
end interface

interface operator(.in.)
  module procedure :: substr_in_str
end interface

character(len=*), parameter :: UPPERCASE_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', &
                               lowercase_chars = 'abcdefghijklmnopqrstuvwxyz'

contains

pure function I_to_str(i) result(str)
  character(len=:), allocatable :: str
  integer(kind=iwp), intent(in) :: i
  character(len=range(i)+2) :: tmp
  write(tmp,'(I0)') I
  str = trim(tmp)
end function I_to_str

pure function R_to_str(x) result(str)
  character(len=:), allocatable :: str
  real(kind=wp), intent(in) :: x
  character(len=range(x)+2) :: tmp
  write(tmp,'(I0)') x
  str = trim(tmp)
end function R_to_str

pure function character_array_to_str(array) result(res)
  character(len=:), allocatable :: res
  character(len=1), intent(in) :: array(:)
  integer(kind=iwp) :: i, L

  L = size(array)
  allocate(character(len=L) :: res)
  do i=1,L
    res(i:i) = array(i)
  end do
end function character_array_to_str

!> Convert C string pointer to Fortran string.
function Cptr_to_str(c_string) result(res)
  character(len=:), allocatable :: res
  type(c_ptr), intent(in) :: c_string
  character, pointer :: string(:)
  integer(kind=iwp) :: i, L
  L = int(strlen_c(c_string))
  allocate(character(len=L) :: res)
  call c_f_pointer(c_string,string,[L])
  do i=1,L
    res(i:i) = string(i)
  end do
end function Cptr_to_str

!> Changes a string to upper case
pure function to_upper(in_str) result(string)
  character(len=*), intent(in) :: in_str
  character(len=len(in_str)) :: string
  integer(kind=iwp) :: ic, i, L

  L = len_trim(in_str)
  do i=1,L
    ic = index(lowercase_chars,in_str(i:i))
    if (ic > 0) then
      string(i:i) = UPPERCASE_chars(ic:ic)
    else
      string(i:i) = in_str(i:i)
    end if
  end do
  string(L+1:) = ' '
end function to_upper

!> Changes a string to lower case
pure function to_lower(in_str) result(string)
  character(len=*), intent(in) :: in_str
  character(len=len(in_str)) :: string
  integer(kind=iwp) :: ic, i, L

  L = len_trim(in_str)
  do i=1,L
    ic = index(UPPERCASE_chars,in_str(i:i))
    if (ic > 0) then
      string(i:i) = lowercase_chars(ic:ic)
    else
      string(i:i) = in_str(i:i)
    end if
  end do
  string(L+1:) = ' '
end function to_lower

pure function substr_in_str(substring,string)
  logical(kind=iwp) :: substr_in_str
  character(len=*), intent(in) :: string, substring
  substr_in_str = index(string,substring) /= 0
end function substr_in_str

!> @brief
!> Split a string at delimiter.
subroutine split(string,delimiter,res)
  character(len=*), intent(in) :: string
  character(len=1), intent(in) :: delimiter
  type(StringWrapper_t), allocatable, intent(out) :: res(:)
  integer(kind=iwp) :: i, n, low

  allocate(res(count_char(string,delimiter)+1))

  ! NOTE: this function is unnecessarily complicated,
  ! because StringWrapper_t cannot have character(len=:), allocatable
  ! components. (Old compilers)

  low = 1; n = 1
  do i=1,len(string)
    if (string(i:i) == delimiter) then
      allocate(character(len=1) :: res(n)%str(i-low))
      res(n)%str(:) = char_array(string(low:i-1))
      n = n+1
      low = i+1
    end if
  end do

  if (n == size(res)) then
    allocate(character(len=1) :: res(n)%str(len(string(low:))))
    res(n)%str(:) = char_array(string(low:))
  end if
end subroutine split

pure function char_array(string) result(res)
  character(len=*), intent(in) :: string
  character(len=1) :: res(len(string))
  integer(kind=iwp) :: i
  do i=1,len(string)
    res(i) = string(i:i)
  end do
end function char_array

!> @brief
!> Count the occurence of a character in a string.
pure function count_char(str,char) result(c)
  integer(kind=iwp) :: c
  character(len=*), intent(in) :: str
  character(len=1), intent(in) :: char
  integer(kind=iwp) :: i
  c = 0
  do i=1,len(str)
    if (str(i:i) == char) c = c+1
  end do
end function count_char

end module fortran_strings
