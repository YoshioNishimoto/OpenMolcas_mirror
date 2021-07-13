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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine page_out(KeyWord,nConf,Vector,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Save any vector for further use by the Davidson diagonalization  *
!     Labels identifying the vectors are kept in a stack and to        *
!     minimize a write through cache strategy is applied               *
!                                                                      *
!     calling arguments:                                               *
!     KeyWord : character*16                                           *
!               record identifier                                      *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     Vector  : array of real*8                                        *
!               any vector of length nConf                             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
real*8 Vector(nConf)
character*16 KeyWord
#include "rasdim.fh"
#include "davctl.fh"
#include "WrkSpc.fh"

! check input arguments
if (nConf < 0) then
  write(6,*) 'page_out: nConf less than 0'
  write(6,*) 'nConf = ',nConf
  call Abend()
end if

! search for a matching record identifier
nStk = 0
do iStk=1,(mxMemStk+mxDiskStk)
  if (LblStk(iStk) == KeyWord) nStk = iStk
end do

! there is a matching record identifier:
! overwrite the current record
if (nStk /= 0) then
  if (nStk <= mxMemStk) then
    iMem = memory_address(nStk)
    call dCopy_(nConf,Vector,1,Work(iMem),1)
  else
    iDisk = disk_address(nStk-mxMemStk)
    call DDaFile(LuDavid,1,Vector,nConf,iDisk)
  end if
end if

! there is no matching record identifier:
! create a new record
if (nStk == 0) then
  if (save_mode == mixed_mode_1) then
    if (KeyWord(1:6) == 'CI_vec') then
      if (save_in_memory) then
        nMemStk = nMemStk+1
        iMem = memory_address(nMemStk)
        call dCopy_(nConf,Vector,1,Work(iMem),1)
        LblStk(nMemStk) = KeyWord
        if (nMemStk == mxMemStk) save_in_memory = .false.
      else
        nMemStk = nMemStk+1
        if (nMemStk > mxMemStk) nMemStk = 1
        iMem = memory_address(nMemStk)
        nDiskStk = nDiskStk+1
        if (nDiskStk > mxDiskStk) nDiskStk = 1
        iDisk = disk_address(nDiskStk)
        call DDaFile(LuDavid,1,Work(iMem),nConf,iDisk)
        call dCopy_(nConf,Vector,1,Work(iMem),1)
        LblStk(mxMemStk+nDiskStk) = LblStk(nMemStk)
        LblStk(nMemStk) = KeyWord
      end if
    else
      nDiskStk = nDiskStk+1
      if (nDiskStk > mxDiskStk) nDiskStk = 1
      iDisk = disk_address(nDiskStk)
      call DDaFile(LuDavid,1,Vector,nConf,iDisk)
      LblStk(mxMemStk+nDiskStk) = KeyWord
    end if
  end if
  if (save_mode == mixed_mode_2) then
    if (save_in_memory) then
      nMemStk = nMemStk+1
      iMem = memory_address(nMemStk)
      call dCopy_(nConf,Vector,1,Work(iMem),1)
      LblStk(nMemStk) = KeyWord
      if (nMemStk == mxMemStk) save_in_memory = .false.
    else
      nMemStk = nMemStk+1
      if (nMemStk > mxMemStk) nMemStk = 1
      iMem = memory_address(nMemStk)
      nDiskStk = nDiskStk+1
      if (nDiskStk > mxDiskStk) nDiskStk = 1
      iDisk = disk_address(nDiskStk)
      call DDaFile(LuDavid,1,Work(iMem),nConf,iDisk)
      call dCopy_(nConf,Vector,1,Work(iMem),1)
      LblStk(mxMemStk+nDiskStk) = LblStk(nMemStk)
      LblStk(nMemStk) = KeyWord
    end if
  end if
end if

return

end subroutine page_out