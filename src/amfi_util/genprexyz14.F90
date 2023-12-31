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

subroutine genprexyz14(icheckz,interxyz)

use AMFI_global, only: Lmax
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: icheckz(0:Lmax,0:Lmax,0:Lmax,0:Lmax), interxyz(16,0:Lmax,0:Lmax,0:Lmax,0:Lmax)
integer(kind=iwp) :: M1, M2, M3, M4
integer(kind=iwp), external :: mcheckz

!bs ####################################################################
!bs   some quick decision for interaction
!bs ####################################################################
do M4=0,Lmax
  do M3=0,Lmax
    do M2=0,Lmax
      do M1=0,Lmax
        icheckz(m1,m2,m3,m4) = mcheckz(m1,m2,m3,m4)
      end do
    end do
  end do
end do
!bs ####################################################################
!bs   there are at most 16 possible combinations of signs (2**4)
!bs ####################################################################
interxyz(:,:,:,:,:) = 0

return

end subroutine genprexyz14
