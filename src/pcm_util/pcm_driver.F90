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

subroutine PCM_Driver(DMat,V,Q,nTs)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTs
real(kind=wp), intent(inout) :: DMat(nTs,nTs)
real(kind=wp), intent(in) :: V(2,nTs)
real(kind=wp), intent(out) :: Q(2,nTs)
integer(kind=iwp) :: iTs, jTs
real(kind=wp) :: tmp

! Computes PCM solvation charges given the nuclear and electronic electrostatic
! potential on each tessera.
! Modifies nuclear repulsion, one-electron and two electron terms.

Q(:,:) = Zero
do iTs=1,nTs
  do jTs=1,nTs
    tmp = Half*(DMat(iTs,jTs)+DMat(jTs,iTs))
    DMat(iTs,jTs) = tmp
    DMat(jTs,iTs) = tmp
  end do
end do

do iTs=1,nTs
  do jTs=1,nTs
    ! Actual nuclear charge
    Q(1,iTs) = Q(1,iTs)+DMat(iTs,jTs)*V(1,jTs)
    ! Effective nuclear charge
    !_rl Q(1,iTs) = Q(1,iTs)+Half*(DMat(iTs,jTs)+DMat(jTs,iTs))*V(1,jTs)
    ! Actual electronic charge
    Q(2,iTs) = Q(2,iTs)+DMat(iTs,jTs)*V(2,jTs)
  end do
end do

!do i=1,nTs   ! yma delete later
!  write(u6,*) ' == V(2,iTs) diff == ',i,V(2,i)
!end do

return

end subroutine PCM_Driver
