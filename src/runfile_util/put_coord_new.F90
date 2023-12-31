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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************
!  Put_Coord_New
!
!> @brief
!>   Write the updated/new symmetry unique Cartesian coordinates of the basis set centers on the runfile
!> @author R. Lindh
!>
!> @details
!> The utility will write updated/new symmetry unique Cartesian coordinates of the basis set centers on the runfile.
!>
!> @param[in] Coord  Array of the updated/new symmetry unique Cartesian coordinates of the basis set centers
!> @param[in] nAtoms Number of symmetry unique basis set centers
!***********************************************************************

subroutine Put_Coord_New(Coord,nAtoms)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: Coord(3,nAtoms)
character(len=*), parameter :: Label = 'GeoNew'

call Put_dArray(Label,Coord,3*nAtoms)

return

end subroutine Put_Coord_New
