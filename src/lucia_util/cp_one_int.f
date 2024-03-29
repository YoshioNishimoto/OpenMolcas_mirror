************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CP_ONE_INT(W1,NDIM)
      use GLBBAS, only: INT1, INT1O
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer nDIM
      Real*8 W1(NDIM)
#include "mxpdim.fh"
#include "orbinp.fh"

      INT1(:)=0.0D0
      INT1(1:NDIM)=W1(1:NDIM)
      INT1O(:)=0.0D0
      INT1O(:)=INT1(:)

      End
