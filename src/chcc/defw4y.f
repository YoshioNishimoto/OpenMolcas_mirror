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
        subroutine  DefW4y (aGrp,bGrp,cGrp,dGrp,w4y)
c
c        this routine do:
c        define w4y key, which indicates, if atleast one
c        W4 file need to be calculated on this node for
c        given a',b',c',d'
c
        implicit none
        integer aGrp,bGrp,cGrp,dGrp,w4y
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "parcc.fh"
c
c        help variables
        integer aSGrp,bSGrp,abSGrp,bSGrpUp,cSGrp,cdSGrp,dSGrp,dSGrpUp
c
        w4y=0
c
c       cycle over a">=b" subgroups for a',b'
        do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
        if (aGrp.eq.bGrp) then
          bSGrpUp=aSGrp
        else
          bSGrpUp=GrpaUp(bGrp)
        end if
        do bSGrp=GrpaLow(bGrp),bSGrpUp
        abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
c
c          cycle over  c">=d" for c',d'
          do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
          if (cGrp.eq.dGrp) then
            dSGrpUp=cSGrp
          else
            dSGrpUp=GrpaUp(dGrp)
          end if
          do dSGrp=GrpaLow(dGrp),dSGrpUp
          cdSGrp=cSGrp*(cSGrp-1)/2+dSGrp
            if (InqW4(abSGrp,cdSGrp)) then
              w4y=w4y+1
            end if
          end do
          end do
c
c       end cycle over a">=b" subgroups
        end do
        end do
c
#else
        w4y=1
c Avoid unused argument warnings
        if (.false.) then
          call Unused_integer(aGrp)
          call Unused_integer(bGrp)
          call Unused_integer(cGrp)
          call Unused_integer(dGrp)
        end if
#endif
c
        return
        end
