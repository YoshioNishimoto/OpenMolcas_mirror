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
        subroutine InsReaW4 (aSGrp,bSGrp,cSGrp,dSGrp,length)
!
!        this routine do:
!        Check which W4 file corresponds to given combination
!        of indexes,
!        - increase the parameter, checking the overal number of
!        required integrals if these block is conted first time
!        - and set corresponding InqW4 to T
!
!
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "parcc.fh"
!
        integer aSGrp,bSGrp,cSGrp,dSGrp,length
!
!        help variables
        integer abSGrp,cdSGrp
        integer abLen,cdLen,abcdLen
        integer dima,dimb,dimc,dimd
        integer pSGrp,qSGrp,rSGrp,sSGrp,dimp,dimq,dimr,dims
        integer pqSGrp,rsSGrp
        integer i
!
!
!        -------- part I - define basic parameters
!
!1        Def dima,dimb,dimc,dimc
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpbe(bSGrp)
        dimc=DimSGrpa(cSGrp)
        dimd=DimSGrpbe(dSGrp)
!
!
!        In steps 2.1 - 2.3 also dimp-dims and pSGrp-sSGrp
!2.1        Def abSGrp
        if (aSGrp.ge.bSGrp) then
          abSGrp=(aSGrp*(aSGrp-1)/2)+bSGrp
          pSGrp=aSGrp
          qSGrp=bSGrp
          dimp=dima
          dimq=dimb
        else
          abSGrp=(bSGrp*(bSGrp-1)/2)+aSGrp
          pSGrp=bSGrp
          qSGrp=aSGrp
          dimp=dimb
          dimq=dima
        end if
!
!2.2        Def cdSGrp
        if (cSGrp.ge.dSGrp) then
          cdSGrp=(cSGrp*(cSGrp-1)/2)+dSGrp
          rSGrp=cSGrp
          sSGrp=dSGrp
          dimr=dimc
          dims=dimd
        else
          cdSGrp=(dSGrp*(dSGrp-1)/2)+cSGrp
          rSGrp=dSGrp
          sSGrp=cSGrp
          dimr=dimd
          dims=dimc
        end if
!
        if (abSGrp.lt.cdSGrp) then
          i=pSGrp
          pSGrp=rSGrp
          rSGrp=i
          i=qSGrp
          qSGrp=sSGrp
          sSGrp=i
          i=dimp
          dimp=dimr
          dimr=i
          i=dimq
          dimq=dims
          dims=i
        end if
!
!
!3.1        Def abLen
        if (aSGrp.eq.bSGrp) then
          abLen=dima*(dima+1)/2
        else
          abLen=dima*dimb
        end if
!
!3.2        Def cdLen
        if (cSGrp.eq.dSGrp) then
          cdLen=dimc*(dimc+1)/2
        else
          cdLen=dimc*dimd
        end if
!
!3.3        Def abcdLen
        abcdLen=abLen*cdLen
!
!
!        -------- part II - read proper integrals (pq|rs) from disc
!
        pqSGrp=(pSGrp*(pSGrp-1)/2)+qSGrp
        rsSGrp=(rSGrp*(rSGrp-1)/2)+sSGrp
!
        if (InqW4(pqSGrp,rsSGrp).eqv..False.) then
          InqW4(pqSGrp,rsSGrp)=.True.
          length=length+abcdLen
        end if
!
!
        return
        end
