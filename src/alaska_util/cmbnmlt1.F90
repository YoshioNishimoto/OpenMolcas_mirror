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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine CmbnMlt1(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,Alpha,Beta,Grad,nGrad,DAO,IfGrad,IndGrd,iStab,jStab,kOp,nOrdOp,Force)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the multipole operator matrix.       *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChBas
use Index_Functions, only: C_Ind
use Constants, only: Zero, Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nGrad, IndGrd(3,2), iStab, jStab, kOp(2), nOrdOp
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1,0:nOrdOp), Zeta(nZeta), Alpha(nZeta), Beta(nZeta), &
                             DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2), Force(*)
real(kind=wp), intent(inout) :: rKappa(nZeta), Grad(nGrad)
real(kind=wp), intent(inout) :: rFinal(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6)
logical(kind=iwp), intent(in) :: IfGrad(3,2)
integer(kind=iwp) :: i1, i2, iCar, iCn, icomp, iGrad, ipa, ipb, iPrint, iRout, ixa, ixb, ixop, iya, iyaMax, iyb, iybMax, iyop, &
                     iza, izb, iZeta, izop, nDAO
real(kind=wp) :: Fact, ff, ps, xa, xb, ya, yb, za, zb
real(kind=wp), parameter :: exp32 = -Three/Two
integer(kind=iwp), external :: iPrmt
real(kind=wp), external :: DDot_
#include "print.fh"

iRout = 134
iPrint = nPrint(iRout)
do iZeta=1,nZeta
  rKappa(iZeta) = rKappa(iZeta)*Zeta(iZeta)**exp32
end do
if (iPrint >= 99) then
  call RecPrt(' In CmbnMlt1: Zeta  ',' ',Zeta,1,nZeta)
  call RecPrt(' In CmbnMlt1: rKappa',' ',rKappa,1,nZeta)
  call RecPrt(' In CmbnMlt1: Alpha ',' ',Alpha,1,nZeta)
  call RecPrt(' In CmbnMlt1: Beta  ',' ',Beta,1,nZeta)
  call RecPrt(' In CmbnMlt1: DAO  ',' ',Dao,nZeta,(la+1)*(la+2)*(lb+1)*(lb+2)/4)
end if
!.... Loop over cartesian components of operator.......
do ixop=0,nOrdOp
  do iyop=0,nOrdOp-ixop
    izop = nOrdOp-ixop-iyop
    icomp = C_Ind(nOrdOp,ixop,izop)
    ff = Force(icomp)
    if (ff == Zero) cycle
    do ixa=0,la
      iyaMax = la-ixa
      do ixb=0,lb
        iybMax = lb-ixb
        do iya=0,iyaMax
          iza = la-ixa-iya
          ipa = C_Ind(la,ixa,iza)
          do iyb=0,iybMax
            izb = lb-ixb-iyb
            ipb = C_Ind(lb,ixb,izb)

            ! Combine overlap integrals

            if (IfGrad(1,1)) then
              if (ixa > 0) then
                xa = real(-ixa,kind=wp)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)* &
                          (Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ixop)+ &
                                         xa*Rnxyz(iZeta,1,ixa-1,ixb,ixop))* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              else
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)* &
                           Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ixop)* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              end if
            end if
            if (IfGrad(1,2)) then
              if (ixb > 0) then
                xb = real(-ixb,kind=wp)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,4) = rKappa(iZeta)* &
                           (Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ixop)+ &
                                         xb*Rnxyz(iZeta,1,ixa,ixb-1,ixop))* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              else
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,4) = rKappa(iZeta)* &
                            Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ixop)* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              end if
            end if
            if (IfGrad(2,1)) then
              if (iya > 0) then
                ya = real(-iya,kind=wp)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,2) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                          (Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iyop)+ &
                                         ya*Rnxyz(iZeta,2,iya-1,iyb,iyop))* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              else
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,2) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                           Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iyop)* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              end if
            end if
            if (IfGrad(2,2)) then
              if (iyb > 0) then
                yb = real(-iyb,kind=wp)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,5) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                           (Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iyop)+ &
                                         yb*Rnxyz(iZeta,2,iya,iyb-1,iyop))* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              else
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,5) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                            Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iyop)* &
                                            Rnxyz(iZeta,3,iza,izb,izop)
                end do
              end if
            end if
            if (IfGrad(3,1)) then
              if (iza > 0) then
                za = real(-iza,kind=wp)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,3) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                          (Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,izop)+ &
                                         za*Rnxyz(iZeta,3,iza-1,izb,izop))
                end do
              else
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,3) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                           Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,izop)
                end do
              end if
            end if
            if (IfGrad(3,2)) then
              if (izb > 0) then
                zb = real(-izb,kind=wp)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,6) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                           (Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,izop)+ &
                                         zb*Rnxyz(iZeta,3,iza,izb-1,izop))
                end do
              else
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,6) = rKappa(iZeta)* &
                                            Rnxyz(iZeta,1,ixa,ixb,ixop)* &
                                            Rnxyz(iZeta,2,iya,iyb,iyop)* &
                            Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,izop)
                end do
              end if
            end if

          end do
        end do
      end do
    end do

    ! Trace the gradient integrals

    nDAO = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
    if (iPrint >= 99) then
      call RecPrt(' S(1)',' ',rFinal,nDAO,6)
      call RecPrt('   D ',' ',DAO,nDAO,1)
    end if
    do iCn=1,2
      do iCar=1,3
        if (IndGrd(iCar,iCn) /= 0) then
          ! Accumulate contribution to the gradient
          iGrad = abs(IndGrd(iCar,iCn))
          if (iCn == 1) then
            i1 = iCar
            i2 = iCar+3
            ps = real(iPrmt(kOp(1),iChBas(1+iCar)),kind=wp)
            Fact = real(iStab,kind=wp)/real(nIrrep,kind=wp)
          else
            i1 = iCar+3
            i2 = iCar
            ps = real(iPrmt(kOp(2),iChBas(1+iCar)),kind=wp)
            Fact = ps*real(jStab,kind=wp)/real(nIrrep,kind=wp)
          end if
          Fact = Fact*ff
          if (IndGrd(iCar,iCn) < 0) then
            ! Gradient via translational invariance.
            Grad(iGrad) = Grad(iGrad)-Fact*DDot_(nDAO,DAO,1,rFinal(1,1,1,i2),1)
          else
            Grad(iGrad) = Grad(iGrad)+Fact*DDot_(nDAO,DAO,1,rFinal(1,1,1,i1),1)
          end if
        end if
      end do
    end do
  end do
end do

return

end subroutine CmbnMlt1
