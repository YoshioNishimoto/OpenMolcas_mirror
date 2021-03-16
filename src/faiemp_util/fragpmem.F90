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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine FragPMem(nHer,MemFrag,la,lb,lr)
!***********************************************************************
!  Object: to compute the number of real*8 the kernel routine will     *
!          need for the computation of a matrix element between two    *
!          cartesian Gaussian functions with the total angular momentum*
!          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
!          lr is the order of the operator (this is only used when the *
!          integrals are computed with the Hermite-Gauss quadrature).  *
!                                                                      *
!  Called from: OneEl                                                  *
!                                                                      *
!                R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: nHer
integer(kind=iwp), intent(out) :: MemFrag
integer(kind=iwp), intent(in) :: la, lb, lr
integer(kind=iwp) :: nExpi, nExpj, i, nElem, maxDensSize, iCnttp, jCnttp, iAng, jAng, iShll, jShll, ip, nac, ncb, MemMlt, nH, &
                     nBasisi, nBasisj

! statement function
nElem(i) = (i+1)*(i+2)/2

MemFrag = 0
maxDensSize = 0
! largest possible fragment energy weighted density matrix
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%nFragType > 0) maxDensSize = max(maxDensSize,dbsc(iCnttp)%nFragDens*(dbsc(iCnttp)%nFragDens+1)/2)
end do

do iCnttp=1,nCnttp
  if (.not. dbsc(iCnttp)%Frag) cycle

  do iAng=0,dbsc(iCnttp)%nVal-1
    iShll = dbsc(iCnttp)%iVal+iAng
    nExpi = Shells(iShll)%nExp
    nBasisi = Shells(iShll)%nBasis
    if (nExpi == 0 .or. nBasisi == 0) cycle

    do jCnttp=iCnttp,nCnttp
      ! still figure out how to loop only over the centers belonging to the
      ! same fragment (keep track of mdc?) ! still to be done !!!
      if (.not. dbsc(jCnttp)%Frag) cycle

      do jAng=0,dbsc(jCnttp)%nVal-1
        jShll = dbsc(jCnttp)%iVal+jAng
        nExpj = Shells(jShll)%nExp
        nBasisj = Shells(jShll)%nBasis
        if (nExpj == 0 .or. nBasisj == 0) cycle
        !go to 1976

        ip = 2*maxDensSize
        nac = nElem(la)*nElem(iAng)
        ip = ip+nExpi*nac
        ip = ip+3*nExpi
        ip = ip+nExpi
        ip = ip+nExpi
        ip = ip+nExpi

        call MltMmP(nH,MemMlt,la,iAng,lr)
        nHer = max(nH,nHer)
        MemFrag = max(MemFrag,ip+nExpi*MemMlt)
        ip = ip-6*nExpi

        ncb = nElem(jAng)*nElem(lb)
        ip = ip+nExpj*ncb
        ip = ip+3*nExpj
        ip = ip+nExpj
        ip = ip+nExpj
        ip = ip+nExpj

        call MltMmP(nH,MemMlt,jAng,lb,lr)
        nHer = max(nH,nHer)
        MemFrag = max(MemFrag,ip+nExpj*MemMlt)
        ip = ip-6*nExpj

        ip = ip+max(nac*max(nExpi,nBasisj),ncb*nBasisj)
        MemFrag = max(MemFrag,ip)
      end do
    end do
  end do
end do

return

end subroutine FragPMem