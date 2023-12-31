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
! Copyright (C) 2023, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine ReadBPT2(LuGAMMA2,B_PT2,nCalAO,kS,lS,MxInShl)
! Read back-transformed density elements of the Kth and Lth shells
! All elements of the Jth shell (auxiliary functions) are read
! Only for C1 symmetry

  use Definitions, only: wp, iwp
  use iSD_data, only: iSD
  use Gateway_global, only: force_part_c
  use pso_stuff, only: nBasA
  use SOAO_Info, only: iAOtSO

  implicit none

  real(kind=wp), intent(out) :: B_PT2(nBasA,MxInShl,MxInShl)
  integer(kind=iwp), intent(in) :: kS, lS, LuGAMMA2, MxInShl
  integer(kind=iwp), intent(inout) :: nCalAO
  integer(kind=iwp) :: i3, i4, kAO, kAOk, kBas, kBsInc, kCmp, kSO, kSOk, kSO0, lAO, lAOl, lBas, lBsInc, lCmp, lSO, lSOl, lSO0

  kCmp = iSD(2,kS)
  kBas = iSD(3,kS)
  kAO  = iSD(7,kS)

  lCmp = iSD(2,lS)
  lBas = iSD(3,lS)
  lAO  = iSD(7,lS)

  if (force_part_c) then
    kBsInc = (kBas+1)/2
    lBsInc = (lBas+1)/2
  else
    kBsInc = kBas
    lBsInc = lBas
  end if

  lSO0 = iAOtSO(lAO+1,0)-1
  kSO0 = iAOtSO(kAO+1,0)-1
  do i4=1,lCmp
    lSO = iAOtSO(lAO+i4,0)
    do i3=1,kCmp
      kSO = iAOtSO(kAO+i3,0)
      do lAOl=0,min(lBsInc,lBas)-1
        lSOl = lSO+lAOl
        do kAOk=0,min(kBsInc,kBas)-1
          kSOk = kSO+kAOk
          nCalAO = nCalAO + 1
!         write (*,'("k,l = ",4i5)') ksok,lsol,ksok-kso0, lsol-lso0
          read(unit=LuGAMMA2,rec=nCalAO) B_PT2(1:nBasA,kSOk-kSO0,lSOl-lSO0)
        end do
      end do
    end do
  end do

  return

end subroutine ReadBPT2
