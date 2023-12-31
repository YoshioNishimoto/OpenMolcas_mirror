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

subroutine GenVoronoi(nR_Eff,Alpha,rm,iNQ)
!***********************************************************************
!                                                                      *
!     This version of GenVoronoi computes the radial quadrature points *
!     and computes data useful for the angular quadrature.             *
!     The angular part is generated by Subblock.                       *
!                                                                      *
!***********************************************************************

use NQ_Structure, only: NQ_Data
use nq_Info, only: L_Quad, lMax_NQ, nR, Quadrature
use stdalloc, only: mma_allocate
use Constants, only: Zero, One, Three, Five, Seven, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: nR_Eff
real(kind=wp), intent(inout) :: Alpha(2), rm(2)
integer(kind=iwp), intent(in) :: iNQ
integer(kind=iwp) :: iANr, l_Max, mR
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iR
#endif
real(kind=wp) :: Dum(2,1), Radius_Max, RBS
logical(kind=iwp) :: Process
real(kind=wp), external :: Bragg_Slater, Eval_RMax

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*) 'nR,L_Quad=',nR,L_Quad
#endif
if (L_Quad > lMax_NQ) then
  call WarningMessage(2,'GenVoronoi: L_Quad > lMax_NQ')
  write(u6,*) 'Redimension lMax_NQ in nq_info'
  write(u6,*) 'lMax_NQ=',lMax_NQ
  write(u6,*) 'L_Quad=',L_Quad
  call Abend()
end if
l_Max = int(rm(1))
Radius_Max = Eval_RMax(Alpha(1),l_Max,rm(2))
#ifdef _DEBUGPRINT_
write(u6,*) 'Alpha(1)=',Alpha(1)
write(u6,*) 'l_max=',l_max
write(u6,*) 'rm(2)=',rm(2)
write(u6,*) 'Radius_Max=',Radius_Max
write(u6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate radial quadrature points. Observe that the integrand
! vanish at (r=0.0).

if (Quadrature == 'MHL') then

  iANr = NQ_Data(iNQ)%Atom_Nr
  RBS = Bragg_Slater(iANr)
  Alpha(1) = RBS
  mR = nR-1
  call mma_allocate(NQ_Data(iNQ)%R_Quad,2,mR,Label='R_Quad')
  NQ_Data(iNQ)%R_Quad(:,:) = Zero
  call GenRadQuad_MHL(NQ_Data(iNQ)%R_Quad,nR,nR_Eff,Alpha(1))
  call Truncate_Grid(NQ_Data(iNQ)%R_Quad,mR,nR_Eff,Radius_Max)
  mR = nR_Eff
  NQ_Data(iNQ)%R_max = NQ_Data(iNQ)%R_Quad(1,mR)

else if (Quadrature == 'LOG3') then

  rm(1) = Three
  ! alpha=5 (alpha=7 for alkali and rare earth metals)
  Alpha(1) = Five
  iANr = NQ_Data(iNQ)%Atom_Nr
  if ((iANr == 3) .or. (iANr == 4) .or. (iANr == 11) .or. (iANr == 12) .or. (iANr == 19) .or. (iANr == 20) .or. (iANr == 37) .or. &
      (iANr == 38) .or. (iANr == 55) .or. (iANr == 56) .or. (iANr == 87) .or. (iANr == 88)) Alpha(1) = Seven
  mR = nR-1
  call mma_allocate(NQ_Data(iNQ)%R_Quad,2,mR,Label='R_Quad')
  NQ_Data(iNQ)%R_Quad(:,:) = Zero
  call GenRadQuad_MK(NQ_Data(iNQ)%R_Quad,nR,nR_Eff,rm(1),Alpha(1))
  call Truncate_Grid(NQ_Data(iNQ)%r_Quad,mR,nR_Eff,Radius_Max)
  mR = nR_Eff
  NQ_Data(iNQ)%R_max = NQ_Data(iNQ)%R_Quad(1,mR)

else if (Quadrature == 'BECKE') then

  iANr = NQ_Data(iNQ)%Atom_Nr
  RBS = Bragg_Slater(iANr)
  if (iANr == 1) then
    Alpha(1) = RBS
  else
    Alpha(1) = Half*RBS
  end if
  mR = nR-1
  call mma_allocate(NQ_Data(iNQ)%R_Quad,2,mR,Label='R_Quad')
  NQ_Data(iNQ)%R_Quad(:,:) = Zero
  call GenRadQuad_B(NQ_Data(iNQ)%R_Quad,nR,nR_Eff,Alpha(1))
  call Truncate_Grid(NQ_Data(iNQ)%R_Quad,mR,nR_Eff,Radius_Max)
  mR = nR_Eff
  NQ_Data(iNQ)%R_max = NQ_Data(iNQ)%R_Quad(1,mR)

else if (Quadrature == 'TA') then

  Alpha(1) = -One
  iANr = NQ_Data(iNQ)%Atom_Nr
  if (iANr == 1) then
    Alpha(1) = 0.8_wp
  else if (iANr == 2) then
    Alpha(1) = 0.9_wp
  else if (iANr == 3) then
    Alpha(1) = 1.8_wp
  else if (iANr == 4) then
    Alpha(1) = 1.4_wp
  else if (iANr == 5) then
    Alpha(1) = 1.3_wp
  else if (iANr == 6) then
    Alpha(1) = 1.1_wp
  else if (iANr == 7) then
    Alpha(1) = 0.9_wp
  else if (iANr == 8) then
    Alpha(1) = 0.9_wp
  else if (iANr == 9) then
    Alpha(1) = 0.9_wp
  else if (iANr == 10) then
    Alpha(1) = 0.9_wp
  else if (iANr == 11) then
    Alpha(1) = 1.4_wp
  else if (iANr == 12) then
    Alpha(1) = 1.3_wp
  else if (iANr == 13) then
    Alpha(1) = 1.3_wp
  else if (iANr == 14) then
    Alpha(1) = 1.2_wp
  else if (iANr == 15) then
    Alpha(1) = 1.1_wp
  else if (iANr == 16) then
    Alpha(1) = 1.0_wp
  else if (iANr == 17) then
    Alpha(1) = 1.0_wp
  else if (iANr == 18) then
    Alpha(1) = 1.0_wp
  else if (iANr == 19) then
    Alpha(1) = 1.5_wp
  else if (iANr == 20) then
    Alpha(1) = 1.4_wp
  else if (iANr == 21) then
    Alpha(1) = 1.3_wp
  else if (iANr == 22) then
    Alpha(1) = 1.2_wp
  else if (iANr == 23) then
    Alpha(1) = 1.2_wp
  else if (iANr == 24) then
    Alpha(1) = 1.2_wp
  else if (iANr == 25) then
    Alpha(1) = 1.2_wp
  else if (iANr == 26) then
    Alpha(1) = 1.2_wp
  else if (iANr == 27) then
    Alpha(1) = 1.2_wp
  else if (iANr == 28) then
    Alpha(1) = 1.1_wp
  else if (iANr == 29) then
    Alpha(1) = 1.1_wp
  else if (iANr == 30) then
    Alpha(1) = 1.1_wp
  else if (iANr == 31) then
    Alpha(1) = 1.1_wp
  else if (iANr == 32) then
    Alpha(1) = 1.0_wp
  else if (iANr == 33) then
    Alpha(1) = 0.9_wp
  else if (iANr == 34) then
    Alpha(1) = 0.9_wp
  else if (iANr == 35) then
    Alpha(1) = 0.9_wp
  else if (iANr == 36) then
    Alpha(1) = 0.9_wp
  else
    call WarningMessage(2,'TA grid not defined')
    write(u6,*) ' TA grid not defined for atom number:',iANR
    call Abend()
  end if
  mR = nR-1
  call mma_allocate(NQ_Data(iNQ)%R_Quad,2,mR,Label='R_Quad')
  NQ_Data(iNQ)%R_Quad(:,:) = Zero
  call GenRadQuad_TA(NQ_Data(iNQ)%R_Quad,nR,nR_Eff,Alpha(1))
  call Truncate_Grid(NQ_Data(iNQ)%R_Quad,mR,nR_Eff,Radius_Max)
  mR = nR_Eff
  NQ_Data(iNQ)%R_max = NQ_Data(iNQ)%R_Quad(1,mR)

else if (Quadrature == 'LMG') then

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate radial quadrature. The first call will generate
  ! the size of the grid.

  nR = 1 ! Dummy size on the first call.
  Process = .false.
  call GenRadQuad_PAM(nR_Eff,rm,Alpha(1),Process,Dum,nR)

  nR = nR_Eff
  call mma_allocate(NQ_Data(iNQ)%R_Quad,2,nR,Label='R_Quad')
  NQ_Data(iNQ)%R_Quad(:,:) = Zero
  Process = .true.
  call GenRadQuad_PAM(nR_Eff,rm,Alpha(1),Process,NQ_Data(iNQ)%R_Quad,nR)
  NQ_Data(iNQ)%R_max = NQ_Data(iNQ)%R_Quad(1,nR)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# ifdef _DEBUGPRINT_
  write(u6,*) 'GenRadQuad_PAM ----> GenVoronoi'
  write(u6,*) 'nR_Eff=',nR_Eff
  write(u6,*) 'rm : ',rm(1),rm(2)
  write(u6,*) 'Alpha : ',Alpha(1),Alpha(2)
# endif
else
  call WarningMessage(2,'Invalid quadrature scheme:'//Quadrature)
  call Quit_OnUserError()
end if

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ******** The radial grid ********'
write(u6,*)
write(u6,*) 'Initial number of radial grid points=',nR
write(u6,*) 'iNQ=',iNQ
write(u6,*) 'Effective number of radial grid points=',nR_Eff
do iR=1,nR_Eff
  write(u6,*) NQ_Data(iNQ)%R_Quad(1,iR),NQ_Data(iNQ)%R_Quad(2,iR)
end do
write(u6,*)
write(u6,*) ' *********************************'
write(u6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine GenVoronoi
