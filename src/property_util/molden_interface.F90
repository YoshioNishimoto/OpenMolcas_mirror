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
! Copyright (C) 1999, Coen de Graaf                                    *
!               1999, Anders Bernhardsson                              *
!               1999, Roland Lindh                                     *
!***********************************************************************

subroutine Molden_Interface(iUHF,FName,filename)
!***********************************************************************
!                                                                      *
!     Object: to generate MOLDEN input file                            *
!                                                                      *
!                                                                      *
!     Authors: Coen de Graaf, Anders Bernardsson and R. Lindh, 1999    *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, MolWgh, nBas, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, lIrrep
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: iUHF
character(len=*), intent(in) :: FName, Filename
integer(kind=iwp) :: i, iAngMx_Valence, iatom, iB, ibas_lab(MxAtom), ic, iCntr, iCnttp, icontr, iD, iData, iDeg, iDummy(1), iErr, &
                     ii, iIrrep, ik, ipc, ipC2, ipC2_ab, ipCent, ipCent2, ipCent3, iPL, ipMull, ipp, ipPhase, iprim, ipV, ipV_ab, &
                     iRc = 0, iS, isegm, ishell, iv, iWFtype, j, jData, jPL, k, kk, kk_Max, l, Lu_, mAdCMO, mAdCMO_ab, mAdEor, &
                     mAdEor_ab, mAdOcc, mAdOcc_ab, mdc, MF, nAtom, nB, nData, nDeg, nOrb(8), nTest, nTot, nTot2
real(kind=wp) :: Check_CMO, Check_Energy, Check_Occupation, coeff, Coor(3,MxAtom), prim, r_Norm(maxbfn), Znuc(MxAtom)
logical(kind=iwp) :: Exists, Found, y_cart, y_sphere
character(len=LenIn) AtomLabel(MxAtom)
character(len=LenIn8+1) :: gtolabel(maxbfn)
character(len=100) :: Supername
character(len=40) :: VTitle
character(len=8) :: Env, MO_Label(maxbfn)
character(len=LenIn8), allocatable :: label(:)
real(kind=wp), parameter :: EorbThr = 50.0_wp
character, parameter :: shelllabel(7) = ['s','p','d','f','g','h','i'], &
                        cNumber(61) = ['1','2','3','4','5','6','7','8','9','0', &
                                       'a','b','c','d','e','f','g','h','i','j', &
                                       'k','l','m','n','o','p','q','r','s','t', &
                                       'u','v','w','x','y','z','A','B','C','D', &
                                       'E','F','G','H','I','J','K','L','M','N', &
                                       'O','P','Q','R','S','T','V','W','X','Y', &
                                       'Z']
integer(kind=iwp), external :: iPrintLevel
real(kind=wp), external :: DblFac
logical(kind=iwp), external :: Reduce_Prt
character(len=100), external :: Get_SuperName
!integer(kind=iwp), parameter :: MaxOrb_Molden = 400
#include "WrkSpc.fh"
! Statement function
integer(kind=iwp) :: ix, iy, iz
real(kind=wp) :: CC
CC(ix,iy,iz) = sqrt(DblFac(2*ix-1)*DblFac(2*iy-1)*DblFac(2*iz-1))

if (iRc == 1) return

! Do nothing within numerical_gradient
SuperName = Get_Supername()
if (SuperName == 'numerical_gradient') then
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the print level

iPL = iPrintLevel(-1)
jPL = iPL
if (Reduce_Prt() .and. (iPL < 3)) jPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
if (MolWgh == 1) then
  if (jPL >= 2) then
    write(u6,*) 'Molden_Interface: Unsupported normalization,Molwgh=1!'
  end if
  iRc = 1
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call GetEnvf('MOLCAS_MOLDEN',Env)
!if ((Env == ' ') .or. (Env == 'OFF')) then
if (Env == 'OFF') then
  if (jPL >= 2) then
    write(u6,*)
    write(u6,*) ' Molden files will not be produced'
    write(u6,*)
  end if
  iRC = 1
  return
end if
!VV: current version of Molden has no clear limit for MaxOrb
!if (MaxOrb > MaxOrb_Molden) then
!   if (jPL >= 2) then
!      write(u6,*)
!      write(u6,*) ' Molden_Interface: W A R N I N G !!!!'
!      write(u6,*)
!      write(u6,*) ' No Molden input file will be generated!'
!      write(u6,*)
!      write(u6,*) ' Calculation exceeds the max number of orbitals allowed for MOLDEN. To change this modify the'
!      write(u6,*) ' parameter MaxOrb_Molden in src/util/molden_interface.f and follow the instructions in Molden'
!      write(u6,*) ' on how to modify the parameter MaxOrb.'
!      write(u6,*)
!   end if
!   iRC = 1
!   return
!end if
!                                                                      *
!***********************************************************************
!                                                                      *
Check_CMO = Zero
Check_Energy = Zero
Check_Occupation = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
call f_Inquire('RUNFILE',Exists)
if (.not. Exists) then
  iRC = 1
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the characteristics of all different basis sets,
! provide each atom with a nuclear charge and establish
! a link between an atom and its basis set ---
!
! NOTICE!!!
! This call will also fill info.fh and the Basis_Info.

call Inter1(AtomLabel,iBas_Lab,Coor,Znuc,nAtom)
call Qpg_iArray('nOrb',Found,nData)
if (Found) then
  call Get_iArray('nOrb',nOrb,nData)
else
  call iCopy(nIrrep,nBas,1,nOrb,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
iAngMx_Valence = 0
do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%Aux) .and. (.not. dbsc(iCnttp)%Frag)) then
    nTest = dbsc(iCnttp)%nVal-1
    iAngMx_Valence = max(iAngMx_Valence,nTest)
  end if
end do
if (iAngMx_Valence > 4) then
  if (jPL >= 2) then
    write(u6,*) 'Sorry, Molden does not know how to handle'
    write(u6,*) 'functions with angular momentum larger than g'
  end if
  Go To 999
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Unnormalize contraction coefficients for the valence shells

do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%Aux) .and. (.not. dbsc(iCnttp)%Frag)) then
    do l=0,dbsc(iCnttp)%nVal-1
      ishell = dbsc(iCnttp)%iVal+l
      if (Shells(ishell)%Transf .and. (.not. Shells(iShell)%Prjct)) then
        if (jPL >= 2) then
          write(u6,*) 'Sorry, Molden does not support contaminants'
        end if
        Go To 999
      end if
      call Unnrmlz(Shells(ishell)%Exp,Shells(ishell)%nExp,Shells(ishell)%pCff,Shells(ishell)%nBasis,l)
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute memory requirements and allocate memory

nB = 0
do iS=0,nirrep-1
  nB = nB+nBas(is)
end do
call GetMem('ICENT','ALLO','INTE',ipCent,8*nB)
call GetMem('IPHASE','ALLO','INTE',ipPhase,8*nB)
call GetMem('nCENT','ALLO','INTE',ipCent2,nB)
call GetMem('ICENTER','ALLO','INTE',ipCent3,nB)
call GetMem('CMO2','ALLO','REAL',ipC2,nB**2)
call GetMem('VECTOR','ALLO','REAL',ipV,nB**2)
call dcopy_(nB**2,[Zero],0,Work(ipV),1)
if (iUHF == 1) then
  call GetMem('CMO2','ALLO','REAL',ipC2_ab,nB**2)
  call GetMem('VECTOR','ALLO','REAL',ipV_ab,nB**2)
  call dcopy_(nB**2,[Zero],0,Work(ipV_ab),1)
else
  ipC2_ab = ip_Dummy
  ipV_ab = ip_Dummy
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Open input file for MOLDEN

MF = 9
call molcas_open(MF,filename)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write the atom information in MOLDEN format to unit MF

y_cart = .false.
y_sphere = .false.
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) Go To 995
  do iCntr=1,dbsc(iCnttp)%nCntr
    do l=0,dbsc(iCnttp)%nVal-1
      ! Test for the appearance of cartesian functions with l=2,3,4
      ishell = dbsc(iCnttp)%iVal+l
      if ((l >= 2) .and. (.not. y_cart)) then
        if (.not. Shells(ishell)%Transf) y_cart = .true.
      end if
      if ((l >= 2) .and. (.not. y_sphere)) then
        if (Shells(ishell)%Transf) y_sphere = .true.
      end if
      if (y_sphere .and. y_cart) then
        if (jPL >= 2) then
          write(u6,*)
          write(u6,*) 'Failed to generate input file to MOLDEN'
          write(u6,*) 'No mixing allowed of spherical and cartesian d, f, g-functions'
        end if
        Go to 991
      end if
    end do
  end do
995 continue
end do
write(MF,'(A)') '[Molden Format]'
!                                                                      *
!***********************************************************************
!                                                                      *
! Write atomic information

write(MF,'(A)') '[N_Atoms]'
write(MF,*) natom
write(MF,'(A)') '[Atoms] (AU)'
do iatom=1,natom
  write(MF,99) AtomLabel(iatom),iatom,int(Znuc(iatom)),(coor(i,iatom),i=1,3)
99  format(A,2(3x,I4),5x,3F16.8,3I4)
end do
if (.not. y_cart) then
  if (S%iAngMx > 1) write(MF,'(A)') '[5D]'
  if (S%iAngMx > 2) write(MF,'(A)') '[7F]'
  if (S%iAngMx > 3) write(MF,'(A)') '[9G]'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out charges and dipole moments

call qpg_dArray('Mulliken Charge',Found,nData)
if (Found) then
  write(MF,'(A)') '[Charge] (Mulliken)'
  call Allocate_Work(ipMull,nData)
  call Get_dArray('Mulliken Charge',Work(ipMull),nData)

  iData = 0
  jData = 0
  do iCnttp=1,nCnttp
    if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%pChrg) Go To 775
    do iCntr=1,dbsc(iCnttp)%nCntr
      iData = iData+1
      mdc = iCntr+dbsc(iCnttp)%mdci
      nDeg = nIrrep/dc(mdc)%nStab
      do iDeg=1,nDeg
        jData = jData+1
        write(MF,*) Work(ipMull+iData-1)
      end do
    end do
775 continue
  end do
  call Free_Work(ipMull)
  if (iData /= nData) then
    write(u6,*) 'Molden_Interface: iData.ne.nData'
    write(u6,*) 'iData,nData=',iData,nData
    call Abend()
  end if
  if (jData /= nAtom) then
    write(u6,*) 'Molden_Interface: jData.ne.nAtom'
    write(u6,*) 'jData,nAtom=',jData,nAtom
    call Abend()
  end if
end if
!write(MF,'(A)') '[NDIPOLE]'
!write(MF,'(A)') '[DIPOLE]'
!                                                                      *
!***********************************************************************
!                                                                      *
! Write Gaussian basis set information to MOLDEN input file.

write(MF,'(A)') '[GTO] (AU)'

! Read exponents and contraction coefficients of each unique basis.
! Write the present basis set (iCnttp) to the molden.input file for
! the appropriate atoms.
! Moreover, a list is constructed which contains a label for each
! GTO (gtolabel). This list follows the MOLDEN order of GTOs.
! Later this list will be used in the transformation of sabf (the
! symmetry adapted basis functions).

iatom = 0
mdc = 0
kk = 0

do iCnttp=1,nCnttp             ! loop over unique basis sets
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) Go To 996

  do iCntr=1,dbsc(iCnttp)%nCntr  ! loop over sym. unique centers
    mdc = mdc+1
    nDeg = nIrrep/dc(mdc)%nStab
    do iDeg=1,nDeg             ! loop over centers
      iAtom = iAtom+1
      write(MF,'(I4)') iAtom

      do l=0,dbsc(iCnttp)%nVal-1
        ishell = dbsc(iCnttp)%iVal+l
        if (Shells(iShell)%nBasis > size(cNumber)) then
          write(u6,*) 'Interf: too many contracted functions!'
          write(u6,*) 'nBasis(iShell)=',Shells(iShell)%nBasis
          call Abend()
        end if

        ! Iterate over each contracted GTO

        do icontr=1,Shells(ishell)%nBasis

          ! Find the number of exponents with non-zero exponents

          isegm = 0
          do iprim=1,Shells(ishell)%nExp
            coeff = Shells(ishell)%pCff(iprim,icontr)
            if (coeff /= Zero) then
              isegm = isegm+1
            end if
          end do

          write(MF,'(3x,A1,I4)') shelllabel(l+1),isegm

          ! Write exponents and contraction coefficients.

          do iprim=1,Shells(ishell)%nExp
            coeff = Shells(ishell)%pCff(iprim,icontr)
            prim = Shells(ishell)%exp(iprim)
            if (coeff /= Zero) then
              write(MF,'(E17.9,E17.9)') prim,coeff
            end if
          end do

          ! Construction of gtolabel
          ! Molden order: for p-functions: px(1), py(1), pz(1),
          !                                px(2), py(2), pz(2), etc.
          ! for d-, and f-functions: similar

          if (l == 0) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'01s     '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
          end if
          if (l == 1) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02px    '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02py    '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02pz    '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
          end if
          if ((l == 2) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d00   '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d01+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d01-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d02+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d02-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
          end if
          if ((l == 2) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d020000 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000200 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000002 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,2)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d010100 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d010001 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,1)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000101 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,1)
            iWork(ipCent3+kk-1) = iAtom
          end if
          if ((l == 3) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f00   '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f01+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f01-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f02+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f02-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f03+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f03-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
          end if
          if ((l == 3) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f030000 '//cNumber(icontr)
            r_Norm(kk) = CC(3,0,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000300 '//cNumber(icontr)
            r_Norm(kk) = CC(0,3,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000003 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,3)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010200 '//cNumber(icontr)
            r_Norm(kk) = CC(1,2,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f020100 '//cNumber(icontr)
            r_Norm(kk) = CC(2,1,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f020001 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,1)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010002 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,2)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000102 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,2)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000201 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,1)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010101 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,1)
            iWork(ipCent3+kk-1) = iAtom
          end if
          if ((l == 4) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g00   '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g01+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g01-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g02+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g02-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g03+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g03-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g04+  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g04-  '//cNumber(icontr)
            r_Norm(kk) = One
            iWork(ipCent3+kk-1) = iAtom
          end if
          if ((l == 4) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g040000 '//cNumber(icontr)
            r_Norm(kk) = CC(4,0,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000400 '//cNumber(icontr)
            r_Norm(kk) = CC(0,4,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000004 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,4)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g030100 '//cNumber(icontr)
            r_Norm(kk) = CC(3,1,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g030001 '//cNumber(icontr)
            r_Norm(kk) = CC(3,0,1)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010300 '//cNumber(icontr)
            r_Norm(kk) = CC(1,3,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000301 '//cNumber(icontr)
            r_Norm(kk) = CC(0,3,1)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010003 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,3)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000103 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,3)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020200 '//cNumber(icontr)
            r_Norm(kk) = CC(2,2,0)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020002 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,2)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000202 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,2)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020101 '//cNumber(icontr)
            r_Norm(kk) = CC(2,1,1)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010201 '//cNumber(icontr)
            r_Norm(kk) = CC(1,2,1)
            iWork(ipCent3+kk-1) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010102 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,2)
            iWork(ipCent3+kk-1) = iAtom
          end if
        end do
      end do
      write(MF,'(A)') ' '
    end do
  end do
996 continue
end do
kk_Max = kk
if (nB > kk_max) then
  if (jPL >= 2) then
    write(u6,*) 'Molden_Interface: nB.gt.kk_max'
    write(u6,*) 'nB,kk_Max=',nB,kk_Max
  end if
  Go To 991
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nTot = 0
nTot2 = 0
do iS=0,nIrrep-1
  nTot = nTot+nBas(iS)
  nTot2 = nTot2+nBas(iS)**2
end do
call GetMem('Occ','Allo','Real',mAdOcc,nTot)
call GetMem('Eor','Allo','Real',mAdEor,nTot)
call GetMem('CMO','Allo','Real',mAdCMO,nTot2)
call FZero(Work(mAdOcc),nTot)
call FZero(Work(mAdEor),nTot)
call FZero(Work(mAdCMO),nTot2)
if (iUHF == 1) then
  call GetMem('Occ','Allo','Real',mAdOcc_ab,nTot)
  call GetMem('Eor','Allo','Real',mAdEor_ab,nTot)
  call GetMem('CMO','Allo','Real',mAdCMO_ab,nTot2)
  call FZero(Work(mAdOcc_ab),nTot)
  call FZero(Work(mAdEor_ab),nTot)
  call FZero(Work(mAdCMO_ab),nTot2)
else
  mAdOcc_ab = ip_Dummy
  mAdEor_ab = ip_Dummy
  mAdCMO_ab = ip_Dummy
end if

! Read HF CMOs from file

Lu_ = 75
call RdVec_(FName,Lu_,'COE',iUHF,nIrrep,nBas,nBas,Work(mAdCMO),Work(mAdCMO_ab),Work(mAdOcc),Work(mAdOcc_ab),Work(mAdEor), &
            Work(mAdEor_ab),iDummy,VTitle,1,iErr,iWFtype)

! Get the coeff. of sym adapted basis functions (ipC2)

call Dens_IF_SCF(Work(ipC2),Work(mAdCMO),'F')
call GetMem('CMO','Free','Real',mAdCMO,nTot2)
if (iUHF == 1) then
  call Dens_IF_SCF(Work(ipC2_ab),Work(mAdCMO_ab),'F')
  call GetMem('CMO','Free','Real',mAdCMO_ab,nTot2)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!  Back 'transformation' of the symmetry adapted basis functions.
!  Probably somewhat clumsy, but it seems to work.If someone
!  knows a more elegant way to do it, please improve this part!
!
!  PART 1: Obtain symmetry information (soout), construct a label
!          for each sabf, which will be used in part 2 to find the
!          corresponding GTO in the MOLDEN list by comparing with
!          gtolabel
!
!  nB       --- Total number of contracted basis functions
!  ipcent2  --- degeneracy of a basis function
!  ipCent   --- centres over which the basis function is
!               delocalized
!  ipPhase  --- phase of the AO in the linear combination

call mma_allocate(label,MaxBfn+MaxBfn_Aux,label='label')
call icopy(8*nB,[0],0,iWork(ipPhase),1)
call icopy(8*nB,[0],0,iWork(ipCent),1)
call SOout(label,iWork(ipCent),iWork(ipPhase))
ipc = 0
do iContr=1,nB
  iWork(ipCent2+iContr-1) = 0
  do k=1,8
    if (iWork(ipCent+ipc) /= 0) iWork(ipcent2+iContr-1) = iWork(ipCent2+iContr-1)+1
    ipc = ipc+1
  end do
  !vv this statement prevents overoptimization
  if (nB < -100) write(u6,*) iWork(ipCent2+iContr-1)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!
! Part 2: -Take a MOLCAS symmetry functions (loop i)
!         -Find the corresponding label in the MOLDEN list (loop j)
!         -Copy the coeff of the sabf in the MOLDEN MO (vector), multiply
!          by the appropriate factor (ipPhase),and divide by the number of
!          centres over which the sabf is delocalized (ipCent3).
!         -The vectors are copied by rows!

! loop over MOLCAS symmetry functions
i = 0
ik = 0
do iIrrep=0,nIrrep-1

  do iB=1,nBas(iIrrep)
    i = i+1

    if (iB == 1) then
      ik = 1
    else
      if (label(i-1) == label(i)) then
        ik = ik+1
      else
        ik = 1
      end if
    end if
    if (ik > size(cNumber)) then
      write(u6,*) 'Molden_Interface: ik > size(cNumber)'
      write(u6,*) 'ik,size(cNumber)=',ik,size(cNumber)
      write(u6,*) 'List of labels'
      do iD=1,nB
        write(u6,'(A)') Label(i)
      end do
      Go To 997
    end if

    write(MO_Label(i),'(I5,A3)') iB,lirrep(iIrrep)

    do j=1,nB

      if (gtolabel(j) == label(i)//cNumber(ik)) then
        do k=1,8
          ipc = (i-1)*8+k-1
          ipp = ipc
          if (iWork(ipCent+ipc) == iWork(ipcent3+j-1)) then
            do ii=1,nB
              ic = (ii-1)*nB+(i-1)
              iv = (ii-1)*nB+(j-1)
              if (MolWgh == 0) then
                Work(ipV+iv) = Work(ipV+iv)+(Work(ipC2+ic)*r_Norm(j))*real(iWork(ipPhase+ipp),kind=wp)/ &
                               real(iWork(ipcent2+i-1),kind=wp)
                if (iUHF == 1) Work(ipV_ab+iv) = Work(ipV_ab+iv)+(Work(ipC2_ab+ic)*r_Norm(j))*real(iWork(ipPhase+ipp),kind=wp)/ &
                                                 real(iWork(ipcent2+i-1),kind=wp)
              else
                Work(ipV+iv) = Work(ipV+iv)+(Work(ipC2+ic)*r_Norm(j))*real(iWork(ipPhase+ipp),kind=wp)/ &
                               sqrt(real(iWork(ipcent2+i-1),kind=wp))
                if (iUHF == 1) Work(ipV_ab+iv) = Work(ipV_ab+iv)+(Work(ipC2_ab+ic)*r_Norm(j))*real(iWork(ipPhase+ipp),kind=wp)/ &
                                                 sqrt(real(iWork(ipcent2+i-1),kind=wp))
              end if
            end do
          end if
        end do
      end if
    end do
  end do
end do
call mma_deallocate(label)
!                                                                      *
!***********************************************************************
!                                                                      *
!  Dump vector in the molden.input file

write(MF,'(A)') '[MO]'
ii = 0
do i=0,nB-1
  if (Work(mAdEOr+i) <= EorbThr) then
    write(MF,'(A,A)') 'Sym= ',MO_Label(i+1)
    write(MF,103) Work(mAdEOr+i)
    write(MF,'(A)') 'Spin= Alpha'
    write(MF,104) Work(mAdOcc+i)
    if (Work(mAdEOr+i) < Zero) then
      Check_Energy = Check_Energy+Work(mAdEOr+i)*real(i,kind=wp)
    end if
    Check_Occupation = Check_Occupation+Work(mAdOcc+i)*real(i,kind=wp)
    do j=1,nB
      write(MF,100) j,Work(ipV+ii+j-1)
      Check_CMO = Check_CMO+Work(ipV+ii+j-1)**2
    end do
  end if

  if (iUHF == 1) then
    if (Work(mAdEOr_ab+i) <= EorbThr) then
      write(MF,'(A,A)') 'Sym= ',MO_Label(i+1)
      write(MF,103) Work(mAdEOr_ab+i)
      write(MF,'(A)') 'Spin= Beta'
      write(MF,104) Work(mAdOcc_ab+i)
      Check_Energy = Check_Energy+Work(mAdEOr_ab+i)*real(i,kind=wp)
      Check_Occupation = Check_Occupation+Work(mAdOcc_ab+i)*real(i,kind=wp)
      do j=1,nB
        write(MF,100) j,Work(ipV_ab+ii+j-1)
        Check_CMO = Check_CMO+Work(ipV_ab+ii+j-1)**2
      end do
    end if
  end if

  ii = ii+nB
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!if (Env /= 'OFF') then
!  if (iUHF < 2) call Add_Info('MOLDEN_CMO',Check_CMO,1,2)
!  call Add_Info('MOLDEN_Occupation',Check_Occupation,1,2)
!  call Add_Info('MOLDEN_Energy',Check_Energy,1,2)
!end if
!                                                                      *
!***********************************************************************
!                                                                      *
do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%Aux) .and. (.not. dbsc(iCnttp)%Frag)) then
    do l=0,dbsc(iCnttp)%nVal-1
      ishell = dbsc(iCnttp)%iVal+l
      call Unnrmlz2(Shells(ishell)%Exp,Shells(ishell)%nExp,Shells(ishell)%pCff,Shells(ishell)%nBasis,l)
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (jPL >= 2) then
  write(u6,*)
  write(u6,*) ' Input file to MOLDEN was generated!'
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
997 continue
call GetMem('Eor','Free','Real',mAdEor,nTot)
call GetMem('Occ','Free','Real',mAdOcc,nTot)
if (iUHF == 1) then
  call GetMem('Eor','Free','Real',mAdEor_ab,nTot)
  call GetMem('Occ','Free','Real',mAdOcc_ab,nTot)
end if
991 continue
call GetMem('ICENT','FREE','INTE',ipCent,8*nB)
call GetMem('IPHASE','FREE','INTE',ipPhase,8*nB)
call GetMem('nCENT','FREE','INTE',ipCent2,nB)
call GetMem('ICENTER','FREE','INTE',ipCent3,nB)
if (iUHF == 1) then
  call GetMem('CMO2','FREE','REAL',ipC2_ab,nB**2)
  call GetMem('VECTOR','FREE','REAL',ipV_ab,nB**2)
end if
call GetMem('CMO2','FREE','REAL',ipC2,nB**2)
call GetMem('VECTOR','FREE','REAL',ipV,nB**2)
close(MF)
999 continue
call ClsSew()

! -------------FORMATS-------------
100 format(I4,3x,F16.8)
103 format('Ene= ',F10.4)
104 format('Occup= ',F10.5)

return

end subroutine Molden_Interface