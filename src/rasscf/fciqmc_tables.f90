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
! Copyright (C) 2014, Giovanni Li Manni                                *
!               2019, Oskar Weser                                      *
!***********************************************************************
module fciqmc_tables
  use stdalloc, only : mma_allocate, mma_deallocate
  implicit none
  private
  public :: FockTable, TwoElIntTable, OrbitalTable, mma_allocate, &
    mma_deallocate, length, print, fill_orbitals, fill_fock, fill_2ElInt, &
    cutoff_default, reorder, get_P_GAS, get_P_inp, unused

  type :: FockTable
    sequence
    real(kind=8), allocatable, dimension(:) :: values ! <i | F | j >
    integer, allocatable, dimension(:, :) :: index ! i, j
    real(kind=8) :: cutoff
    integer :: length
  end type FockTable

  type :: TwoElIntTable
    sequence
    real(kind=8), allocatable, dimension(:) :: values ! <ij| 1/r_{12} |kl>
    integer, allocatable, dimension(:, :) :: index ! i, j, k, l
    real(kind=8) :: cutoff
    integer :: length
  end type TwoElIntTable

  type :: OrbitalTable
    sequence
    real(kind=8), allocatable, dimension(:) :: values ! <i| F |i>
    integer, allocatable, dimension(:) :: index ! i
  end type OrbitalTable

  real(kind=8), parameter :: cutoff_default = 1.0d-11

  interface mma_allocate
    module procedure FockTable_allocate, TwoElIntTable_allocate, &
      OrbitalTable_allocate
  end interface

  interface mma_deallocate
    module procedure FockTable_deallocate, TwoElIntTable_deallocate, &
      OrbitalTable_deallocate
  end interface

  interface reorder
    module procedure FockTable_reorder, TwoElIntTable_reorder, &
      OrbitalTable_reorder, ALL_reorder
  end interface

  interface length
    module procedure FockTable_length, TwoElIntTable_length, &
      OrbitalTable_length
  end interface

  interface print
    module procedure FockTable_print, TwoElIntTable_print, OrbitalTable_print
  end interface

  interface unused
    module procedure FockTable_unused, TwoElIntTable_unused, OrbitalTable_unused
  end interface
  save
contains
  subroutine OrbitalTable_allocate(orbital_table, n)
    implicit none
    integer, intent(in) :: n
    type(OrbitalTable), intent(inout) :: orbital_table
    call mma_allocate(orbital_table%values, n)
    call mma_allocate(orbital_table%index, n)
  end subroutine OrbitalTable_allocate

  subroutine OrbitalTable_deallocate(orbital_table)
    implicit none
    type(OrbitalTable), intent(inout) :: orbital_table
    call mma_deallocate(orbital_table%values)
    call mma_deallocate(orbital_table%index)
  end subroutine OrbitalTable_deallocate

!>  @brief
!>    Fill in orbital energies
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The orbitals table gets filled with the orbital energies.
!>
!>  @param[in,out] orbitals Core
!>  @param[in] DIAF
!>  @param[in] iter
  subroutine fill_orbitals(orbitals, DIAF, iter)
    use general_data, only : nBas, nSym, nAsh, nFro, nIsh
    use rasscf_data, only : FDIAG

    implicit none
    type(OrbitalTable), intent(inout) :: orbitals
    integer, intent(in) :: iter
    real(kind=8), intent(in) :: DIAF(:)
    real(kind=8), allocatable :: EOrb(:)
    integer :: i, n, iSym, ioff, iDummy, itotnbas, Dummy, iErr
    character(*), parameter ::  FnInpOrb = 'INPORB'
    character(80) :: VecTit
    integer :: LuInpOrb = 10
    logical :: okay

    integer, parameter :: max_test = 20
    integer :: l_orb_test

    if (iter .eq. 1) then
      itotnbas = sum(nBas(:nSym))
      call f_Inquire(FnInpOrb,okay)
      if (okay) then
        call mma_allocate(EOrb, itotnbas)
        call RdVec(FnInpOrb,LuInpOrb,'E',nSym,nBas,nBas, &
          Dummy, Dummy, EOrb, iDummy, &
          VecTit, 0, iErr)
      else
        Write (6,*) 'RdCMO: Error finding MO file'
        call QTrace()
        call Abend()
      end if

      iOff = 0
      n = 1
      do iSym = 1, nSym
        if (nAsh(iSym) .gt. 0) then
          do i = 1, nAsh(isym)
            orbitals%index(n) = n
            orbitals%values(n) = EOrb(ioff + nFro(iSym) + nIsh(iSym) + i)
            n = n + 1
          enddo
        end if
        iOff   = iOff + nBas(iSym)
      end do
      call mma_deallocate(EOrb)
    else
      ioff = 0
      n = 1
      do iSym = 1, nSym
        do i = 1, nAsh(iSym)
          orbitals%index(n) = n
          orbitals%values(n) = DIAF(ioff + nFro(iSym) + nIsh(iSym) + i)
          n = n + 1
        end do
        ioff = ioff + nBas(iSym)
      end do
    end if

! ========== For testing purposes FROM HERE =============
    l_orb_test = min(max_test, length(orbitals))
    call Add_Info('Orbital Energy Input', &
      orbitals%values(:l_orb_test), l_orb_test, 8)
! ========== For testing purposes TO HERE ===============
  end subroutine fill_orbitals

  subroutine OrbitalTable_reorder(orbitals, P)
    type(OrbitalTable), intent(inout) :: orbitals
    integer, intent(in) :: P(:)
    integer :: i
    do i = 1, length(orbitals)
      orbitals%index(i) = P(orbitals%index(i))
    end do
  end subroutine

  pure integer function OrbitalTable_length(table)
    type(OrbitalTable), intent(in) :: table
    OrbitalTable_length = size(table%values)
  end function OrbitalTable_length

  subroutine OrbitalTable_print(table)
    type(OrbitalTable), intent(in) :: table
    integer :: i
    do i = 1, length(table)
      write(6, '(E15.7, I7)') table%values(i), table%index(i)
    end do
  end subroutine OrbitalTable_print

  subroutine OrbitalTable_unused(table)
    type(OrbitalTable), intent(in) :: table
    integer :: n
    if (.false.) n = length(table)
  end subroutine OrbitalTable_unused

  subroutine FockTable_allocate(fock_table, n)
    implicit none
    integer, intent(in) :: n
    type(FockTable), intent(inout) :: fock_table
    call mma_allocate(fock_table%values, n)
    call mma_allocate(fock_table%index, 2, n)
  end subroutine FockTable_allocate

  subroutine FockTable_deallocate(fock_table)
    implicit none
    type(FockTable), intent(inout) :: fock_table
    call mma_deallocate(fock_table%values)
    call mma_deallocate(fock_table%index)
  end subroutine FockTable_deallocate

!>  @brief
!>    Fill in fock matrix elements.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The fock_table gets filled with the Fock matrix elements
!>  whose absolute value is larger than cutoff.
!>  The values are given by
!>  \f[ < i | F | j > \f]
!>  The index is given by i and j.
!>
!>  @param[in,out] fock_table
!>  @param[in] Fock
!>  @param[in] EMY
!>  @param[in] cutoff Optional parameter that is set by default to
!>    fciqmc_tables::cutoff_default.
  subroutine fill_fock(fock_table, Fock, Emy, cutoff)
    use general_data, only : nActEl
    implicit none
    real(kind=8), intent(in) :: Fock(:), Emy
    type(FockTable), intent(inout) :: fock_table
    integer :: i, n, iOrb, jOrb
    real(kind=8) :: Emyn
    real(kind=8), optional :: cutoff
    real(kind=8) :: cutoff_

    integer, parameter :: max_test = 20
    integer :: l_fock_test

    cutoff_ = merge(cutoff, cutoff_default, present(cutoff))


    if (nActEl .ne. 0) then
      Emyn = Emy / dble(nActEl)
    else
      Emyn = 0.0d0
    end if

    n = 0
    do i = 1, size(Fock)
      if (abs(Fock(i)) .ge. cutoff_) then
        n = n + 1
        iOrb = ceiling(-0.5d0 + sqrt(2.0d0 * i))
        jOrb = i - (iOrb - 1) * iOrb / 2
        fock_table%index(:, n) = [iOrb, jOrb]
        if (iOrb .eq. jOrb) then
          fock_table%values(n) = Fock(i) - Emyn
        else
          fock_table%values(n) = Fock(i)
        end if
      end if
    end do
    fock_table%length = n
    fock_table%cutoff = cutoff_

! ========== For testing purposes FROM HERE =============
    l_fock_test = min(max_test, length(fock_table))
    call Add_Info('Fock element Input', Fock(:l_fock_test), l_fock_test, 8)
! ========== For testing purposes TO HERE ===============
  end subroutine fill_fock

  subroutine FockTable_reorder(fock, P)
    type(FockTable), intent(inout) :: fock
    integer, intent(in) :: P(:)
    integer :: i, j
    do j = 1, length(fock)
      do i = 1, 2
        fock%index(i, j) = P(fock%index(i, j))
      end do
    end do
  end subroutine

  pure integer function FockTable_length(table)
    type(FockTable), intent(in) :: table
    FockTable_length = table%length
  end function FockTable_length

  subroutine FockTable_print(table)
    type(FockTable), intent(in) :: table
    integer :: i, j
    do j = 1, length(table)
      write(6, '(E15.7, I7, I7)') table%values(j), (table%index(i, j), i=1, 2)
    end do
  end subroutine FockTable_print

  subroutine FockTable_unused(table)
    type(FockTable), intent(in) :: table
    integer :: n
    if (.false.) n = length(table)
  end subroutine FockTable_unused

  subroutine TwoElIntTable_allocate(table, n)
    implicit none
    integer, intent(in) :: n
    type(TwoElIntTable), intent(inout) :: table
    call mma_allocate(table%values, n)
    call mma_allocate(table%index, 4, n)
  end subroutine TwoElIntTable_allocate

  subroutine TwoElIntTable_deallocate(table)
    implicit none
    type(TwoElIntTable), intent(inout) :: table
    call mma_deallocate(table%values)
    call mma_deallocate(table%index)
  end subroutine TwoElIntTable_deallocate

!>  @brief
!>    Fill in two electron integrals
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The two_el_table gets filled with those two electron integrals
!>  whose absolute value is larger than cutoff.
!>  The values are given by
!>  \f[ < i, j | \frac{1}{r_{1,2}} | k, l > \f]
!>  The index is given by i, j, k, and l.
!>
!>  @param[in,out] two_el_table
!>  @param[in] TUVX
!>  @param[in] cutoff Optional parameter that is set by default to
!>    fciqmc_tables::cutoff_default.
  subroutine fill_2ElInt(two_el_table, TUVX, cutoff)
    implicit none
    real(kind=8), intent(in) :: TUVX(:)
    type(TwoElIntTable), intent(inout) :: two_el_table
    integer :: i, n, ijidx, klidx, iorb, jorb, korb, lorb
    real(kind=8), optional :: cutoff
    real(kind=8) :: cutoff_

    integer, parameter :: max_test = 20
    integer :: l_twoel_test

    cutoff_ = merge(cutoff, cutoff_default, present(cutoff))

    n = 0
    do i = 1, size(TUVX)
      if (abs(TUVX(i)) .ge. cutoff_) then
        n = n + 1
        ijidx = ceiling(-0.5d0 + sqrt(2.0d0 * i))
        klidx = i - (ijidx - 1) * ijidx / 2
        iorb = ceiling(-0.5d0 + sqrt(2.0d0 * ijidx))
        jorb = ijidx - (iorb - 1) * iorb / 2
        korb = ceiling(-0.5d0 + sqrt(2.0d0 * klidx))
        lorb = klidx - (korb - 1) * korb / 2
        two_el_table%index(:, n) = [iorb, jorb, korb, lorb]
        two_el_table%values(n) = TUVX(i)
      end if
    end do
    two_el_table%length = n
    two_el_table%cutoff = cutoff_

! ========== For testing purposes FROM HERE =============
    l_twoel_test = min(max_test, length(two_el_table))
    call Add_Info('TwoEl Integral element Input', &
      TUVX(:l_twoel_test), l_twoel_test, 8)
! ========== For testing purposes TO HERE ===============
  end subroutine fill_2ElInt

  subroutine TwoElIntTable_reorder(two_el_table, P)
    type(TwoElIntTable), intent(inout) :: two_el_table
    integer, intent(in) :: P(:)
    integer :: i, j
    do j = 1, length(two_el_table)
      do i = 1, 4
        two_el_table%index(i, j) = P(two_el_table%index(i, j))
      end do
    end do
  end subroutine

  pure integer function TwoElIntTable_length(table)
    type(TwoElIntTable), intent(in) :: table
    TwoElIntTable_length = table%length
  end function TwoElIntTable_length

  subroutine TwoElIntTable_print(table)
    type(TwoElIntTable), intent(in) :: table
    integer :: i, j
    do j = 1, length(table)
      write(6, '(E15.7, I7, I7, I7, I7)') &
        table%values(j), (table%index(i, j), i=1, 4)
    end do
  end subroutine TwoElIntTable_print

  subroutine TwoElIntTable_unused(table)
    type(TwoElIntTable), intent(in) :: table
    integer :: n
    if (.false.) n = length(table)
  end subroutine TwoElIntTable_unused

  function get_P_GAS(ngssh) result(P)
    use sorting, only : argsort
    implicit none
    integer, intent(in) :: ngssh(:, :)
    integer :: P(sum(ngssh)), X(sum(ngssh))
    integer :: iGAS, iSym, iOrb, bounds(2)
    bounds = shape(ngssh)
    iOrb = 1
    do iSym = 1, bounds(2)
      do iGAS = 1, bounds(1)
        X(iOrb : iOrb + ngssh(iGAS, iSym)) = iGAS
        iOrb = iOrb + ngssh(iGAS, iSym) + 1
      end do
    end do
    P = argsort(X)
  end function

  function get_P_inp(ReOrInp) result(P)
    use sorting, only : sort
    use general_data, only : nAsh
    implicit none
    integer, intent(in) :: ReOrInp(:)
    integer :: P(sum(nAsh)), change_idx(size(ReOrInp)), i
    P = [(i, i = 1, size(P))]
    change_idx = ReOrInp
    call sort(change_idx)
    P(change_idx) = ReOrInp
  end function

  subroutine ALL_reorder(orbitals, fock, two_el_table, P)
    type(OrbitalTable), intent(inout) :: orbitals
    type(FockTable), intent(inout) :: fock
    type(TwoElIntTable), intent(inout) :: two_el_table
    integer, intent(in) :: P(:)

    call reorder(orbitals, P)
    call reorder(fock, P)
    call reorder(two_el_table, P)
  end subroutine
end module fciqmc_tables
