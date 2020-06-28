************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICEnSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Oskar Weser                                      *
************************************************************************
#include "intent.h"
      module orthonormalization
        use definitions, only: wp
        use stdalloc, only : mma_allocate, mma_deallocate
        use fortran_strings, only : to_upper, str
        use blockdiagonal_matrices, only : t_blockdiagonal, new, delete,
     &    from_raw, to_raw, from_symm_raw, blocksizes
        use sorting, only : argsort
        use sorting_funcs, only : ge_r
        use linalg_mod, only: operator(.mult.), operator(.Tmult.),
     &      operator(.multT.), diagonalize, assert_

        implicit none
        save
        private
        public ::
     &    t_ON_scheme, ON_scheme, ON_scheme_values, orthonormalize

! TODO: Should be changed to default construction in the future.
! As of July 2019 the Sun and PGI compiler have problems.
        type :: t_ON_scheme_values
          integer :: no_ON, Gram_Schmidt, Lowdin, Canonical
        end type
        type(t_ON_scheme_values), parameter ::
! TODO: Dear fellow MOLCAS developer of the future:
! Please replace the following explicit constructur
! with the default constructor
!     instance = type()
! as soon as possible.
! As of July 2019 the Sun compiler requires explicit construction
! for parameter variables. (Which is wrong IMHO.)
     &    ON_scheme_values =
     &    t_ON_scheme_values(no_ON = 1, Gram_Schmidt = 2,
     &                       Lowdin = 3, Canonical = 4)

        type :: t_ON_scheme
          integer :: val = ON_scheme_values%Gram_Schmidt
        end type
        type(t_ON_scheme) :: ON_scheme = t_ON_scheme()


!>  @brief
!>    Orthonormalize a basis.
!>
!>  @author
!>    Oskar Weser
!>
!>  @details
!>  Reads in the overlap matrix and orthonormalizes accordingly.
!>
!>  @param[in] basis Can be raw memory i.e. a 1D double array
!>    or a t_blockdiagonal matrix.
!>    return type depends on input type.
!>  @param[in] scheme Optional argument. The orthonormalization scheme to use.
!>    The possibilities are Gram_Schmidt, Lowdin, Canonical, or no_ON
!>    (no_orthonormalization)
!>    as given by orthonormalization::ON_scheme_values.
!>    For a detailed explanation see \cite szabo_ostlund (p. 143).
        interface orthonormalize
          module procedure orthonormalize_raw, orthonormalize_blocks
        end interface

        interface Gram_Schmidt
          module procedure Gram_Schmidt_Array, Gram_Schmidt_Blocks
        end interface

        interface Lowdin
          module procedure Lowdin_Array, Lowdin_Blocks
        end interface

        interface Canonical
          module procedure Canonical_Array, Canonical_Blocks
        end interface

      contains

      subroutine orthonormalize_blocks(basis, scheme, ONB)
        use general_data, only : nSym, nBAs, nDel, nDelt, nSSH, nOrb
        use rasscf_data, only : nSec, nOrbt, nTot3, nTot4
        type(t_blockdiagonal), intent(in) :: basis(:)
        type(t_ON_scheme), intent(in) :: scheme
        type(t_blockdiagonal), intent(_OUT_) :: ONB(:)

        type(t_blockdiagonal), allocatable :: S(:)

        integer :: n_to_ON(nSym), n_new(nSym)

        call new(S, blocksizes=blocksizes(basis))
        call read_S(S)

        select case (scheme%val)
          case(ON_scheme_values%no_ON)
            continue
          case(ON_scheme_values%Lowdin)
            call Lowdin(basis, S, ONB)
          case(ON_scheme_values%Canonical)
            n_to_ON(:) = nBas(:nSym) - nDel(:nSym)
            call Canonical(basis, S, n_to_ON, ONB, n_new)
            call update_orb_numbers(n_to_ON, n_new,
     &          nDel, nSSH, nOrb, nDelt, nSec, nOrbt, nTot3, nTot4)
          case(ON_scheme_values%Gram_Schmidt)
            n_to_ON(:) = nBas(:nSym) - nDel(:nSym)
            call Gram_Schmidt(basis, S, n_to_ON, ONB, n_new)
            call update_orb_numbers(n_to_ON, n_new,
     &          nDel, nSSH, nOrb, nDelt, nSec, nOrbt, nTot3, nTot4)
        end select

        call delete(S)
      end subroutine orthonormalize_blocks

      subroutine orthonormalize_raw(CMO, scheme, ONB_v)
        use general_data, only : nBas, nSym
        real(wp), intent(in) :: CMO(:)
        type(t_ON_scheme), intent(in) :: scheme
        real(wp), intent(out) :: ONB_v(:)

        type(t_blockdiagonal), allocatable :: basis(:), ONB(:)

        call new(basis, blocksizes=nBAS(:nSym))
        call new(ONB, blocksizes=nBAS(:nSym))

        call from_raw(CMO, basis)
        call orthonormalize(basis, scheme, ONB)
        call to_raw(ONB, ONB_v)

        call delete(ONB)
        call delete(basis)
      end subroutine

! TODO: It would be nice, to use `impure elemental`
! instead of the manual looping.
! As of July 2019 some compilers don't support it.
      subroutine Lowdin_Blocks(basis, S, ONB)
        type(t_blockdiagonal), intent(in) :: basis(:), S(:)
        type(t_blockdiagonal), intent(_OUT_) :: ONB(:)

        integer :: i

        do i = 1, size(basis)
          call Lowdin(basis(i)%block, S(i)%block, ONB(i)%block)
        end do
      end subroutine Lowdin_Blocks


      subroutine Lowdin_Array(basis, S, ONB)
        real(wp), intent(in) :: basis(:, :), S(:, :)
        real(wp), intent(out) :: ONB(:, :)

        integer :: i
        real(wp), allocatable :: U(:, :), s_diag(:), X(:, :),
     &        tmp(:, :)

        call mma_allocate(U, size(S, 1), size(S, 2))
        call mma_allocate(X, size(S, 1), size(S, 2))
        call mma_allocate(tmp, size(S, 1), size(S, 2))
        call mma_allocate(s_diag, size(S, 2))

! Transform AO-overlap matrix S to the overlap matrix of basis.
! We search X that diagonalizes basis^T S basis
        call diagonalize(basis .Tmult. S .mult. basis , U, s_diag)

        call assert_(all(s_diag > 1.0d-10),
     &      "Linear dependency detected. "//
     &      "Lowdin can't cure it. Please use other ORTH keyword "//
     &      "from {Gram_Schmidt, Canonical}.")

! X = U s_diag^{-1/2} U^T
        do i = 1, size(tmp, 2)
          tmp(:, i) = U(:, i) / sqrt(s_diag(i))
        end do
        X(:, :) = tmp .multT. U
! With this X the overlap matrix has diagonal form.
! X^T basis^T S basis X = 1
! We finally have to convert to get the form:
! ONB^T S ONB = 1
        ONB = basis .mult. X

        call mma_deallocate(tmp)
        call mma_deallocate(s_diag)
        call mma_deallocate(X)
        call mma_deallocate(U)
      end subroutine Lowdin_Array


! TODO: It would be nice, to use `impure elemental`
! instead of the manual looping.
! As of July 2019 some compilers don't support it.
      subroutine Canonical_Blocks(basis, S, n_to_ON, ONB, n_new)
        type(t_blockdiagonal), intent(in) :: basis(:), S(:)
        integer, intent(in) :: n_to_ON(:)
        type(t_blockdiagonal), intent(_OUT_) :: ONB(:)
        integer, intent(out) :: n_new(:)

        integer :: i

        do i = 1, size(basis)
          call Canonical(basis(i)%block, S(i)%block, n_to_ON(i),
     &                   ONB(i)%block, n_new(i))
        end do
      end subroutine Canonical_Blocks


      subroutine Canonical_Array(basis, S, n_to_ON, ONB, n_new)
        real(wp), intent(in) :: basis(:, :), S(:, :)
        integer, intent(in) :: n_to_ON
        real(wp), intent(out) :: ONB(:, :)
        integer, intent(out) :: n_new

        logical :: lin_dep_detected
        integer :: i
        integer, allocatable :: idx(:)
        real(wp), allocatable :: U(:, :), s_diag(:),
     &      X(:, :), tmp(:, :)

        call mma_allocate(U, size(S, 1), size(S, 2))
        call mma_allocate(s_diag, size(S, 2))
        call mma_allocate(X, size(S, 1), size(S, 2))
        call mma_allocate(idx, size(S, 1))
        call mma_allocate(tmp, size(S, 1), size(S, 2))

! Transform AO-overlap matrix S to the overlap matrix of basis.
! We search X that diagonalizes basis^T S basis
        call diagonalize(basis .Tmult. S .mult. basis , U, s_diag)

        idx(:) = argsort(s_diag, ge_r)
        U(:, :) = U(:, idx)
        s_diag(:) = s_diag(idx)

        i = 0
        lin_dep_detected = .false.
        do while(.not. lin_dep_detected .and. i < n_to_ON)
          if (s_diag(i + 1) < 1.0d-10) then
            n_new = i
            lin_dep_detected = .true.
          end if
          i = i + 1
        end do
        if (.not. lin_dep_detected) n_new = n_to_ON

! X = U s_diag^{-1/2}
        do i = 1, n_new
          X(:, i) = U(:, i) / sqrt(s_diag(i))
        end do
! With this X the transformed overlap matrix has diagonal form.
! X^T basis^T S basis X = 1
! We finally have to convert to get the form:
! ONB^T S ONB = 1
        ONB(:, n_new + 1:) = basis(:, n_new + 1 :)
        ONB(:, :n_new) = basis .mult. X(:, :n_new)

        call mma_deallocate(tmp)
        call mma_deallocate(X)
        call mma_deallocate(idx)
        call mma_deallocate(s_diag)
        call mma_deallocate(U)
      end subroutine Canonical_Array

! TODO: It would be nice, to use `impure elemental`
! instead of the manual overloading.
! As of July 2019 some compilers don't support it.
      subroutine Gram_Schmidt_Blocks(basis, S, n_to_ON, ONB, n_new)
        type(t_blockdiagonal), intent(in) :: basis(:), S(:)
        integer, intent(in) :: n_to_ON(:)
        type(t_blockdiagonal), intent(_OUT_) :: ONB(:)
        integer, intent(out) :: n_new(:)

        integer :: i

        do i = 1, size(basis)
          call Gram_Schmidt(basis(i)%block, S(i)%block,
     &                      n_to_ON(i), ONB(i)%block, n_new(i))
        end do
      end subroutine Gram_Schmidt_Blocks



      subroutine Gram_Schmidt_Array(basis, S, n_to_ON, ONB, n_new)
        real(wp), intent(in) :: basis(:, :), S(:, :)
        integer, intent(in) :: n_to_ON
        real(wp), target, intent(out) :: ONB(:, :)
        integer, intent(out) :: n_new

        real(wp) :: L
        integer :: i, j
        logical :: lin_dep_detected, improve_solution

        real(wp), allocatable :: previous(:), correction(:), v(:)
        real(wp), pointer :: curr(:)

        call mma_allocate(previous, size(S, 1))
        call mma_allocate(correction, size(S, 1))
        call mma_allocate(v, size(S, 1))

        n_new = 0
        ONB(:, n_to_ON + 1 :) = basis(:, n_to_ON + 1 :)
        do i = 1, n_to_ON
          curr => ONB(:, n_new + 1)
          curr = basis(:, i)

          improve_solution = .true.
          lin_dep_detected = .false.
          do while (improve_solution .and. .not. lin_dep_detected)
            correction = 0._wp
            v(:) = S .mult. curr

            do j = 1, n_new
              correction(:) = correction(:)
     &                     + ONB(:, j) * dot_product(ONB(:, j), v)
            end do
            curr = curr - correction
            improve_solution = norm(correction, S=S) > 0.1_wp
            L = norm(curr, S=S)
            lin_dep_detected = L < 1.0e-10_wp
            if (.not. lin_dep_detected) then
              curr = curr / L
            end if
            if (.not. (improve_solution .or. lin_dep_detected)) then
              n_new = n_new + 1
            end if
          end do
        end do
        ONB(:, n_new + 1 : n_to_ON) = basis(:, n_new + 1 : n_to_ON)

        call mma_deallocate(v)
        call mma_deallocate(correction)
        call mma_deallocate(previous)
      end subroutine Gram_Schmidt_Array

      subroutine update_orb_numbers(
     &    n_to_ON, nNew,
     &    nDel, nSSH, nOrb, nDelt, nSec, nOrbt, nTot3, nTot4)
      use general_data, only : nSym
#include "warnings.fh"
#include "output_ras.fh"
      integer, intent(in) :: n_to_ON(:), nNew(:)
      integer, intent(inout) :: nDel(:), nSSH(:), nOrb(:),
     &  nDelt, nSec, nOrbt, nTot3, nTot4
      parameter(ROUTINE='update_orb_numbe')

      integer :: iSym, remove(nSym), total_remove

      call qEnter(ROUTINE)

      remove = n_to_ON(:nSym) - nNew(:nSym)
      total_remove = sum(remove(:nSym))
      if (total_remove > 0) then
        do iSym = 1, nSym
          if (nSSH(iSym) < remove(iSym)) then
            call WarningMessage(2,'Orthonormalization Error')
            Write(LF,*) 'Exact or very near linear dependence '
            Write(LF,*) 'forces RASSCF to stop execution.'
            Write(LF,*) 'Symmetry block:', iSym
            Write(LF,*) 'Effective NR of orthonormal orbs:', nNew(iSym)
            Write(LF,*) 'Earlier number of deleted orbs:', nDel(iSym)
            Write(LF,*) 'Earlier number of secondary orbs:', nSSH(iSym)
            Write(LF,*) 'New number of deleted orbs:',
     &                  nDel(iSym) + remove(iSym)
            Write(LF,*) 'New number of secondary orbs:',
     &                  nSSH(iSym) - remove(iSym)
            call quit(_RC_GENERAL_ERROR_)
          else if (iPRLoc(1) >= usual) then
            call WarningMessage(1,'Orthonormalization Warning')
            Write(LF,*) 'Exact or very near linear dependence'
            Write(LF,*) 'forces RASSCF to delete additional orbitals.'
            Write(LF,*) 'Symmetry block:', iSym
            Write(LF,*) 'Earlier number of deleted orbs =', nDel(iSym)
            Write(LF,*) 'New number of deleted orbs =',
     &                  nDel(iSym) + remove(iSym)
          end if
        end do
        nDel(:nSym) = nDel(:nSym) + remove(:nSym)
        nSSH(:nSym) = nSSH(:nSym) - remove(:nSym)
        nOrb(:nSym) = nOrb(:nSym) - remove(:nSym)
        nDelt = nDelt + total_remove
        nSec = nSec - total_remove
        nOrbt = nOrbt - total_remove
        nTot3 = sum((nOrb(:nSym) + nOrb(:nSym)**2) / 2)
        nTot4 = sum(nOrb(:nSym)**2)
      end if
      call qExit(ROUTINE)
      end subroutine update_orb_numbers


      subroutine read_raw_S(S_buffer)
        real(wp), intent(inout) :: S_buffer(:)
        integer :: i_Rc, i_Opt, i_Component, i_SymLbl
#include "warnings.fh"

        i_Rc = 0
        i_Opt = 2
        i_Component = 1
        i_SymLbl = 1
        Call RdOne(i_Rc, i_Opt, 'Mltpl  0', i_Component,
     &             S_buffer, i_SymLbl)
        if ( i_rc /= 0 ) then
          write(6,*)' RASSCF is trying to orthonormalize orbitals but'
          write(6,*)' could not read overlaps from ONEINT. Something'
          write(6,*)' is wrong with the file, or possibly with the'
          write(6,*)' program. Please check.'
          call quit(_RC_IO_ERROR_READ_)
        end if
      end subroutine


      subroutine read_S(S)
        use general_data, only : nBas, nSym, nActEl
        use rasscf_data, only : nFr, nIn, Tot_Nuc_Charge
#include "warnings.fh"
#include "output_ras.fh"
        type(t_blockdiagonal) :: S(nSym)

        parameter(ROUTINE='update_orb_numbe')
        integer :: size_S_buffer
        real(wp) :: Mol_Charge
        real(wp), allocatable :: S_buffer(:)

        call qEnter(ROUTINE)

        size_S_buffer = sum(nBas(:nSym) * (nBas(:nSym) + 1) / 2)
        call mma_allocate(S_buffer, size_S_buffer + 4)
        call read_raw_S(S_buffer)
        Tot_Nuc_Charge = S_buffer(size_S_buffer + 4)
        call from_symm_raw(S_buffer, S)
        call mma_deallocate(S_buffer)

        Mol_Charge = Tot_Nuc_Charge - dble(2 * (nFr + nIn) + nActEl)
        call put_dscalar('Total Charge    ', Mol_Charge)
        if (IPRLOC(1) >= usual) then
          write(6,*)
          write(6,'(6x,A,f8.2)') 'Total molecular charge',Mol_Charge
        end if

        call qExit(ROUTINE)
      end subroutine

!>
!>  @brief
!>    Calculates v1^T S v2.
!>
!>  @author
!>    Oskar Weser
!>
!>  @details
!>  Calculates \f[ v1^T S v2 \f]
!>  S has to be a symmetric positive definite matrix, which is not
!>  tested.
!>  If S is ommited, it defaults to the unit matrix,
!>  i.e. the Euclidean dot-product.
      function dot_product_(v1, v2, S) result(dot)
        real(wp), intent(in) :: v1(:), v2(:)
        real(wp), intent(in), optional :: S(:, :)
        real(wp) :: dot

        if (present(S)) then
          dot = dot_product(v1, S .mult. v2)
        else
          dot = dot_product(v1, v2)
        end if
      end function

      function norm(v, S) result(L)
        real(wp), intent(in) :: v(:)
        real(wp), intent(in), optional :: S(:, :)
        real(wp) :: L
        if (present(S)) then
          L = sqrt(dot_product_(v, v, S))
        else
! One could use norm2 here, but Sun and PGI compilers don't know this.
          L = sqrt(sum(v**2))
        end if
      end function

      end module orthonormalization
