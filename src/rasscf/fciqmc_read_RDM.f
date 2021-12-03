************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2016,2017, Giovanni Li Manni                           *
*               2019-2021, Oskar Weser                                 *
*               2021, Werner Dobrautz                                  *
************************************************************************

      module fciqmc_read_RDM
      use fortran_strings, only: str
      use definitions, only: wp, u6
      use stdalloc, only: mma_allocate, mma_deallocate
      use para_info, only: myRank
      use rasscf_data, only : NRoots
      use general_data, only : nActEl
      ! Note that two_el_idx_flatten has also out parameters.
      use index_symmetry, only : two_el_idx_flatten
      use CI_solver_util, only: CleanMat
      use linalg_mod, only: abort_, verify_

      implicit none

      private
      public :: read_neci_RDM, cleanup
      contains

!>  @brief
!>    Read NECI RDM files
!>
!>  @author Werner Dobrautz, Oskar Weser
!>
!>  @details
!>  Read the spin-free TwoRDM file written by the GUGA-NECI
!>  implementation and transfer them to Molcas.
!>  The spin density is set to zero, because spin projection
!>  is not properly defined in the GUGA framework.

!>  @paramin[out] DMAT Average spin-free 1 body density matrix
!>  @paramin[out] DSPN spin-dependent 1-RDM (set to zero)
!>  @paramin[out] PSMAT Average spin-free 2 body density matrix
!>  @paramin[out] PAMAT 'fake' Average antisymm. 2-dens matrix


      subroutine read_neci_RDM(
     &    iroot, weight, tGUGA, DMAT, DSPN, PSMAT, PAMAT
     &  )
          ! wrapper around `read_single_neci_(GUGA)_RDM` to average
          ! normal and GUGA density matrices for stochastic SA-RASSCF.
          integer(wp), intent(in) :: iroot(:)
          real(wp), intent(in) :: weight(:)
          logical, intent(in) :: tGUGA
          real(wp), intent(out) :: DMAT(:), DSPN(:), PSMAT(:), PAMAT(:)
          integer :: i, j
          real(wp), allocatable :: temp_DMAT(:), temp_DSPN(:),
     &                             temp_PSMAT(:), temp_PAMAT(:)

          ! prevent stackoverflow
          call mma_allocate(temp_DMAT, size(DMAT))
          call mma_allocate(temp_DSPN, size(DSPN))
          call mma_allocate(temp_PSMAT, size(PSMAT))
          call mma_allocate(temp_PAMAT, size(PAMAT))

          DMAT = 0.0_wp; DSPN = 0.0_wp
          PSMAT = 0.0_wp; PAMAT = 0.0_wp

          do i = 1, NRoots
              do j = 1, size(iroot)
                  if (iroot(j) == i) then
                      if (tGUGA) then
                          call read_single_neci_GUGA_RDM(
     &                        iroot(j), temp_DMAT, temp_DSPN,
     &                        temp_PSMAT, temp_PAMAT
     &                    )
                      else
                          call read_single_neci_RDM(
     &                        iroot(j), temp_DMAT, temp_DSPN,
     &                        temp_PSMAT, temp_PAMAT
     &                    )
                      end if

                      DMAT = DMAT + weight(j) * temp_DMAT
                      DSPN = DSPN + weight(j) * temp_DSPN
                      PSMAT= PSMAT + weight(j) * temp_PSMAT
                      PAMAT= PAMAT + weight(j) * temp_PAMAT
                  end if
              end do
          end do
          call mma_deallocate(temp_DMAT)
          call mma_deallocate(temp_DSPN)
          call mma_deallocate(temp_PSMAT)
          call mma_deallocate(temp_PAMAT)
      end subroutine read_neci_RDM


      subroutine read_single_neci_GUGA_RDM(
     &  iroot, DMAT, DSPN, PSMAT, PAMAT
     & )
          integer, intent(in) :: iroot
          real(wp), intent(out) :: DMAT(:), DSPN(:), PSMAT(:), PAMAT(:)
          integer :: file_id, isfreeunit, i
          logical :: tExist
          real(wp) :: RDMval

          if (myRank /= 0) then
              call bcast_2RDM("PSMAT." // str(iroot))
              call bcast_2RDM("PAMAT." // str(iroot))
              call bcast_2RDM("DMAT."  // str(iroot))
          end if

          call f_Inquire('PSMAT.' // str(iroot), tExist)
          call verify_(tExist, 'PSMAT.' // str(iroot) //
     &                 ' does not exist')
          call f_Inquire('PAMAT.' // str(iroot), tExist)
          call verify_(tExist, 'PAMAT.' // str(iroot) //
     &                 ' does not exist')
          call f_Inquire('DMAT.'  // str(iroot), tExist)
          call verify_(tExist, 'DMAT.'  // str(iroot) //
     &                 ' does not exist')

          PSMAT(:) = 0.0_wp; PAMAT(:) = 0.0_wp;
          DMAT(:) = 0.0_wp; DSPN(:) = 0.0_wp;

          file_id = IsFreeUnit(11)
          call Molcas_Open(file_id, 'PSMAT.' // str(iroot))
              do while (read_line(file_id, i, RDMval))
                psmat(i) = RDMval
              end do
          close(file_id)

          file_id = IsFreeUnit(11)
          call Molcas_Open(file_id, 'PAMAT.' // str(iroot))
              do while (read_line(file_id, i, RDMval))
                pamat(i) = RDMval
              end do
          close(file_id)

          file_id = IsFreeUnit(11)
          call Molcas_Open(file_id, 'DMAT.'  // str(iroot))
              do while (read_line(file_id, i, RDMval))
                dmat(i) = RDMval
              end do
          close(file_id)

          ! Clean evil non-positive semi-definite matrices,
          ! by clamping the occupation numbers between 0 and 2.
          ! DMAT is intent(inout)
          call cleanMat(DMAT)

      contains

          !> Read a line from GUGA NECI RDM file with
          !> `file_id` as unit and parse it into i and RDMval.
          !> Return true, if line was successfully read and false
          !> if EOF reached.
          !> Aborts if IO error happens.
          !> **i and RDMval are undefined, if functions returns false**
          logical function read_line(file_id, i, RDMval)
              ! changed variable names to prevent masking the parent
              ! scope
              integer, intent(in) :: file_id
              integer, intent(out) :: i
              real(wp), intent(out) :: RDMval
              integer :: iread
              read_line = .false.
              read(file_id, "(I6,G25.17)", iostat=iread) i, RDMval
              if (iread > 0) then
                  call abort_('Error in read_next')
              else if (is_iostat_end(iread)) then
                  ! Let's try to cause an error, if code uses undefined
                  ! values for i and RDMval
                  i = huge(i)
                  RDMval = huge(RDMval)
              else
                  read_line = .true.
              end if
          end function

      end subroutine read_single_neci_GUGA_RDM

!>  @brief
!>    Start and control FCIQMC.
!>
!>  @author Giovanni Li Manni, Oskar Weser
!>
!>  @details
!>  Read TwoRDM files written by NECI and transfer them to Molcas.
!>  Neci can have some intermediate spin-resolved/spin-free RDMs where
!>  basically aaaa contains average of aaaa and bbbb, abab contains
!>  average of abab and baba... This is ok for CASSCF but not ok
!>  for spin-resolved properties, in which case the completely
!>  spin-resolved RDMs need to be read-in. In principle, NECI could
!>  also evaluate and store completely spin-free matrices. In that
!>  case only a reordering following Molcas convention is necessary.
!>
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] DSPN Average spin 1-dens matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
!>
!>  This function is unmaintainable, someone has to rewrite it.
      subroutine read_single_neci_RDM(iroot, DMAT, DSPN, PSMAT, PAMAT)
          use Para_Info, only: MyRank
#include "output_ras.fh"
          integer, intent(in) :: iroot
          real*8, intent(out) :: DMAT(:), DSPN(:), PSMAT(:), PAMAT(:)
          integer :: iUnit, isfreeunit, p, q, r, s, pq, rs, ps, rq,
     &               psrq, pqrs, iread, norb, iprlev
          logical :: tExist, switch
          real*8 :: fac, RDMval, fcnacte
          real*8 :: D_alpha(size(DMAT)), D_beta(size(DMAT))

          iprlev = iprloc(1)
          if(iprlev == debug) then
            write(u6,*) 'Rank of process: ', MyRank
          end if
        ! TODO: Does it really matter, can we not just read all spin
        ! density matrices? Currently, this just adds another level of
        ! unnecessary nesting.
          switch = .true.  ! spin-resolved

          if(myRank /= 0) then
            call bcast_2RDM("TwoRDM_aaaa." // str(iroot))
            call bcast_2RDM("TwoRDM_abab." // str(iroot))
            call bcast_2RDM("TwoRDM_abba." // str(iroot))
            call bcast_2RDM("TwoRDM_bbbb." // str(iroot))
            call bcast_2RDM("TwoRDM_baba." // str(iroot))
            call bcast_2RDM("TwoRDM_baab." // str(iroot))
          end if
          call f_Inquire('TwoRDM_aaaa.' // str(iroot),tExist)
          call verify_(tExist, 'TwoRDM_aaaa.' // str(iroot) //
     &                 ' does not exist')
          call f_Inquire('TwoRDM_abab.'// str(iroot),tExist)
          call verify_(tExist, 'TwoRDM_abab.' // str(iroot) //
     &                 ' does not exist')
          call f_Inquire('TwoRDM_abba.'// str(iroot),tExist)
          call verify_(tExist, 'TwoRDM_abba.' // str(iroot) //
     &                 ' does not exist')
          if(switch) then
            call f_Inquire('TwoRDM_bbbb.' // str(iroot),tExist)
            call verify_(tExist, 'TwoRDM_bbbb.' // str(iroot) //
     &                   ' does not exist')
            call f_Inquire('TwoRDM_baba.' // str(iroot),tExist)
            call verify_(tExist, 'TwoRDM_baba.' // str(iroot) //
     &                   ' does not exist')
            call f_Inquire('TwoRDM_baab.' // str(iroot),tExist)
            call verify_(tExist, 'TwoRDM_baab.' // str(iroot) //
     &                   ' does not exist')
          end if

        D_alpha(:) = 0.0d0
        D_beta(:) = 0.0d0
        PSMAT(:) = 0.0d0
        PAMAT(:) = 0.0d0

        fac = merge(0.5d0, 1.0d0, switch)
        fcnacte = 1.0d0 / dble(nactel - 1)

        iUnit = IsFreeUnit(11)
        call Molcas_Open(iUnit, 'TwoRDM_aaaa.' // str(iroot))
        Rewind(iUnit)
        IF(IPRLEV >= DEBUG) THEN
          write(u6,*) 'p   q   r   s  pq  rs  pqrs',
     &    'RDMval                PSMAT                   PAMAT'
          write(u6,*) '******* AAAA *******'
        end if
        do

            read(iUnit, "(4I6,G25.17)", iostat=iread)
     &          s, q, r, p, RDMval
            if(iread /= 0) exit
            pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
            ! Contribution to PSMAT and PAMAT:
            PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
            if (r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs)
     &          + fac * RDMval
            if (p > q.and.r > s) PAMAT(pqrs) = PAMAT(pqrs)
     &          + fac * RDMval
            if (p > q.and.r < s) PAMAT(pqrs) = PAMAT(pqrs)
     &          - fac * RDMval
            IF(IPRLEV >= DEBUG) THEN
               write(u6,'(7I6,3G25.17)')
     &         p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
            END IF
            ! Contribution to D_alpha (not final):
            if (p == q) D_alpha(rs) = D_alpha(rs) + RDMval
            if (r == s) D_alpha(pq) = D_alpha(pq) + RDMval
            psrq = two_el_idx_flatten(p, s, r, q, ps, rq)
            ! Contribution to PSMAT and PAMAT:
            if (r <= q) then
              PSMAT(psrq) = PSMAT(psrq) - fac * RDMval
              if (r /= q) PAMAT(psrq) = PAMAT(psrq) + fac * RDMval
            end if
            if (r > q) then
              PSMAT(psrq) = PSMAT(psrq) - fac * RDMval
              PAMAT(psrq) = PAMAT(psrq) - fac * RDMval
            end if
            IF(IPRLEV >= DEBUG) THEN
              write(u6,'(7I6,3G25.17)')
     &          p,s,r,q,ps,rq,psrq, RDMval,PSMAT(psrq),PAMAT(psrq)
            END IF
            ! Contribution to D_alpha (not final): The minus sign comes
            ! from the fact that in NECI these elements have opposite
            ! sign compared to the element in normal order, that is
            ! d_pqrs = -d_psrq.
            if (p == s) D_alpha(rq) = D_alpha(rq) - RDMval
            if (r == q) D_alpha(ps) = D_alpha(ps) - RDMval
        end do
        close(iunit)

        ! Processing TwoRDM-BBBB
        if (switch) then
          iUnit=IsFreeUnit(11)
          Call Molcas_Open(iUnit,'TwoRDM_bbbb.' // str(iroot))
          Rewind(iUnit)
          if (IPRLEV >= DEBUG) then
           write(u6,*) '  p   q   r   s  pq  rs  pqrs',
     &     'RDMval                PSMAT                   PAMAT'
           write(u6,*) '******* BBBB *******'
          end if
          do
            ! processing as PQRS
            read(iUnit,"(4I6,G25.17)",iostat=iread) s, q, r, p, RDMval
            if(iread /= 0) exit
            pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
            ! Contribution to PSMAT and PAMAT:
            PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
            if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
            if(p > q.and.r > s) PAMAT(pqrs) = PAMAT(pqrs) + fac * RDMval
            if(p > q.and.r < s) PAMAT(pqrs) = PAMAT(pqrs) - fac * RDMval
            IF(IPRLEV >= DEBUG) THEN
              write(u6,'(7I6,3G25.17)')
     &           p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
            END IF
            ! Contribution to D_beta (not final):
            if (p == q) D_beta(rs) = D_beta(rs) + RDMval
            if (r == s) D_beta(pq) = D_beta(pq) + RDMval
            ! processing as PSRQ
            psrq = two_el_idx_flatten(p, s, r, q, ps, rq)
            ! Contribution to PSMAT and PAMAT:
            if(r <= q) then
              PSMAT(psrq) = PSMAT(psrq) - fac*RDMval
              if(r /= q) PAMAT(psrq) = PAMAT(psrq) + fac*RDMval
            end if
            if(r > q) then
              PSMAT(psrq) = PSMAT(psrq) - fac*RDMval
              PAMAT(psrq) = PAMAT(psrq) - fac*RDMval
            end if
            IF(IPRLEV >= DEBUG) THEN
              write(u6,'(7I6,3G25.17)')
     &          p,s,r,q,ps,rq,psrq, RDMval,PSMAT(psrq),PAMAT(psrq)
            END IF
            ! Contribution to D_beta (not final): The minus sign comes
            ! from the fact that in NECI these elements have opposite
            ! sign compared to the element in normal order, that is
            ! d_pqrs = -d_psrq.
              if(p == s) D_beta(rq)=D_beta(rq)-RDMval
              if(r == q) D_beta(ps)=D_beta(ps)-RDMval
            end do
            close(iunit)
        end if ! End statement for spin-resolved RDMs.
        ! Processing TwoRDM-ABAB
        iUnit=IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_abab.' // str(iroot))
        Rewind(iUnit)
        IF(IPRLEV >= DEBUG) THEN
         write(u6,*) '******* ABAB *******'
        END IF
        do
          read(iUnit,"(4I6,G25.17)",iostat=iread) s,q,r,p,RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
        ! Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
        if(r > s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
        if(r < s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
        if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
        if (IPRLEV >= DEBUG) then
          write(u6,'(7I6,3G25.17)')
     &        p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        end if
        ! Contribution to D_alpha and D_beta (not final):
        if (p == q) D_alpha(rs) = D_alpha(rs) + RDMval
        if (r == s .and. p /= r) D_beta(pq) = D_beta(pq) + RDMval
      end do

        ! Copy D_beta to D_alpha and clean D_beta again for further use:
      if (.not. switch) then
        D_alpha(:) = D_beta(:) + D_alpha(:)
        D_beta(:) = 0.0d0
      end if
      close(iunit)
      ! Processing TwoRDM-BABA
      if (switch) then
        iUnit = IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_baba.' // str(iroot))
        rewind(iUnit)
        if (IPRLEV >= DEBUG) then
           write(u6,*) '******* BABA *******'
        end if
        do
          read(iUnit,"(4I6,G25.17)",iostat=iread) s,q,r,p,RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
          ! Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
          if(r > s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
          if(r < s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
          if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
          IF(IPRLEV >= DEBUG) THEN
            write(u6,'(7I6,3G25.17)')
     &          p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
          ! Contribution to D_alpha (not final):
          if(p == q) D_beta(rs) = D_beta(rs) + RDMval
          if(r == s.and.p /= r) D_alpha(pq) = D_alpha(pq)+RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
      ! Processing TwoRDM-ABBA
      iUnit=IsFreeUnit(11)
      Call Molcas_Open(iUnit,'TwoRDM_abba.' // str(iroot))
      Rewind(iUnit)
      IF(IPRLEV >= DEBUG) THEN
        write(u6,*) '******* ABBA *******'
      END IF
      do
        read(iUnit,"(4I6,G25.17)",iostat=iread) q,s,r,p,RDMval
        if(iread /= 0) exit
        pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
        ! Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) - fac*RDMval
        if(r < s) then
          PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
        end if
        if(r > s) then
          PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
        end if
        IF(IPRLEV >= DEBUG) THEN
          write(u6,'(7I6,3G25.17)')
     &        p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        END IF
          ! Contribution to D_alpha (not final):
          if(r == s) D_alpha(pq)=D_alpha(pq)-RDMval
          if(p == q) D_beta(rs)=D_beta(rs)-RDMval
      end do
      if (.not. switch) then
        D_alpha(:) = D_beta(:) + D_alpha(:)
        D_beta(:) = 0.0d0
      end if
      close(iunit)
      ! Processing TwoRDM-BAAB
      if(switch) then
        iUnit=IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_baab.' // str(iroot))
        Rewind(iUnit)
        IF(IPRLEV >= DEBUG) THEN
          write(u6,*) ' ******* BAAB *******'
        END IF
        do
          read(iUnit,"(4I6,G25.17)",iostat=iread) q,s,r,p,RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
          ! Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) - fac*RDMval
          if(r < s) then
            PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
          end if
          if(r > s) then
            PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
          end if
          IF(IPRLEV >= DEBUG) THEN
            write(u6,'(7I6,3G25.17)')
     &          p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
          ! Contribution to D_alpha (not final):
          if(p == q) D_alpha(rs)=D_alpha(rs)-RDMval
          if(r == s) D_beta(pq)=D_beta(pq)-RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
      ! Final Updates to RDMs
      if (.not.switch) D_beta(:) = D_alpha(:)
      D_alpha(:) = fcnacte * D_alpha(:)
      D_beta(:) = fcnacte * D_beta(:)
      DSPN(:) = D_Beta(:) - D_alpha(:)
      DMAT(:) = D_Beta(:) + D_alpha(:)

      ! Clean evil non-positive semi-definite matrices. DMAT is input
      ! and output.
      call cleanMat(DMAT)

      IF(IPRLEV >= DEBUG) THEN
          norb  = (int(sqrt(dble(1 + 8 * size(DMAT)))) - 1) / 2
          call triprt('D_alpha in neci2molcas',' ',D_alpha,norb)
          call triprt('D_beta  in neci2molcas',' ',D_beta ,norb)
          call triprt('DMAT in neci2molcas',' ',DMAT,norb)
          call triprt('DSPN in neci2molcas',' ',DSPN,norb)
      END IF
      Return
      end subroutine read_single_neci_RDM


      subroutine bcast_2RDM(InFile)
        use filesystem, only : symlink_, strerror_, get_errno_
        character(len=*), intent(in) :: InFile
        character(len=1024) :: master
        integer :: lmaster1, err

        call prgmtranslate_master(InFile, master, lmaster1)
        call symlink_(trim(master), trim(InFile), err)
        if (err == 0) write(u6, *) strerror_(get_errno_())
      end subroutine bcast_2RDM


      ! Add your deallocations here. Called when exiting rasscf.
      subroutine cleanup()
        continue
      end subroutine


      end module fciqmc_read_RDM
