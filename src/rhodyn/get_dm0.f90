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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine get_dm0()
!***********************************************************************
!
!  Purpose : to get the initial density matrix in CSF basis
!            which used in propagation!
!
!***********************************************************************
  use rhodyn_data
  use rhodyn_utils, only: transform, dashes
  use definitions, only: wp, iwp, u6
  use stdalloc, only: mma_allocate, mma_deallocate
  use mh5, only: mh5_put_dset
  implicit none

  complex(kind=wp),dimension(:,:),allocatable:: DET2CSF, DM0_bas
  complex(kind=wp) :: Z
  complex(kind=wp), dimension(:), allocatable:: temp_dm
  integer(kind=iwp) :: lu !temporary io unit
  integer(kind=iwp) :: i, j, k, ii, jj, kk
  integer(kind=iwp), external :: isFreeUnit

  if (ipglob>3) write(u6,*) 'Begin get_dm0'

  if (p_style=='SO_THERMAL'.and.T==0) then
    p_style='SO'
  endif
  if (p_style=='SF_THERMAL'.and.T==0) then
    p_style='SF'
  endif

  if (runmode/=4) then
    call mma_allocate(DM0,nconftot, nconftot)
  else ! in charge migration case prepare DM in SF/SO basis
    call mma_allocate(DM0,lrootstot,lrootstot)
  endif
  DM0=zero

! If populate DET, read the transform matrix from DET to CSF to DET2CSF
  if (p_style=='DET') then
    call mma_allocate(DET2CSF,ndet_tot,nconftot)
    DET2CSF=zero
    ii=0
    jj=0
    kk=0
    do i=1,N
      if (i==1) then
        ii=0
      else
        ii=nconf(i-1)
      endif
      do j=1,ndet(i)
        jj=jj+1
        do k=1,nconf(i)
          kk=ii+kk
          DET2CSF(jj,kk)=DTOC(j,k,i)
        enddo
      enddo
    enddo
  endif

! Construct the initial density matrix for populating

  if (flag_so) then
    select case (p_style)
    case ('CSF')
      if ((N_Populated>nconftot).or.(N_populated<=0)) then
        call dashes()
        write(u6,*)'WARNING!!! Nr of the populated states', &
                N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      endif
      if (runmode/=4) then
        DM0(N_Populated,N_Populated)=one
      else
        call mma_allocate(DM0_bas,nconftot,nconftot)
        DM0_bas=zero
        DM0_bas(N_Populated,N_Populated)=one
        write(u6,*) 'DM element ',N_populated,' set to 1'
        call transform(DM0_bas,CSF2SO,DM0)
      endif

    case ('DET')
      call mma_allocate(DM0_bas,ndet_tot,ndet_tot)
      DM0_bas=zero
      if ((N_Populated>ndet_tot).or.(N_populated<=0)) then
        call dashes()
        write(u6,*)'WARNING!!! Nr of the populated states', &
                N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      endif
      DM0_bas(N_Populated,N_Populated)=one
      call transform(DM0_bas,DET2CSF,DM0)
      call mma_deallocate(DET2CSF)

    case ('SF')
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas=zero
      if ((N_Populated>lrootstot).or.(N_populated<=0)) then
        call dashes()
        write (u6,*) 'WARNING!!! Nr of the populated is wrong'
        call dashes()
        call abend()
      endif
      DM0_bas(N_Populated,N_Populated)=one
      call transform(DM0_bas,dcmplx(U_CI,0d0),DM0,.False.)

    case ('SO')
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas=zero
      if ((N_Populated>lrootstot).or.(N_populated<=0)) then
        call dashes()
        write(u6,*)'WARNING!!! Nr of the populated states', &
                  N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      endif
      DM0_bas(N_Populated,N_Populated)=one
      call transform(DM0_bas,CSF2SO,DM0,.False.)

    case ('SO_THERMAL')
      if (ipglob>3) then
        call dashes()
        write(u6,*) 'Printout the SO eigenvalues'
        write(u6,sint)'Number of the SO states:', lrootstot
        call dashes()
        do i=1,10
          write(u6,*) E_SO(i)
        enddo
      endif
      Z=zero
      do i=1,lrootstot
        Z=Z+exp(-(E_SO(i)-E_SO(1))/(k_B*T))
      enddo
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas=zero
      if (ipglob>3) then
        call dashes()
        write(u6,*) 'I            E_SO(I)-E_SO(1)               '// &
                   'exp(-(E_SO(I)-E_SO(1))/(k_B*T)),   '// &
                   'exp(-(E_SO(I)-E_SO(1))/(k_B*T))/Z'
        call dashes()
        do i=1,20
          write(u6,*) i,E_SO(i)-E_SO(1), &
                   exp(-(E_SO(i)-E_SO(1))/(k_B*T)), &
                   exp(-(E_SO(i)-E_SO(1))/(k_B*T))/Z
        enddo
      endif
      do i=1,lrootstot
        DM0_bas(i,i)=exp(-(E_SO(i)-E_SO(1))/(k_B*T))/Z
      enddo
      if (runmode/=4) then
        call transform(DM0_bas,CSF2SO,DM0,.False.)
      else
        DM0(:,:) = DM0_bas
      endif
    case default
      write(6,*)'Population style ',p_style, ' is not recognized'
      call abend()
    end select

  else
!     SO is off
!     has not checked yet
    select case (p_style)
    case ('CSF')
      if ((N_Populated>nconftot).or.(N_populated<=0)) then
        call dashes()
        write(u6,*)'WARNING!!! Nr of the populated states', &
                 N_populated,'read from the input file is wrong'
        call dashes()
        call abend()
      endif
      if (runmode/=4) then
        DM0(N_Populated,N_Populated)=one
      else
        call mma_allocate(DM0_bas,nconftot,nconftot)
        DM0_bas=zero
        DM0_bas(N_Populated,N_Populated)=one
        write(u6,*) 'DM element ',N_populated,' set to 1'
        ! transform right to SF basis for propagation
        call transform(DM0_bas,dcmplx(U_CI,0d0),DM0)
      endif

    case ('DET')
      call mma_allocate(DM0_bas,ndet_tot,ndet_tot)
      DM0_bas=zero
      if ((N_Populated>NDET_TOT).or.(N_populated<=0)) then
        call dashes()
        write(u6,*)'WARNING!!! Nr of the populated states', &
                  N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      endif
      DM0_bas(N_Populated,N_Populated)=one
      call transform(DM0_bas,DET2CSF,DM0)
      call mma_deallocate(DET2CSF)

    case ('SF')
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas=zero
      if ((N_Populated>lrootstot).or.(N_populated<=0)) then
        call dashes()
        write(u6,*)'WARNING!!! Nr of the populated states', &
                N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      endif
      DM0_bas(N_Populated,N_Populated)=one
      call transform(DM0_bas,dcmplx(U_CI,0d0),DM0,.False.)

    case ('SO')
      call dashes()
      write(u6,*) 'WARNING! IF IFSO=OFF, CAN NOT POPULATE SO STATE'
      call dashes()
      call abend()

    case ('SF_THERMAL')
!       Construct the total eigenvalue of SF states to E_SF
      ii=0
      do i=1,N
        do j=1,lroots(i)
          ii=ii+1
          E_SF(ii)=E(i,j)
        enddo
      enddo
      Z=0d0
      do i=1,nconftot
        Z=Z+exp(-(E_SF(i)-E_SF(1))/(k_B*T))
      enddo
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas=zero
      call dashes()
      write(u6,*)' I               '// &
                 '(E_SF(I)-E_SF(1))/(k_B*T)              '// &
                 'exp(-(E_SF(I)-E_SF(1))/(k_B*T))  '// &
                 'exp(-(E_SF(I)-E_SF(1))/(k_B*T))/Z'
      call dashes()
      do i=1,nconftot
        write(u6,*)i,(E_SF(i)-E_SF(1))/(k_B*T), &
                   exp(-(E_SF(i)-E_SF(1))/(k_B*T)), &
                   exp(-(E_SF(i)-E_SF(1))/(k_B*T))/Z
        DM0_bas(i,i)=exp(-(E_SF(i)-E_SF(1))/(k_B*T))/Z
      enddo
      call transform(DM0_bas,dcmplx(U_CI,0d0),DM0,.False.)
    case default
      write(6,*)'Population style ',p_style, ' is not recognized'
      call abend()
    end select
  endif !ifso

! it will be transformed as usual further
  if (p_style=='FROMFILE') then
   call mma_allocate(temp_dm, NState)
   lu = isFreeUnit(11)
   call molcas_open(lu,'INDENS')
   do i=1,NState
     do j=1,Nstate
       read(lu,'(E16.8)',advance='no') temp_dm(j)
     enddo
     DM0(i,:) = dcmplx(temp_dm)
     read(lu,*)
   enddo
   close(lu) ! close INDENS file
   call mma_deallocate(temp_dm)
 endif

  if (allocated(DM0_bas)) call mma_deallocate(DM0_bas)

! Printout the initial density matrix to SDPREP in CSF basis
  ! not CM case
  if (runmode/=4) then
    if (ipglob>3) then
      call dashes()
      write(u6,*) 'The initial density matrix is saved in SDPREP'
      call dashes()
    endif
    call mh5_put_dset(prep_dm_r, dble(DM0))
    call mh5_put_dset(prep_dm_i, aimag(DM0))
  endif

  if (ipglob>3) write(u6,*) 'End get_dm0'

end
