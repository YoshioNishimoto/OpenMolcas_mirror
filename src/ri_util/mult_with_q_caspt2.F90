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
! Copyright (C) Jonas Bostrom                                          *
!***********************************************************************

subroutine Mult_with_Q_CASPT2(nBas_aux,nBas,nIrrep,SubAux)
!***********************************************************************
! Author: Jonas Bostrom
!
! Purpose: Multiply MP2 A~_sep and B~_sep with inverse cholesky factors.
!
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Cholesky, only: nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nIrrep, nBas_Aux(1:nIrrep), nBas(1:nIrrep)
logical(kind=iwp), intent(in) :: SubAux
integer(kind=iwp) :: i, iAdrQ, id, iOffQ1, iOpt, iost, ip_B, ip_B2, iSym, j, jSym, jVec, kSym, kVec, l_A, l_A_ht, l_A_t, l_B_t, &
                     l_Q, lRealName, Lu_Q, LUGAMMA, LuGamma2, LUAPT2, lVec, MaxMem, nBas2, nBasTri, nLR, nLRb(8), nseq, NumAux, &
                     NumCV, NumVecJ, NumVecK, nVec
real(kind=wp) :: aaa, Fac, TotCPU0, TotCPU1, TotWall0, TotWall1
logical(kind=iwp) :: is_error
character(len=4096) :: RealName
character(len=6) :: Name_Q
real(kind=wp), allocatable :: A(:), A_ht(:), A_t(:), B_t(:), QVec(:)
character(len=*), parameter :: SECNAM = 'Mult_with_Q_MP2'
integer(kind=iwp), external :: IsFreeUnit
return

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TotCPU0,TotWall0)
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym
  NumCV = NumCho(iSym)
  NumAux = nBas_Aux(iSym)-1
  if (SubAux) NumAux = nBas_Aux(iSym)-1
  nLR = 0
  do jSym=1,nSym
    kSym = Mul(iSym,jSym)
    nLR = nLR+nBas(jSym)*nBas(kSym)
  end do
  nLRb(iSym) = nLR
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym

  nBas2 = nLRb(iSym)
  !write(u6,*) 'nBas2 = ',nBas2
  NumCV = NumCho(iSym)
  NumAux = nBas_Aux(iSym)-1
  if (SubAux) NumAux = nBas_Aux(iSym)-1

  ! Get Q-vectors from disk
  ! -----------------------

  l_Q = NumCV*NumAux
  call mma_allocate(QVec,l_Q,Label='Q_Vector')

  Lu_Q = IsFreeUnit(7)
  write(Name_Q,'(A4,I2.2)') 'QVEC',iSym-1
  call DaName_MF_WA(Lu_Q,Name_Q)

  iOpt = 2
  iAdrQ = 0
  call dDaFile(Lu_Q,iOpt,QVec,l_Q,iAdrQ)

  ! Get MP2 A-tilde vectors from disk
  ! ---------------------------------

  l_A_t = NumCV*NumCV
  l_A = NumAux*NumAux
  l_A_ht = NumAux*NumCV
  call mma_allocate(A_t,l_A_t,Label='A_t')
  call mma_allocate(A,l_A,Label='A')
  call mma_allocate(A_ht,l_A_ht,Label='A_ht')

  !LUCMOPT2 = 61
  !call PrgmTranslate('CMOPT2',RealName,lRealName)
  !call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.false.,1,'OLD',is_error)

  !read(LuCMOPT2,iostat=iost) A_t
  !if (iost < 0) then
  !  write (6,*) 'Maybe, you did not add GRAD or GRDT keyword in &CASPT2?'
  !  write (6,*) 'Please add either one, if this is single-point gradient calculation.'
  !  write (6,*) 'Otherwise, something is wrong...'
  !  call abend()
  !end if

  ! Read A_PT2 from LUAPT2
  LuAPT2 = isFreeUnit(68)
  call daname_mf_wa(LUAPT2,'A_PT2')
  id = 0
  call ddafile(LUAPT2,2,A_t,l_A_t,id)

  ! Symmetrized A_PT2
  do i=1,NumCV
    do j=1,i
      aaa = Half*(A_t((j-1)*NumCV+i)+A_t((i-1)*NumCV+j))
      A_t((j-1)*NumCV+i) = aaa
      A_t((i-1)*NumCV+j) = aaa
    end do
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'Q-vectors'
  do i=1,l_Q
    write(u6,*) QVec(i)
  end do

  write(u6,*) 'A-vectors'
  do i=1,l_A_t
    write(u6,*) A_t(i)
  end do
# endif

  ! Make first halftransformation to cholesky-base
  ! ----------------------------------------------

  call dGemm_('N','N',NumAux,NumCV,NumCV,One,QVec,NumAux,A_t,NumCV,Zero,A_ht,NumAux)

  call dGemm_('N','T',NumAux,NumAux,NumCV,One,A_ht,NumAux,QVec,NumAux,Zero,A,NumAux)

  ! Put transformed A-vectors back on disk

  !rewind(LuCMOPT2)
  !write(LuCMOPT2) A
  !close(LuCMOPT2)

  ! write A_PT2 to LUAPT2
  id = 0
  call ddafile(LUAPT2,1,A,l_A,id)

  call mma_deallocate(A_t)
  call mma_deallocate(A)
  call mma_deallocate(A_ht)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_maxDBLE(MaxMem)
  MaxMem = 9*(MaxMem/10)
  call mma_allocate(B_t,MaxMem,Label='B_t')

  ! Applicable to C1 only
  nBasTri = nTri_Elem(nBas(1))
  nVec = MaxMem/(2*nBasTri+1) ! MaxMem/(2*nLRb(iSym)+1)
  nVec = min(max(NumCV,NumAux),nVec)
  if (nVec < 1) call SysAbendMsg(SecNam,'nVec is non-positive','[1]')

  l_B_t = nBasTri*nVec ! nLRb(iSym)*nVec
  ip_B = 1+l_B_t
  ip_B2 = ip_B+l_B_t

  LuGAMMA = IsFreeUnit(65)
  call PrgmTranslate('GAMMA',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nBas2*8,'OLD',is_error)

  LuGamma2 = IsFreeUnit(67)
  call PrgmTranslate('GAMMA2',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuGamma2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,NumAux*8,'REPLACE',is_error)

  ! The B-vectors should be read one batch at the time
  ! --------------------------------------------------

  do kVec=1,NumAux,nVec
    NumVecK = min(nVec,NumAux-(kVec-1))

    do jVec=1,NumCV,nVec
      NumVecJ = min(nVec,NumCV-(jVec-1))

      do lVec=1,NumVecJ
        read(Unit=LuGAMMA,Rec=jVec+lVec-1) B_t(ip_B2:ip_B2+nBas2-1)
        ! symmetrize (mu nu | P)
        ! only the lower triangle part is used
        nseq = 0
        do i=1,nBas(1)
          do j=1,i
            nseq = nseq+1
            aaa = B_t(ip_B2+i-1+nBas(1)*(j-1))+B_t(ip_B2+j-1+nBas(1)*(i-1))
            B_t(nseq+nBasTri*(lVec-1)) = aaa
          end do
        end do
      end do

      Fac = Zero
      if (jVec /= 1) Fac = One
      iOffQ1 = kVec+NumAux*(jVec-1)
      call dGemm_('N','T',nBasTri,NumVecK,NumVecJ,One,B_t,nBasTri,QVec(iOffQ1),NumAux,Fac,B_t(ip_B),nBasTri)
    end do

    ! Because scaling with 0.5 is omitted when symmetrized
    call DScal_(nBasTri*NumVecK,Half,B_t(ip_B),1)

    ! (mu nu | P) --> (P | mu nu)
    if (max(NumCV,NumAux) == nVec) then
      nseq = 0
      do lVec=1,nBas(1)
        do jVec=1,lVec
          nseq = nseq+1
          call DCopy_(NumAux,B_t(ip_B+nseq-1),nBasTri,B_t,1)
          write(unit=LuGAMMA2,rec=nseq) B_t(1:NumAux)
        end do
      end do
    else
      nseq = 0
      do lVec=1,nBas(1)
        do jVec=1,lVec
          nseq = nseq+1
          if (kVec /= 1) read(unit=LuGAMMA2,rec=nseq) B_t(1:kVec-1)
          call DCopy_(NumVecK,B_t(ip_B+nseq-1),nBasTri,B_t(kVec),1)
          write(unit=LuGAMMA2,rec=nseq) B_t(1:kVec+NumVecK-1)
        end do
      end do
    end if
  end do

  close(LuGAMMA2)

  call mma_deallocate(B_t)
  call mma_deallocate(QVec)

  call DaClos(Lu_Q)
  call DaClos(LUAPT2)
  close(LuGAMMA)

end do ! iSym

call CWTime(TotCPU1,TotWall1)
!write(u6,*) 'CPU/Wall Time for mult_with_q_caspt2:',totcpu1-totcpu0,totwall1-totwall0

return

end subroutine Mult_with_Q_CASPT2

subroutine Mult_with_Q_CASPT22(nAct)
!***********************************************************************
! Author: Jonas Bostrom
!
! Purpose: Multiply MP2 A~_sep and B~_sep with inverse cholesky factors.
!
!***********************************************************************

use Basis_Info, only: nBas, nBas_Aux
use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, Two
use Definitions, only: wp, iwp
use RI_glob, only: LuCVector, nIJR

use pso_stuff, only: V_k, D0, AOrb, nnP, nV_K, Txy, n_Txy, G2

implicit none
!#include "cholesky.fh"
!#include "etwas.fh"
integer(kind=iwp) :: iSym, jSym, kSym, lSym, i, j, k, l, nseq, NumAux, LuGamma, lRealName, iost, is_error, &
                     iVec, iMO1, iMO2, iVec_, kAct, lAct, iSO, iMOleft, iMOright, lCVec
character(len=4096) :: RealName
real(kind=wp) :: aaa, fact, TotCPU0, TotCPU1, TotWall0, TotWall1, tmp
real(kind=wp), allocatable :: B_t(:), RDMWRK(:,:,:,:), g2sq(:,:,:,:), &
  trfwrk(:,:,:,:), trfwrk2(:,:,:,:)

integer(kind=iwp), intent(in) :: nAct
integer(kind=iwp), external :: IsFreeUnit
integer(kind=iwp) :: nBasSq,it,iu,iv,ix,nseq1,nseq2,nseq3,nseq4,nseq5,nseq6,itu,ivx,ituvx,iact,jact,nsym,itx, &
                     ivu,itv,ixu

!real(kind=wp), pointer :: Xki(:), Xli(:)
!type(V1) :: Xki2(2), Xki3(2), Xli2(2), Xli3(2)

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TotCPU0,TotWall0)
!                                                                      *
!***********************************************************************
!                                                                      *
NumAux = nBas_Aux(0)-1
nBasSq = nBas(0)**2
LuGAMMA = IsFreeUnit(65)
call PrgmTranslate('GAMMA',RealName,lRealName)
call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nBasSq*8,'OLD',is_error)
!write (*,*) "numbas = ", nbas(0)
!write (*,*) "numaux = ", nbas_aux(0)-1

call mma_allocate(RDMWRK,nBas(0),nBas(0),nBas(0),nBas(0),Label='RDMWRK')
call mma_allocate(B_t,nBasSq,Label='B_t')

if (is_error) then
  write (6,*) "reading RDMWRK failed"
  RDMWRK = 0.0d+00
else
  do i = 1, nBasSq
    Read (LuGamma,Rec=i) (B_t(k),k=1,nBasSQ)
    call dcopy_(nBasSq,B_t,1,RDMWRK(1,1,i,1),1)
  end do

  !! for the moment, work with the non-perturbed density only
  !! call dcopy_(nbassq*nbassq,[0.0d+00],0,RDMWRK,1)

  close (LuGamma)
  ! write (6,*) ":: two-electron initial"
  ! do i = 1, nBas(0)
  !   do j = 1, nBas(0)
  !     do k = 1, nBas(0)
  !       do l = 1, nBas(0)
  ! write (6,'(4i3,f20.10)') i,j,k,l,RDMWRK(i,j,k,l)
  !       end do
  !     end do
  !   end do
  ! end do
end if

call mma_deallocate(B_t)

!lda = size(CMOi(1)%SB(1)%A2,1)
!ik = 1+lda*(kSO-1)
!il = 1+lda*(lSO-1)
!Xki(1:) => CMOi(1)%SB(1)%A1(ik:)
!Xli(1:) => CMOi(1)%SB(1)%A1(il:)

!write (*,*) "D0(:,1)" !! inactive density
!call prtril(d0(1,1),nbas(0))
!write (*,*) "D0(:,2)" !! (total - inactive/2) = (active + inactive/2)  density
!call prtril(d0(1,2),nbas(0))
!write (*,*) "D0(:,3)" !! averaged active density
!call prtril(d0(1,3),nbas(0))
!write (*,*) "D0(:,4)"
!call prtril(d0(1,4),nbas(0))
!write (*,*) "D0(:,5)"
!call prtril(d0(1,5),nbas(0))

nSym = 1

do iSym=0,nSym-1

  do k = 1, nBas(iSym)
  do l = 1, nBas(iSym)
    nseq2 = max(k,l)*(max(k,l)-1)/2 + min(k,l)

    do i = 1, nBas(iSym)
    do j = 1, nBas(iSym)
      nseq1 = max(i,j)*(max(i,j)-1)/2 + min(i,j)

      nseq3 = max(i,k)*(max(i,k)-1)/2 + min(i,k)
      nseq4 = max(j,l)*(max(j,l)-1)/2 + min(j,l)
      nseq5 = max(i,l)*(max(i,l)-1)/2 + min(i,l)
      nseq6 = max(j,k)*(max(j,k)-1)/2 + min(j,k)

      !! Coulomb and exchange
      aaa = D0(nseq1,1)*D0(nseq2,2) &
          + D0(nseq1,2)*D0(nseq2,1) &
          - D0(nseq3,1)*D0(nseq4,2)*0.25d+00 &
          - D0(nseq3,2)*D0(nseq4,1)*0.25d+00 &
          - D0(nseq5,1)*D0(nseq6,2)*0.25d+00 &
          - D0(nseq5,2)*D0(nseq6,1)*0.25d+00 &
          + D0(nseq1,3)*D0(nseq2,4) &
          + D0(nseq1,4)*D0(nseq2,3) &
          - D0(nseq3,3)*D0(nseq4,4)*0.25d+00 &
          - D0(nseq3,4)*D0(nseq4,3)*0.25d+00 &
          - D0(nseq5,3)*D0(nseq6,4)*0.25d+00 &
          - D0(nseq5,4)*D0(nseq6,3)*0.25d+00 &
          + D0(nseq1,1)*D0(nseq2,5) &
          + D0(nseq1,5)*D0(nseq2,1) &
          - D0(nseq3,1)*D0(nseq4,5)*0.25d+00 &
          - D0(nseq3,5)*D0(nseq4,1)*0.25d+00 &
          - D0(nseq5,1)*D0(nseq6,5)*0.25d+00 &
          - D0(nseq5,5)*D0(nseq6,1)*0.25d+00
      RDMWRK(i,j,k,l) = RDMWRK(i,j,k,l) + aaa*0.5d+00

    end do ! j
    end do ! i
  end do ! l
  end do ! k
end do ! iSym

call mma_allocate(G2SQ,nAct,nAct,nAct,nAct,Label="G2SQ")

!write (*,*) "MO coefficients(1)"
!do i = 1, nact
!do j = 1, nbas(0)
!write (*,'(2i3,f20.10)') i,j,AOrb(1)%SB(1)%A2(i,j)
!end do
!end do
!write (*,*) "MO coefficients(2)"
!do i = 1, nact
!do j = 1, nbas(0)
!write (*,'(2i3,f20.10)') i,j,AOrb(2)%SB(1)%A2(i,j)
!end do
!end do

do ivec = 1, 2
  do it = 1, nact
  do iu = 1, nact
  itu = iTri(it,iu)
  do iv = 1, nact
  do ix = 1, nact
  ivx = iTri(ix,iv)
    ituvx = iTri(ivx,itu)
    aaa = g2(ituvx,ivec)
    if (itu.ge.ivx .and. iv.eq.ix) aaa = aaa+aaa
    if (itu.lt.ivx .and. it.eq.iu) aaa = aaa+aaa
    g2sq(it,iu,iv,ix) = aaa*0.5d+00
!   write (*,'(4i3,f20.10)') it,iu,iv,ix, g2sq(it,iu,iv,ix)
  end do
  end do
  end do
  end do

  call mma_allocate(trfwrk ,nact,nact,nact   ,nbas(0),label='trfwrk')
  call mma_allocate(trfwrk2,nact,nact,nbas(0),nbas(0),label='trfwrk2')
  call dgemm_('N','N',nact**3,nbas(0),nact, &
              1.0d+00,g2sq,nact**3,aorb(ivec)%sb(1)%a2(1,1),nact, &
              0.0d+00,trfwrk,nact**3)
  do i = 1, nbas(0)
    call dgemm_('N','N',nact**2,nBas(0),nact, &
              1.0D+00,trfwrk(1,1,1,i),nact**2,aorb(ivec)%sb(1)%a2(1,1),nact, &
              0.0D+00,trfwrk2(1,1,1,i),nact**2)
  end do
  call mma_deallocate(trfwrk)

  call mma_allocate(B_t,nact*nbas(0),label="B_t")

  do j = 1, nbas(0)
  do i = 1, nbas(0)
    call dgemm_('N','N',nact,nbas(0),nact, &
                1.0d+00,trfwrk2(1,1,i,j),nact,aorb(ivec)%sb(1)%a2(1,1),nact, &
                0.0d+00,b_t,nact)
    call dgemm_('T','N',nbas(0),nbas(0),nact, &
                1.0d+00,aorb(ivec)%sb(1)%a2(1,1),nact,b_t,nact, &
                1.0d+00,RDMWRK(1,1,i,j),nbas(0))
  end do
  end do
  call mma_deallocate(trfwrk2)
  call mma_deallocate(B_t)
end do

call mma_deallocate(G2SQ)

write (6,*) ": one-electron"
do i = 1, nbas(0)
do j = 1, nbas(0)
nseq1 = max(i,j)*(max(i,j)-1)/2 + min(i,j)
write (6,'(2i3,f20.10)') i,j,d0(nseq1,1)*0.5d+00 + d0(nseq1,2)
end do
end do

do i = 1, nbas(0)
do j = 1, nbas(0)
  do k= 1, nbas(0)
  do l = 1, nbas(0)
    aaa = rdmwrk(i,j,k,l) + rdmwrk(i,j,l,k) + rdmwrk(j,i,k,l) + rdmwrk(j,i,l,k) &
        + rdmwrk(k,l,i,j) + rdmwrk(l,k,i,j) + rdmwrk(k,l,j,i) + rdmwrk(l,k,j,i)
    aaa = aaa/8.0d+00
    rdmwrk(i,j,k,l) = aaa
    rdmwrk(i,j,l,k) = aaa
    rdmwrk(j,i,k,l) = aaa
    rdmwrk(j,i,l,k) = aaa
    rdmwrk(k,l,i,j) = aaa
    rdmwrk(l,k,i,j) = aaa
    rdmwrk(k,l,j,i) = aaa
    rdmwrk(l,k,j,i) = aaa
  end do
  end do
end do
end do

write (6,*) ":: two-electron"
do i = 1, nBas(0)
  do j = 1, nBas(0)
    do k = 1, nBas(0)
      do l = 1, nBas(0)
write (6,'(4i3,f20.10)') i,j,k,l,RDMWRK(i,j,k,l)
      end do
    end do
  end do
end do

call mma_deallocate(RDMWRK)

call CWTime(TotCPU1,TotWall1)
!write(u6,*) 'CPU/Wall Time for mult_with_q_caspt2:',totcpu1-totcpu0,totwall1-totwall0

return

end subroutine Mult_with_Q_CASPT22
