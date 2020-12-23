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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
Subroutine Hessian_Kriging(x0_,ddy_,ndimx)
  use kriging_mod
  Implicit None
  Integer ndimx
  Real*8 x0_(ndimx),ddy_(ndimx,ndimx)
!
!#define _Hess_Test
#ifdef _Hess_Test
  Real*8 Scale,Delta,Fact,tgrad(ndimx),thgrad(ndimx)
  Real*8 HessT, tmp
  Integer i, j
  HessT = 1.0D-3
#endif
!
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
  x0(:) = x0_(:)
!
  call covarvector(2) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(2)
  ddy_(:,:) = hpred(:,:)
!
#ifdef _Hess_Test
! Numerical Hessian of GEK
  write(6,*) 'Begining Numerical Hessian'
!
  hpred(:,:) = 0.0D0
  Scale=0.01D0
  write(6,*) 'Hess Threshold',HessT
!
  do i = 1,nInter
    tmp=x0(i)
!
    Delta = 1.0D-5!Max(Abs(x_(i,1)),1.0D-5)*Scale
!
    x0(i) = tmp + Delta
    call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
    call predict(1)
    tgrad=gpred(:)
!
    x0(i) = tmp - Delta
    call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
    call predict(1)
    thgrad=gpred(:)
!
    do j=1,nInter
      hpred(i,j) = (tgrad(j)-thgrad(j))/(2.0D0*Delta)
    enddo
    x0(i) = tmp
  enddo
! Comparing Analytical solution with Numerical
  do i = 1,nInter
    do j = 1,nInter
      write(6,*) 'i,j',i,j
      write(6,*) 'hpred, ddy_',hpred(i,j),ddy_(i,j)
      if (abs(ddy_(i,j)-hpred(i,j)).gt.HessT) then
        Write(6,*) 'Error in entry',i,',',j,'of the hessian matrix'
        Call RecPrt('Anal Hess',' ',ddy_,nInter,nInter)
        Call RecPrt('Num Hess',' ',hpred,nInter,nInter)
        Write(6,*) 'abs(ddy_(i,j)+ HessT)',abs(ddy_(i,j)+ HessT)
        Write(6,*) 'abs(ddy_(i,j)- HessT)',abs(ddy_(i,j)- HessT)
        Write(6,*) 'abs(hpred(i,j))',abs(hpred(i,j))
        Call Abend()
      endif
    enddo
  enddo
#endif
!
  return
End Subroutine Hessian_Kriging
