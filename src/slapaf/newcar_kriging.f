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
* Copyright (C) 2019, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine NewCar_Kriging(kIter,nAtom,nInter,SaveBMx,RefIter,
     &                          Error)
      use Slapaf_Info, only: Cx, BMx
      use Slapaf_Parameters, only: PrQ
      Implicit None
#include "info_slapaf.fh"
#include "db.fh"
#include "stdalloc.fh"
      Integer :: kIter,nAtom,nInter,RefIter

      Real*8, Allocatable :: Coor(:,:), BMx_Tmp(:,:)
      Logical :: Numerical,PrQ_Save,Error,SaveBMx
*
      Call mma_allocate(Coor,3,nAtom,Label='Coor')
      Coor(:,:) = Cx(:,:,kIter)

      Call mma_allocate(BMx_tmp,SIZE(BMx,1),SIZE(BMx,2),Label='BMx_tmp')
      BMx_tmp(:,:) = BMx(:,:)
*
      Numerical=.False.
      PrQ_Save=PrQ
      PrQ=.False.
*
      Force_dB=SaveBMx
*
      Call NewCar(kIter,nAtom,nInter,Coor,iSym,mTtAtm,
     &            Numerical,RefIter,Error)
*
      PrQ=PrQ_Save
      Force_dB=.False.
      Call mma_deallocate(Coor)

      If (.NOT.SaveBMx) Then
         Call mma_deallocate(BMx)
         Call mma_allocate(BMx,SIZE(BMx_tmp,1),SIZE(BMx_Tmp,2),
     &                     Label='BMx')
         BMx(:,:) = BMx_tmp(:,:)
      End If

      Call mma_deallocate(BMx_tmp)
*
      End Subroutine NewCar_Kriging
