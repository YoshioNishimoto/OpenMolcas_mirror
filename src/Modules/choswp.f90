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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
Module ChoSwp
Implicit none
Private
Public:: iQuAB, iQuAB_L, iQuAB_Hidden, iQuAB_L_Hidden, pTemp, iQuAB_here, &
         nnBstRSh, nnBstRSh_Hidden, nnBstRSh_G, nnBstRSh_L_Hidden, pTemp3, &
         iiBstRSh, iiBstRSh_Hidden, iiBstRSh_G, iiBstRSh_L_Hidden


Integer, Allocatable, Target:: iQuAB_Hidden(:,:), iQuAB_L_Hidden(:,:), iQuAB_here(:,:)
Integer, Pointer:: iQuAB(:,:)=>Null() , iQuAB_L(:,:)=>Null(), pTemp(:,:)=>Null()

Integer, Allocatable, Target:: nnBstRSh_Hidden(:,:,:), nnBstRSh_L_Hidden(:,:,:)
Integer, Allocatable, Target:: iiBstRSh_Hidden(:,:,:), iiBstRSh_L_Hidden(:,:,:)
Integer, Pointer:: nnBstRSh(:,:,:)=>Null(), nnBstRSh_G(:,:,:)=>Null(), pTemp3(:,:,:)=>Null()
Integer, Pointer:: iiBstRSh(:,:,:)=>Null(), iiBstRSh_G(:,:,:)=>Null()
End Module ChoSwp
