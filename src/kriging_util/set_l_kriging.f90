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
SUBROUTINE set_l_kriging(lv,nInter_In)
  use kriging_mod
  integer nInter_In
  real*8 lv(nInter_In)
!
! Set the characteristic length of all the components of the coordintes.
!
  If (nInter_In.eq.nInter) Then
    l(:)=lv(:)
  Else If (nInter.eq.1) Then
    l(:)=lv(1)
  Else
    Write (6,*) "setlkriging: illegal nInter value."
    Call Abend()
  End If
!
! Generate the covariance matrix
!
  call covarMatrix()
!
! Form the inverse of the covariance matrix times the generalized value vector.
!
  call kriging_model()

END SUBROUTINE set_l_kriging
