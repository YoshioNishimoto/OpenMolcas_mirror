************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NIntX(AOInt,mGrid,TabAO1,TabAO2,nBfn,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_info.fh"
      Real*8 AOInt(nBfn,nBfn,nD), TabAO1(nFn,mGrid,nBfn,nD),
     &       TabAO2(mAO,mGrid,nBfn)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Functional_type.eq.LDA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1 = ToAdd1
     &                + TabAO1(1,iGrid,iCB,iD)*TabAO2(1,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
     &               + TabAO1(5,iGrid,iCB,iD) *(TabAO2(5,iGrid,jCB)
     &                                         +TabAO2(8,iGrid,jCB)
     &                                        +TabAO2(10,iGrid,jCB))
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Write (6,*) 'DFT_Int: Illegal functional type!'
         Call Abend()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End
