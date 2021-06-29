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
*----------------------------------------------------------------------*
* This subroutine reads from the formatted file DIFFPR coming from     *
* generated by MpProp the Slater Exponents, Factors and Nuclear Charges*
*----------------------------------------------------------------------*
      Subroutine Get_Slater(SlExpQ,LMltSlQ,outxyz,nAt)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "warnings.h"

      Dimension CoordTest(3),SlFactQ(6)
      Dimension SlExpQ(MxMltp+1,MxQCen),outxyz(MxQCen,3)
      Logical Exist,lCheck


*Open the file
      Lu=40
      Lu=IsFreeUnit(40)
      Call Opnfl('DIFFPR',Lu,Exist)
      If(.not.Exist) then
        Write(6,*)
        Write(6,*)' Can not locate output file DiffPr. '
        Call Quit(_RC_IO_ERROR_READ_)
      Endif
      Rewind(Lu)

*-- Read Number of centers and angular momentums.
*
      Read(Lu,101)nSlCentQ
      Read(Lu,101)LMltSlQ

* A first test
        nTestjhr=nAt*(nAt+1)/2
        If(nSlCentQ.ne.nTestjhr) then
        Write(6,*)'ERROR! Number of centers in DiffPr file',nSlCentQ
     &  ,' is different from number of centers obtained from RUNFILE'
     &  ,nTestjhr,' Check your files.'
        Call Quit(_RC_GENERAL_ERROR_)
        Endif

*-- Read Exponentials for the Centers
      Do iC=1,nSlCentQ
        lCheck=.false.
        Read(Lu,103)(CoordTest(k),k=1,3)
        ind=0
        Do jhr=1,nSlCentQ
          If(abs(CoordTest(1)-outxyz(jhr,1)).lt.1.0d-4) then
            If(abs(CoordTest(2)-outxyz(jhr,2)).lt.1.0d-4) then
              If(abs(CoordTest(3)-outxyz(jhr,3)).lt.1.0d-4) then
                lCheck=.true.
                ind=jhr
              Endif
            Endif
          Endif
        Enddo
        If(.not.lCheck) then
          write(6,*)'ERROR. Something is very wrong, coordinates'
     &//' of DiffPr and MpProp files do not match.'
     &//' DiffPr center',iC
        Endif
        Do l=0,LMltSlQ
          nS=l*(l+1)*(l+2)/6
          nT=(l+1)*(l+2)*(l+3)/6
          Read(Lu,104)SlExpQ(l+1,ind)
          Read(Lu,105)(SlFactQ(kk),kk=nS+1,nT)
*         Read(Lu,105)(SlFactQ(kk,ind),kk=nS+1,nT)
        End do
*Jose. No read nuclear charge
*       Read(Lu,104)PointP(ind)
        Read(Lu,*)
      End do

      Close(Lu)

101   Format(I5)
103   Format(3(F20.14))
104   Format(F20.14)
105   Format(3(F20.14))

      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_real_array(SlFactQ)
#endif
      End
