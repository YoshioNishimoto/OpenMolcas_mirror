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
* Copyright (C) Ben Swerts                                             *
*               2020, Roland Lindh                                    *
************************************************************************
      Subroutine FragExpand(nInfo,LuRd,DInf,nDInf)
************************************************************************
*                                                                      *
*    Objective: To expand the data for the fragments and append them   *
*               to regular arrays in info.fh                           *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : GetBS                                                   *
*                                                                      *
*     Author: Ben Swerts                                               *
*                                                                      *
************************************************************************
      Use Basis_Info
      Implicit None
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "print.fh"
      integer     nDInf, nInfo, storageSize, LineWords
      parameter(  storageSize = 200, LineWords=storageSize/8)
      Real*8      DInf(nDInf), eqBasis(LineWords)
      Integer     BasisTypes(4), nDel(MxAng),
     &            LenLbl, LuRd, iAtom, ib, iBas, iCnttp, iCntr,
     &            idummy, ii, Indx, iOptn, iSh, iShll, jShll,
     &            lAng, Last, LenBSL, lSTDINP, mCnttp, mdc, nAIMP, ndc,
     &            ninfo_stupid, nVal, nPrj, nSRO, nSOC, nPP, nProj,
     &            StayAlone,
     &            ipVal_, ipPrj_, ipSRO_, ipSOC_, ipPP_
      Real*8      x1, y1, z1
      Character*4  label
      Character*13 DefNm
      Character*80 Ref(2)
      Character*(storageSize) sBasis
      Equivalence( sBasis, eqBasis)
      Character *256 Basis_lib, Fname
      Character*180  STDINP(mxAtom*2)
!#define _DEBUG_
#ifdef _DEBUG_
      Integer i
#endif
* external functions and procedures
      Integer     iMostAbundantIsotope, iCLast
      Real*8      NucExp, rMass
      External    NucExp, rMass, iMostAbundantIsotope, iCLast
      Data DefNm/'basis_library'/
*                                                                      *
************************************************************************
*                                                                      *
#include "getbs_interface.fh"
*                                                                      *
************************************************************************
*                                                                      *

      UnNorm = .False.
      LenLbl=0
      mdc = mdciCnttp(nCnttp)+dbsc(nCnttp)%nCntr
      BasisTypes(1)=0
      BasisTypes(2)=0
      BasisTypes(3)=0
      BasisTypes(4)=0
      iShll = Mx_Shll-1
      lSTDINP=0
#ifdef _DEBUG_
      write(6,*) 'nCnttp, iShll, mdc = ',nCnttp,iShll,mdc
      write(6,*) 'ipExp(iShll+1) = ',ipExp(iShll+1)
#endif
*
      mCnttp = nCnttp
      ndc = 0

#ifdef _DEBUG_
       write(6,'(A,i6)')'FragExpand: just before the ''Do 1000 iCnttp'''
       write(6,'(A,i6)') 'FragExpand:       mdc          = ',mdc
       write(6,'(A,i6)') 'FragExpand:    mCnttp          = ',mCnttp
       write(6,'(A,10i6)') 'FragExpand: mdciCnttp(nCnttp)  = ',
     &                                    (mdciCnttp(i),i=1,mCnttp)
       write(6,'(A,10i6)') 'FragExpand:     dbsc(nCnttp)%nCntr  = ',
     &                              (dbsc(i)%nCntr,i=1,mCnttp)
       write(6,'(A,10I6)') 'FragExpand: nFragType(nCnttp)  = ',
     &                              (dbsc(i)%nFragType,i=1,mCnttp)
       write(6,'(A,10I6)') 'FragExpand: nFragCoor(nCnttp)  = ',
     &                              (dbsc(i)%nFragCoor,i=1,mCnttp)
#endif

      Do iCnttp = 1, mCnttp
        If(dbsc(iCnttp)%nFragType.le.0) Then
          ndc = ndc + dbsc(iCnttp)%nCntr
          Cycle
        End If
        Do iCntr = 1, dbsc(iCnttp)%nCntr
          ndc = ndc + 1
          Do iAtom = 1, dbsc(iCnttp)%nFragCoor
*
* create a new basis set center
*
            nCnttp = nCnttp + 1
            If (nCnttp.gt.Mxdbsc) Then
              Write (6,*) ' Increase Mxdbsc'
              Call ErrTra
              Call Quit_OnUserError()
            End If
* read the associated basis set in sBasis
            iBas = int(dbsc(iCnttp)%FragCoor(1,iAtom))
            call dcopy_(LineWords,
     &       dbsc(iCnttp)%FragType(1,iBas), 1, eqBasis, 1)
* get the basis set directory
            LenBSL = Len(sBasis)
            Last = iCLast(sBasis,LenBSL)
            Indx = Index(sBasis,'/')
            If(Indx.eq.0) Then
              Call WhichMolcas(Basis_lib)
              If(Basis_lib(1:1).ne.' ') then
                StayAlone = 1
                ib = index(Basis_lib,' ')-1
                If(ib.lt.1) Call SysAbendMsg('fragexpand',
     &                      'Too long PATH to MOLCAS',' ')
                Fname=Basis_lib(1:ib)//'/basis_library'
              Else
                Fname=DefNm
              Endif
              Indx = Last+1
              Bsl(nCnttp) = Trim(sBasis)
            Else
              Fname = sBasis(Indx+2:Last)
              If (Fname.eq.' ') Then
                Write (6,*) ' No basis set library specified for'
                Write (6,'(A,A)') 'Fname=',Fname
                Call Quit_OnUserError()
              End If
 1001         If (Fname(1:1).eq.' ') Then
                Fname(1:79)=Fname(2:80)
                Fname(80:80) = ' '
                Go To 1001
              End If
              Bsl(nCnttp)=sBasis(1:Indx-1)
            Endif
c           write(*,*) 'Setting Bsl(',nCnttp,') to ',Bsl(nCnttp)
c           write(*,*) 'Fname = ',Fname
* now Fname contains the basis set directory and Bsl(nCnttp) contains
* the basis set label
            jShll = iShll
            ExpNuc(nCnttp)=-One
            SODK(nCnttp)=.False.
            mdciCnttp(nCnttp)=mdc
            Call GetBS(Fname,sBasis(1:Indx-1),Indx-1,lAng,ipExp,
     &                 ipCff,ipCff_Cntrct,ipCff_Prim,ipFockOp,nExp,
     &                 nBasis,nBasis_Cntrct,MxShll,iShll,MxAng,
     &                 Charge(nCnttp),iAtmNr(nCnttp),BLine,Ref,
     &                 PAM2(nCnttp),FockOp(nCnttp),
     &                 ECP(nCnttp),NoPairL(nCnttp),SODK(nCnttp),ipBk,
     &                 CrRep(nCnttp),nProj,nAIMP,ipAkl,ip_Occ,iOptn,
     &                 UnNorm,nDel,
     &                 nVal,   nPrj,   nSRO,   nSOC,  nPP,
     &                 ipVal_, ipPrj_, ipSRO_, ipSOC_,ipPP_,
     &                 LuRd,BasisTypes,AuxCnttp(nCnttp),idummy,
     &                 STDINP,lSTDINP,.False.,.true.,' ',
     &                 DInf,nDInf,nCnttp)
            iAngMx=Max(iAngMx,lAng)
            Transf(jShll+1)=.False.
            Prjct(jShll+1)=.False.
            Transf(jShll+2)=.False.
            Prjct(jShll+2)=.False.
            pChrg(nCnttp)=.False.
            Fixed(nCnttp)=.True.
            FockOp(nCnttp)=.False.
            nOpt(nCnttp) = iOptn
            ipVal(nCnttp) = ipVal_
            ipPrj(nCnttp) = ipPrj_
            ipSRO(nCnttp) = ipSRO_
            ipSOC(nCnttp) = ipSOC_
            ipPP(nCnttp)  = ipPP_
            nVal_Shells(nCnttp) = nVal
            nPrj_Shells(nCnttp) = nPrj
            nSRO_Shells(nCnttp) = nSRO
            nSOC_Shells(nCnttp) = nSOC
            nPP_Shells(nCnttp)  = nPP
            nTot_Shells(nCnttp) = nVal+nPrj+nSRO+nSOC+nPP
            CntMass(nCnttp) = rMass(iAtmNr(nCnttp))
            Do iSh = jShll+1,iShll
              FragShell(iSh)=.True.
            End Do
            If (INDEX(sBasis(1:Indx-1),'6-31G').ne.0) Then
              Do iSh = jShll+3, iShll
                Prjct(iSh)=.False.
                Transf(iSh)=.False.
              End Do
            End If
            FragCnttp(nCnttp)=.True.
* add the coordinates (1 atom / basis set center)
            dbsc(nCnttp)%nCntr = 1
            mdc = mdc + 1
            If (mdc.gt.Mxdc) Then
              Write (6,*) ' FragExpand: Increase Mxdc'
              Write (6,*) '        Mxdc=',Mxdc
              Call ErrTra
              Call Quit_OnUserError()
            End If
* get the relative coordinates
            x1 = dbsc(iCnttp)%FragCoor(2,iAtom)
            y1 = dbsc(iCnttp)%FragCoor(3,iAtom)
            z1 = dbsc(iCnttp)%FragCoor(4,iAtom)
* make them absolute
            x1 = x1 + dbsc(iCnttp)%Coor(1,iCntr)
            y1 = y1 + dbsc(iCnttp)%Coor(2,iCntr)
            z1 = z1 + dbsc(iCnttp)%Coor(3,iCntr)
#ifdef _DEBUG_
            write(6,'(a,i3,3(a,F12.7))') 'FragExpand: Center ',nCnttp,
     &      ' Coordinates:  x =',x1,' y=',y1,' z=',z1
#endif
* store them
            Call mma_allocate(dbsc(nCnttp)%Coor,3,1,Label='dbsc:C')
            dbsc(nCnttp)%Coor(1,1) = x1
            dbsc(nCnttp)%Coor(2,1) = y1
            dbsc(nCnttp)%Coor(3,1) = z1
* store the Mulliken charge
            FragCharge(nCnttp) = dbsc(iCnttp)%FragCoor(5,iAtom)
* create custom (hopefully) unique label
            LenLbl = Index(sBasis,'.') - 1
            label = sBasis(1:LenLbl)
            Do ii = LenLbl+1,4
              label(ii:ii) = '_'
            End Do
c LENIN possible BUG
            LblCnt(mdc) = label
            If(mdc.lt.10) then
              write(label,'(a3,i1)') '___',mdc
            Else If(mdc.lt.100) then
              write(label,'(a2,i2)') '__',mdc
            Else If(mdc.lt.1000) then
              write(label,'(a1,i3)') '_',mdc
            Else
              write(label,'(i4)') mdc
            End If
c            LblCnt(mdc)(LENIN1:LENIN4) = label
            LblCnt(mdc)(5:LENIN2) = label
            Call ChkLbl(LblCnt(mdc),LblCnt,mdc-1)
* store a reference to the originating fragment placeholder
* misuse nFragCoor for this purpose: it will not overwrite anything, but
* beware of redimensioning this array to anything other than Mxdbsc
* To signify this we store the negative value such that we can identify
* the that the actual number of centers is 0 and that the corresponding
* size of dbsc()%FragCoor is 0 and nothing else.
            dbsc(mdc)%nFragCoor = -ndc
*
            nInfo = ipExp(iShll+1) - 1
            ninfo_stupid = nInfo
            If (ExpNuc(nCnttp).lt.Zero) ExpNuc(nCnttp) =
     &        NucExp(iMostAbundantIsotope(iAtmNr(nCnttp)))
          End Do
        End Do
      End Do ! iCnttp
*
c      write(6,'(A,i6)') 'FragExpand: after the main DO loop'
c      write(6,'(A,i5)') 'FragExpand:  nCnttp = ',nCnttp
c      write(6,'(A,10i6)') 'FragExpand: mdciCnttp(nCnttp)  = ',
c     &                                    (mdciCnttp(i),i=1,nCnttp)
c      write(6,'(A,10i6)') 'FragExpand:     dbsc(nCnttp)%nCntr  = ',
c     &                                  (dbsc(i)%nCntr,i=1,nCnttp)
c      write(6,'(A,10I6)') 'FragExpand: nFragType(nCnttp)  = ',
c     &                                  (dbsc(i)%nFragType,i=1,nCnttp)
c
c      write(6,'(A,10I6)') 'FragExpand: nFragCoor(nCnttp)  = '
c     &                                  (dbsc(i)%nFragCoor,i=1,nCnttp)
*                                                                      *
************************************************************************
*                                                                      *
      Mx_Shll=iShll+1
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
