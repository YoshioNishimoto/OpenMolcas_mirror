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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               2017, Roland Lindh                                     *
!***********************************************************************
      SubRoutine RdInp_scf()
!***********************************************************************
!                                                                      *
!     purpose: Read input options                                      *
!                                                                      *
!     called from: ReadIn                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!     modified by M.Schuetz @teokem.lu.se, 1995                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: UHF- V.Veryazov, 2003                                   *
!                                                                      *
!***********************************************************************
      use OccSets, only: OccSet_e, OccSet_M, nOccSet_e, nOccSet_m
      use KSDFT_Info, only: CoefR, CoefX
      use OFembed, only: OFE_KSDFT, dfmd, Do_OFemb, KEonly, ThrFThaw, XSigma
      use Functionals, only: Custom_File, Custom_Func
      use IOBuf, only: lDaRec,nSect!,DiskMx_MByte
      use InfSO, only: DltnTh, QNRTh
#ifdef _HDF5_
      use mh5, only: mh5_is_hdf5, mh5_open_file_r
      use InfSCF, only: FileOrb_ID
#endif
      use Fock_util_global, only: Deco, DensityCheck, Estimate, Update
!
      use SpinAV, only: Do_SpinAV
      use InfSCF, only: nIter, nAufb, AddFragments, Aufb, C1DIIS, Damping, DDnOff, DelThr, DIIS,      &
                        DIISTh, DoCholesky, DoHLgap, DSCF, DThr, EThr, FThr, Falcon, FckAuf,              &
                        FlipThr, ExFac, FckAuf, HLgap, iAu_ab, iCoCo, iDKeep, InVec, iPrForm, iPrint, iPrOrb,    &
                        isHDF5, iStatPRN, Iter2run, IterPrlv, nD, ivvloop, jPrint, jVOut, kIVO, klockan,       &
                        kOptim_Max, KSDFT, LKon, MaxFlip, MiniDn, MSYMON, MxConstr, nCore, nDisc, Neg2_Action,   &
                        NoExchange, NoProp, nSym, nTit, One_Grid, OnlyProp, PmTime, PreSch, QudThr, RFPert,      &
                        RGEK, RotFac, RotLev, RotMax, RSRFO, RTemp, SCF_FileOrb, ScrFac, Scrmbl, Teee, TemFac,   &
                        Thize, ThrEne, Tot_Charge, Tot_Nuc_Charge, TStop, WrOutD, nConstr, nOrb, nBas,           &
                        Tot_El_Charge, Title, nOcc, nDel, nFro, LstVec, indxc
      use Cholesky, only: ChFracMem, timings
      use ChoSCF, only: ALGO, dmpk, nScreen, ReOrd
      use MxDM, only: MxIter, MxOptm
      use AddCorr, only: Addc_KSDFT, Do_Addc, Do_Tw
      use ChoAuf, only: Cho_Aufb
      use Constants, only: Zero, Half, One, Two
      use stdalloc, only: mma_allocate
      Implicit None
!
#include "hfc_logical.fh"
!
!---- Define local variables
      Integer i, iAuf, iD, iFroz, iOccu, iOrbi, iPri, iStatus, iSym, j, KeywNo, lthSet_a, lthSet_b, LuCF, LuSpool, nFunc, &
              nnn, nSqrSum
      Real*8 Tot_Ml_Charge
      Integer, External:: Allocdisk, IsFreeUnit
      Real*8, External:: Get_ExFac
      Character(LEN=180)  Key, Line
      Character(LEN=180), External :: Get_Ln
      Integer iArray(32)
      Logical lTtl, IfAufChg,OccSet,FermSet,CharSet,UHFSet,SpinSet
      Logical Chol
      Integer Mode(1)
      character Method*8
      Logical TDen_UsrDef

!     copy input from standard input to a local scratch file
!
      Call SpoolInp(LuSpool)
!
      OccSet=.false.
      FermSet=.false.
      CharSet=.false.
      SpinSet=.false.
!
      Neg2_Action='STOP'
!
! Some initialization
!
      Chol=.false.
      ALGO  = 4
      REORD =.false.
      DECO  =.true.
      DensityCheck=.false.
      timings=.false.
      UHFSet=.false.
      Nscreen = 10    ! default screening interval (# of red sets)
      dmpk = 0.045d0   ! default damping of the screening threshold
      Estimate=.false.
      Update=.true.
#ifdef _MOLCAS_MPP_
      ChFracMem=0.3d0
#else
      ChFracMem=half
#endif
      SCF_FileOrb='INPORB'
      isHDF5=.False.
! Constrained SCF initialization
      Do i=1,nSym
         nConstr(i)=0
      End Do
      MxConstr=0
      klockan=1
      Do_Addc=.false.
      iTer2run=2
! Delta_Tw correlation energy calculation
      Do_Tw=.false.
! Read Cholesky info from runfile and save in infscf.fh
      Call DecideOnCholesky(DoCholesky)
      if (DoCholesky) then
         Chol=.True.
         Call Cho_scf_rdInp(.True.,LuSpool)
         if(ALGO.ge.2)then
            DDnOFF = .True. !do not use diffential density
            MiniDn = .False.
         endif
!tbp, may 2013: no thr modification with Cholesky
!tbp     Call Get_dScalar('Cholesky Threshold',ThrCom)
!tbp     EThr = Max(EThr,ThrCom)
      else
         DDnOFF = .false. !default for conventional
         MiniDn = .true.
      endif
      TDen_UsrDef=.False.

!---- Set up number of orbitals
      Do iSym = 1, nSym
         nOrb(iSym) = nBas(iSym)
      End Do
!---- Set up some counters
      nSqrSum=0
      Do iSym = 1, nSym
         nSqrSum=nSqrSum+nBas(iSym)*nBas(iSym)
      End Do
!
      iOrbi = 0
      iFroz = 0
      iOccu = 0
      nTit  = 0
      ivvloop=0
      iPrForm=-1
      iterprlv=0
      ScrFac=Zero
!
!---- Parameters that control how new orbitals
!---- are generated by neworb.
!
      RotLev=Zero
      RotFac=One
      RotMax=1.0d1
!     HLgap=-One
      HLgap= 0.2d0
      DoHLgap=.false.
      MaxFlip=10
      FlipThr=0.1
!
! Skip exchange when building Fock matrix
! (for debugging purposes)
!
      NoExchange=.False.
!
!---- Default value for starting orbitals
!---- Invec=-1 indicate decision in sorb
!
      InVec=-1
!
!---- Default to aufbau for neutral species
!
      Aufb   = .True.
      Teee   = .True.
      RTemp  = 0.500d0
      TemFac = 0.460d0
      TStop  = 0.010d0
      nAufb(1) = -1
      nAufb(2) = -1
      Call Get_dScalar('Total Nuclear Charge',Tot_Nuc_Charge)
      Tot_El_Charge=Zero
      Tot_Ml_Charge=Zero
      Tot_Charge=Zero
      iAu_ab=0
      iStatPRN=0
!
      IfAufChg=.False.
!
      Falcon=.False.
      MSYMON=.False.
!
      nD = 1
!
!---- Locate "start of input"
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'SCF')
!
      KeywNo=0
 1000 lTtl=.False.
  999 Continue
      Key = Get_Ln(LuSpool)
      Line = Key
      Call UpCase(Line)
      KeywNo=KeywNo+1
!
      If (Line(1:4).eq.'TITL') Go To 1100
      If (Line(1:4).eq.'ITER') Go To 1200
      If (Line(1:4).eq.'OCCU') Go To 1300
      If (Line(1:4).eq.'ORBI') Go To 1400
      If (Line(1:4).eq.'FROZ') Go To 1500
      If (Line(1:4).eq.'OVLD') Go To 1700
      If (Line(1:4).eq.'PRLS') Go To 1800
      If (Line(1:4).eq.'PROR') Go To 1900
      If (Line(1:4).eq.'KEEP') Go To 2000
      If (Line(1:4).eq.'STAR') Go To 2100
!
! Section for Cholesky input
      If (Line(1:4).eq.'CHOL') Go To 20000
      If (Line(1:4).eq.'CHOI') Go To 20001
!
      If (Line(1:4).eq.'CONS') Go To 20002
      If (Line(1:4).eq.'CORE') Go To 21000
      If (Line(1:4).eq.'NDDO') Go To 21001
      If (Line(1:4).eq.'LUMO') Go To 21002
      If (Line(1:4).eq.'GSSR') Go To 21003
      If (Line(1:4).eq.'REST') Go To 21004
      If (Line(1:4).eq.'THRE') Go To 2200
      If (Line(1:4).eq.'NODI') Go To 2300
      If (Line(1:4).eq.'DIIS') Go To 2400
      If (Line(1:4).eq.'OCCN') Go To 2500
      If (Line(1:4).eq.'MCCN') Go To 2510
      If (Line(1:4).eq.'IVO ') Go To 2600
      If (Line(1:4).eq.'UHF ') Go To 2700
      If (Line(1:4).eq.'HFC ') Go To 2701
      If (Line(1:4).eq.'NODA') Go To 2900
      If (Line(1:4).eq.'CONV') Go To 3000
      If (Line(1:4).eq.'DISK') Go To 3100
      If (Line(1:4).eq.'THIZ') Go To 3200
      If (Line(1:4).eq.'SIMP') Go To 3300
      If (Line(1:4).eq.'NOMI') Go To 3400
      If (Line(1:4).eq.'TDEN') Go To 3401
      If (Line(1:4).eq.'WRDE') Go To 3500
      If (Line(1:4).eq.'C1DI') Go To 3600
      If (Line(1:4).eq.'QUAD') Go To 3700
      If (Line(1:4).eq.'RS-R') Go To 3750
      If (Line(1:4).eq.'S-GE') Go To 3751
      If (Line(1:4).eq.'SCRA') Go To 3800
      If (Line(1:4).eq.'EXTR') Go To 3900
      If (Line(1:4).eq.'RFPE') Go To 4000
      If (Line(1:4).eq.'QNRT') Go To 4100
      If (Line(1:4).eq.'AUFB') Go To 4200
      If (Line(1:4).eq.'FERM') Go To 4250
      If (Line(1:4).eq.'TEEE') Go To 4300
      If (Line(1:4).eq.'CHAR') Go To 4400
      If (Line(1:4).eq.'NOTE') Go To 4500
      If (Line(1:4).eq.'KSDF') Go To 4600
      If (Line(1:4).eq.'DFCF') Go To 4605
      If (Line(1:4).eq.'OFEM') Go To 4650
      If (Line(1:4).eq.'FTHA') Go To 4655
      If (Line(1:4).eq.'DFMD') Go To 4656
      If (Line(1:4).eq.'KEON') Go To 4660
      If (Line(1:4).eq.'TWCO') Go To 4661
      If (Line(1:4).eq.'ADDC') Go To 4662
      If (Line(1:4).eq.'SAVE') Go To 4663
      If (Line(1:4).eq.'DEBU') Go To 4700
      If (Line(1:4).eq.'ZSPI') Go To 4800
      If (Line(1:4).eq.'SPIN') Go To 4850
      If (Line(1:4).eq.'EXFA') Go To 4900
      If (Line(1:4).eq.'ONEG') Go To 4901
      If (Line(1:4).eq.'ROTP') Go To 5000
      If (Line(1:4).eq.'HLGA') Go To 5002
      If (Line(1:4).eq.'FLIP') Go To 5020
      If (Line(1:4).eq.'PMTI') Go To 6000
      If (Line(1:4).eq.'STAT') Go To 6010
      If (Line(1:4).eq.'ADDF') Go To 6020
      If (Line(1:4).eq.'FILE') Go To 6030
      If (Line(1:4).eq.'ITPR') Go To 7100
      If (Line(1:4).eq.'PROP') Go To 7200
      If (Line(1:4).eq.'NOPR') Go To 7201
      If (Line(1:4).eq.'NOX ') Go To 8700
      If (Line(1:4).eq.'PSDC') Go To 8900
      If (Line(1:4).eq.'USEX') Go To 8901
      If (Line(1:4).eq.'NEG2') Go To 8903
      If (Line(1:4).eq.'MSYM') Go To 8904
      If (Line(1:4).eq.'ITDI') Go To 8905
      If (Line(1:4).eq.'FCKA') Go To 8906
      If (Line(1:4).eq.'DEPT') Go To 8907
!
      If (Line(1:4).eq.'FALC') Go To 30000
!
      If (Line(1:4).eq.'END ') Go To 9000

      If (lTtl) Go To 1101
      Write (6,*) 'Unidentified key word:', Key
      call FindErrorLine
      Call Quit_OnUserError()
!
!>>>>>>>>>>>>> TITL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1100 Continue
      Line = Get_Ln(LuSpool)
      lTtl=.True.
 1101 Continue
      nTit = nTit + 1
      if(nTit.eq.1) Then
        Title(1)=Line(1:72)
      else
        if(nTit.eq.2) then
        call WarningMessage(1,'More than one title line!')
        endif
      endif
!      If (nTit.le.MxTit) Then
!         Title(nTit) = Line(1:72)
!      End If
      GoTo 999
!
!>>>>>>>>>>>>> ITER <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1200 Continue
      Line=Get_Ln(LuSpool)
      Line(179:180)='-1'
      Call Put_Ln(Line)
      Call Get_I(1,nIter,2)
      If(nIter(1).eq.-1) nIter(1)=nIter(0)
      GoTo 1000
!
!>>>>>>>>>>>>> OCCU <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1300 Continue
      If(FermSet) Then
         call WarningMessage(2,'Options OCCUpied and FERMi are mutually exclusive')
         Call Abend()
      End If
      If(SpinSet) Then
         call WarningMessage(2,'Keyword OCCUpied and ZSPIn/SPIN are mutually exclusive')
         Call Abend()
      End If
      If(CharSet) Then
         call WarningMessage(2,'Options OCCUpied and CHARGE are mutually exclusive')
         Call Abend()
      End If
      Do iD = 1, nD
         Line=Get_Ln(LuSpool)
         Call Get_I(1,nOcc(1,iD),nSym)
         If (iD.eq.1) Then
            Call Put_iArray('nIsh',nOcc(1,iD),nSym)
         Else
            Call Put_iArray('nIsh beta',nOcc(1,iD),nSym)
         End If
      End Do
!
      iOccu = 1
      if(nD==1) then
      Do i = 1, nSym
         Tot_El_Charge=Tot_El_Charge-DBLE(2*nOcc(i,1))
      End Do
      else
      Do i = 1, nSym
         Tot_El_Charge=Tot_El_Charge-DBLE(nOcc(i,1)+nOcc(i,2))
      End Do
      endif
      Aufb=.False. ! Disable default action.
      Teee=.False.
      OccSet=.true.
      Cho_Aufb=.false.
      UHFSet=.true.
      GoTo 1000
!
!>>>>>>>>>>>>> ORBI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1400 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nOrb,nSym)
      iOrbi = 1
      Do iSym=1,nSym
         nDel(iSym)=nBas(iSym)-nOrb(iSym)
      End Do
      GoTo 1000
!
!>>>>>>>>>>>>> FROZ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1500 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nFro,nSym)
      Call Put_iArray('nFro',nFro,nSym)
      iFroz = 1
      GoTo 1000
!
!>>>>>>>>>>>>> OVLD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1700 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,DelThr)
      GoTo 1000
!
!>>>>>>>>>>>>> PRLS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1800 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,iPri)
      iPrint = Max(iPri,iPrint)
      GoTo 1000
!
!>>>>>>>>>>>>> PROR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1900 Continue
      Line=Get_Ln(LuSpool)
      Line(179:180)='-1'
      Call Put_Ln(Line)
      Call Get_I1(1,iPrOrb)
      If (iPrOrb.ge.2) then
           Call Get_F1(2,ThrEne)
           Call Get_I1(3,iPrForm)
      else
            Call Get_I1(2,iPrForm)
      endif
      GoTo 1000
!
!>>>>>>>>>>>>> KEEP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,iDKeep)
      GoTo 1000
!
!>>>>>>>>>>>>> STAR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2100 Continue
      Line=Get_Ln(LuSpool)
      Line(157:180)='-1 -1 -1 -1 -1 -1 -1 -1'
      Call Put_Ln(Line)
      Call Get_I(1,LstVec,7)
! temporary hack to use density
      If(LstVec(1).eq.3) Then
         InVec=3
         RTemp  = 0.100d0
         TemFac = 0.100d0
         TStop  = 0.005d0
!        Aufb=.false.
!        Teee=.false.
      End If
      GoTo 1000
!
!...... Explicit Start Options
!
!>>>>>>>>>>>>> CHOL <<<<<<< Cholesky Default Input <<<<<<<<<<<<<<<<<
!
20000 Continue
      Chol=.True.
      Call Cho_scf_rdInp(.True.,LuSpool)
       if(ALGO.ge.2)then
         DDnOFF = .True. !do not use diffential density
         MiniDn = .False.
       endif
      GoTo 1000
!
!>>>>>>>>>>>>> CHOI <<<<<<< Cholesky Custom  Input <<<<<<<<<<<<<<<<<
!
!     Cholesky with user-defined settings.
!
20001 Continue
      Chol=.True.
      DDnOFF = .False. ! reset to default value
      MiniDn = .True.  ! reset to default value
      Call Cho_scf_rdInp(.False.,LuSpool)
       if(ALGO.ge.2)then
         DDnOFF = .True. !do not use diffential density
         MiniDn = .False.
       else
         DECO=.false.
       endif
      GoTo 1000
!
!>>>>>>>>>>>>> CONS <<<<<<<<<<<< Constrained SCF <<<<<<<<<<<
20002 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nConstr,nSym)
      Call Izero(indxC,16*2*8)
      Do i=1,nSym
         MxConstr=Max(MxConstr,nConstr(i))
         Line=Get_Ln(LuSpool)
         Call Get_I(1,indxC(1,1,i),nConstr(i))
         If (nConstr(i).gt.16) Then
            write(6,*) ' Max nConstr=16. Increase 1st dim of indxC and recompile'
            Call Abend()
         EndIf
         Do j=1,nConstr(i)
            If (indxC(j,1,i).ne.1 .and. indxC(j,1,i).ne.-1) Then
               write(6,*) ' Only 1 and -1 are accepted values'
               Call Abend()
            ElseIf (indxC(j,1,i).eq.1) Then
               indxC(j,2,i)=2
            Else
               indxC(j,1,i)=2
               indxC(j,2,i)=1
            EndIf
         End Do
      End Do
      GoTo 1000
!
!>>>>>>>>>>>>> CORE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21000 Continue
      InVec=0
      LstVec(1)=4
      LstVec(2)=-1
      GoTo 1000
!>>>>>>>>>>>>> NDDO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21001 Continue
      If (Chol) Then
         call WarningMessage(1,'RdInp: NDDO and Cholesky are incompatible!!; NDDO Option ignored')
      Else
      InVec=1
      LstVec(1)=5
      LstVec(2)=-1
      EndIf
      GoTo 1000
!>>>>>>>>>>>>> LUMO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21002 Continue
      InVec=2
      One_Grid=.True.
      LstVec(1)=2
      LstVec(2)=-1
      GoTo 1000
!>>>>>>>>>>>>> FILE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
6030  Continue
      InVec=2
      One_Grid=.True.
      LstVec(1)=2
      LstVec(2)=-1
      Line=Get_Ln(LuSpool)
      call fileorb(Line,SCF_FileOrb)
#ifdef _HDF5_
      If (mh5_is_hdf5(SCF_FileOrb)) Then
         isHDF5=.True.
         fileorb_id=mh5_open_file_r(SCF_FileOrb)
      End If
#endif
      goto 1000
!>>>>>>>>>>>>> GSSR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21003 Continue
      InVec=9
      One_Grid=.True.
      LstVec(1)=1
      LstVec(2)=-1
      GoTo 1000
!>>>>>>>>>>>>> REST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21004 Continue
      Write (6,*) 'WARNING: Option REST is redundant.'
      GoTo 1000
!>>>>>>>>>>>>> THRE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2200 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,EThr)
      Call Get_F1(2,DThr)
      Call Get_F1(3,FThr)
      Call Get_F1(4,DltNTh)
!tbp, may 2013: no thr modification with Cholesky
!tbp  If (DoCholesky) then
!tbp     write(ww,'(a,es20.8)')'Detected Cholesky or RI/DF calculation BUT user specified EThr will be used. Ethr = ',EThr
!tbp     call WarningMessage(1,ww)
!tbp  EndIf
!     If (  DThr*1.D-2.gt.EThr) Then
!        Write (6,*)
!        Write (6,*) '----> WARNING! <----'
!        Write (6,*)
!        Write (6,*) ' The value of DThr is inconsistent with'
!        Write (6,*) ' with the value of EThr. The code will'
!        Write (6,*) ' automatically reset the value to something'
!        Write (6,*) ' more reasonable based on the requested'
!        Write (6,*) ' threshold of the energy.'
!        DThr=100.D+00*EThr
!     End If
!     If (  FThr*1.D-2.gt.EThr) Then
!        Write (6,*)
!        Write (6,*) '----> WARNING! <----'
!        Write (6,*)
!        Write (6,*) ' The value of FThr is inconsistent with'
!        Write (6,*) ' with the value of EThr. The code will'
!        Write (6,*) ' automatically reset the value to something'
!        Write (6,*) ' more reasonable based on the requested'
!        Write (6,*) ' threshold of the energy.'
!        FThr=100.D+00*EThr
!     End If
!     If (DltNTh*1.D-2.gt.Sqrt(EThr)) Then
!        Write (6,*)
!        Write (6,*) '----> WARNING! <----'
!        Write (6,*)
!        Write (6,*) ' The value of DltNTh is inconsistent with'
!        Write (6,*) ' with the value of EThr. The code will'
!        Write (6,*) ' automatically reset the value to something'
!        Write (6,*) ' more reasonable based on the requested'
!        Write (6,*) ' threshold of the energy.'
!        DltNTh=100.D+00*Sqrt(EThr)
!     End If
      GoTo 1000
!
!>>>>>>>>>>>>> NODI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2300 Continue
      Diis = .False.
      GoTo 1000
!
!>>>>>>>>>>>>> DIIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2400 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,DiisTh)
      GoTo 1000
!
!>>>>>>>>>>>>> OCCN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2500 Continue
      If (iOccu.ne.1) Then
         Call WarningMessage(2,' Input Error!;The OCCNumber option works only if  the OCCUpied option has been specified!')
         Call Abend()
      End If
      lthSet_a = 0
      lthSet_b = 0
      Do 2501 iSym = 1, nSym
         lthSet_a = lthSet_a + nOcc(iSym,1)
         If(nD==2) Then
            lthSet_b = lthSet_b + nOcc(iSym,2)
         End If
 2501 Continue
      nOccSet_e=max(lthSet_a,lthSet_b)
!     Write(6,'(a,i5)') 'rdinp: lthset_a ',lthset_a
!     Write(6,'(a,i5)') 'rdinp: lthset_b ',lthset_b
!
!---- Note, it is dangerous to read Line first. There may be many
!     lines with occupation numbers.
!
      Call mma_allocate(OccSet_e,nOccSet_e,nD,Label='OccSet_e')
      Call FZero(OccSet_e,nOccSet_e*nD)
      Read(LuSpool,*,End=902,Err=903) (OccSet_e(i,1),i=1,lthSet_a)
      If(nD==2) then
        Read(LuSpool,*,End=902,Err=903) (OccSet_e(i,2),i=1,lthSet_b)
      End If

      Tot_El_Charge=Zero
      Do iD = 1, nD
         Do i = 1, nOccSet_e
            Tot_El_Charge=Tot_El_Charge-OccSet_e(i,iD)
         End Do
      End Do
      UHFSet=.true.
      iCoCo=1
      GoTo 1000
!
!>>>>>>>>>>>>> MCCN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     Just as OCCN, but for muons
 2510 Continue
      If (iOccu.ne.1) Then
         Call WarningMessage(2,' Input Error!;The OCCNumber option works only if the OCCUpied option has been specified!')
         Call Abend()
      End If
      lthSet_a = 0
      lthSet_b = 0
      Do 2511 iSym = 1, nSym
         lthSet_a = lthSet_a + nOcc(iSym,1)
         If(nD==2) Then
            lthSet_b = lthSet_b + nOcc(iSym,2)
         End If
 2511 Continue
      nOccset_m=max(lthSet_a,lthSet_b)
!     Write(6,'(a,i5)') 'rdinp: lthset_a ',lthset_a
!     Write(6,'(a,i5)') 'rdinp: lthset_b ',lthset_b
!
!---- Note, it is dangerous to read Line first. There may be many
!     lines with occupation numbers.
!
      Call mma_allocate(OccSet_m,nOccSet_m,nD,Label='OccSet_m')
      Call FZero(OccSet_m,nOccSet_m*nD)
      Read(LuSpool,*,End=902,Err=903) (OccSet_m(i,1),i=1,lthSet_a)
      If(nD==2) then
        Read(LuSpool,*,End=902,Err=903) (OccSet_m(i,2),i=1,lthSet_b)
      End If

      Tot_Ml_Charge=Zero
      Do iD = 1, nD
         Do i = 1, nOccSet_m
            Tot_Ml_Charge=Tot_Ml_Charge-OccSet_m(i,iD)
         End Do
      End Do
      UHFSet=.true.
      iCoCo=1
      GoTo 1000
!
!>>>>>>>>>>>>> IVO  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2600 Continue
      kIvo = 1
      GoTo 1000
!
!>>>>>>>>>>>>> UHF  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2700 Continue
      if(UHFSet) then
      call sysAbendMsg('rdinp','incorrect input','UHF keyword should be placed before others')
      endif
      nD       = 2
      MiniDn = .False.
      GoTo 1000
!
!>>>>>>>>>>>>> HFC  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2701 Continue
      UHF_HFC     = .True.
      GoTo 1000
!
!>>>>>>>>>>>>> NODA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2900 Continue
      Damping=.False.
      GoTo 1000
!
!>>>>>>>>>>>>> CONV <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3000 Continue
      DSCF = .False.
      GoTo 1000
!
!>>>>>>>>>>>>> DISK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3100 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,nDisc)
      Call Get_I1(2,nCore)
      GoTo 1000
!
!>>>>>>>>>>>>> THIZ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3200 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,Thize)
      GoTo 1000
!
!>>>>>>>>>>>>> SIMP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3300 Continue
      PreSch = .True.
      GoTo 1000
!
!>>>>>>>>>>>>> NOMI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3400 Continue
      MiniDn = .False.  !do not use minimized density diff
      GoTo 1000
!
!>>>>>>>>>>>>> TDEN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3401 Continue
      DDnOFF = .True. !do not use diffential density
      MiniDn = .False.
      TDen_UsrDef=.True.
      GoTo 1000
!
!>>>>>>>>>>>>> WRDE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3500 Continue
      WrOutD = .True.
      GoTo 1000
!
!>>>>>>>>>>>>> C1DI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3600 Continue
      c1Diis = .True.
      GoTo 1000
!
!>>>>>>>>>>>>> QUAD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3700 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,QudThr)
      GoTo 1000
!
!>>>>>>>>>>>>> RS-R <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3750 Continue
      RSRFO = .True.
      RGEK  = .False.
      GoTo 1000
!
!>>>>>>>>>>>>> S-GE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3751 Continue
      RGEK  = .True.
      RSRFO = .False.
      GoTo 1000
!
!>>>>>>>>>>>>> SCRA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3800 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,ScrFac)
      Scrmbl = .True.
      GoTo 1000
!
!>>>>>>>>>>>>> EXTR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 3900 Continue
      call WarningMessage(1,'EXTRACT option is redundant and is ignored!')
      Goto 1000
!
!>>>>>>>>>>>>> RFPE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4000 Continue
      RFpert = .True.
      GoTo 1000
!
!>>>>>>>>>>>>> QNRT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4100 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,QNRTh)
      GoTo 1000
!
!>>>>>>>>>>>>> AUFB <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4200 Continue
      call WarningMessage(2,' Error: Keyword AUFBau is now obsolete!;Use keyword CHARge')
      Call Abend()
      if(Chol)then
        DDnOFF = .true.
        DECO   = .true.
        MiniDn = .false.
      endif
      Line=Get_Ln(LuSpool)
      Line(180:180)='2'
      Call Put_Ln(Line)
      Call Get_I(1,iArray,nD+1)
      Do i=1,nD
         nAufb(i)=iArray(i)
      EndDo
      iAuf=iArray(nD+2)
      If (IfAufChg) Then
      call WarningMessage(2,'Option AUFBau is mutually exclusive CHARge')
         Call Abend()
      End If
      if(nD==1) then
       Tot_El_Charge=-Two*DBLE(nAufb(1))
      else
       Tot_El_Charge=-DBLE(nAufb(1)+nAufb(2))
      endif
      IfAufChg=.True.
 4210 If(iAuf.eq.0) Then
         Teee = .False.
      Else If(iAuf.eq.1) Then
         RTemp  = 0.500d0
         TemFac = 0.400d0
         TStop  = 0.005d0
      Else If(iAuf.eq.2) Then
         RTemp  = 0.500d0
         TemFac = 0.460d0
         TStop  = 0.010d0
      Else If(iAuf.eq.3) Then
         RTemp  = 0.500d0
         TemFac = 0.610d0
         TStop  = 0.025d0
      Else If(iAuf.eq.4) Then
         RTemp  = 1.000d0
         TemFac = 0.730d0
         TStop  = 0.060d0
      Else If(iAuf.eq.5) Then
         RTemp  = 1.000d0
         TemFac = 0.870d0
         TStop  = 0.150d0
      Else
         call WarningMessage(1,' RdInp: Aufbau case must be in the range 0-5;Using case 5!')
         RTemp  = 1.000d0
         TemFac = 0.870d0
         TStop  = 0.150d0
      End If
      UHFSet=.true.
      GoTo 1000
!>>>>>>>>>>>>> FERM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4250 Continue
      If(OccSet) Then
         call WarningMessage(2,'Options OCCUpied and FERMi are mutually exclusive')
         Call Abend()
      End If
      if(Chol)then
        DDnOFF = .true.
        DECO   = .true.
        MiniDn = .false.
      endif
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,iAuf)
      FermSet=.true.
      GoTo 4210
!>>>>>>>>>>>>> TEEE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4300 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,RTemp)
      Call Get_F1(2,TemFac)
      Call Get_F1(3,TStop)
      If (TStop.lt.Zero) Then
      call WarningMessage(2,'Input Error!; End temperture < 0.0 ')
         Call Abend
      End If
      If (TemFac.lt.Zero) Then
      call WarningMessage(2,'Input Error!; Temperture factor < 0.0 ')
         Call Abend
      End If
      If (RTemp.lt.Zero) Then
      call WarningMessage(2,'Input Error!; Start temperature < 0.0 ')
         Call Abend
      End If
      If (RTemp.lt.TStop) Then
      call WarningMessage(2,'Input Error!; End temperture > start temperature ')
         Call Abend
      End If
      GoTo 1000
!>>>>>>>>>>>>> CHAR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4400 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,iArray,1)
      Tot_Charge=DBLE(iArray(1))
      nAufb(1)=-1
      nAufb(2)=-1
      If (IfAufChg) Then
      call WarningMessage(2,'Input Error!; Option CHARge is mutually exclusive to AUFBau')
         Call Abend()
      End If
      If (OccSet) Then
      call WarningMessage(2,'Input Error!; Option CHARge is mutually exclusive to OCCUpied')
         Call Abend()
      End If
      IfAufChg=.True.
      CharSet=.True.
      GoTo 1000
!>>>>>>>>>>>>> NOTE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4500 Continue
      Teee=.False.
      GoTo 1000
!
!>>>>>>>>>>>>> KSDF <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4600 Continue
      Line=Get_Ln(LuSpool)
      Call UpCase(Line)
      Line = adjustl(Line)
      KSDFT=Line(1:80)
      nFunc = 0
      Read(Line,*,iostat=istatus) nFunc
      If ((istatus == 0) .and. (nFunc > 0)) Then
        KSDFT = Custom_Func
        LuCF = IsFreeUnit(10)
        Call molcas_open(LuCF,Custom_File)
        Write(LuCF,*) Trim(KSDFT),nFunc
        Do i=1,nFunc
          Line=Get_Ln(LuSpool)
          Write(LuCF,*) Trim(Line)
        End Do
        Close(LuCF)
      End If
      GoTo 1000
!
!>>>>>>>>>>>>> DFCF <<<< Factors to scale exch. and corr. <<
 4605 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,CoefX)
      Call Get_F1(2,CoefR)
!      Call put_dscalar('DFT exch coeff',CoefX)
!      Call put_dscalar('DFT corr coeff',CoefR)
      GoTo 1000
!
!>>>>>>>>>>>>> OFEM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4650 Continue
      Do_OFemb=.true.
      Line=Get_Ln(LuSpool)
      Call UpCase(Line)
      Line = adjustl(Line)
      OFE_KSDFT=Line(1:16)
      write(6,*)  '  --------------------------------------'
      write(6,*)  '   Orbital-Free Embedding Calculation'
      write(6,*)  '  --------------------------------------'
      If (OFE_KSDFT(1:4).eq.'LDTF') Then
       write(6,*) '    T_nad potential   : Thomas-Fermi    '
      Else
       write(6,*) '    T_nad potential   : ',OFE_KSDFT(1:4)
      EndIf
      If (KEonly) Then
       write(6,*) '    Exc_nad potential :  None           '
      Else
       write(6,*) '    Exc_nad potential : ',OFE_KSDFT(6:10)
      EndIf
      write(6,*)  '  --------------------------------------'
      write(6,*)
      GoTo 1000
!
!>>>>>>>>>>>>> FTHA <<<< threshold for Freeze-n-Thaw <<<<<<<
 4655 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,ThrFThaw)
      GoTo 1000
!
!>>>>>>>>>>>>> DFMD <<<< fraction of correlation potential <
 4656 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,dFMD)
      Call Get_F1(2,Xsigma)
      If (dFMD+Xsigma.lt.Zero) Then
       write(6,*)' *** Warning: arguments to DFMD must be nonnegative!'
       write(6,*)' ***          I will take their ABS !!! '
       dFMD=abs(dFMD)
       Xsigma=abs(Xsigma)
      EndIf
      GoTo 1000
!
!>>>>>>>>>>>>> KEON <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4660 Continue
      KEonly=.true.
      If (.not.Do_OFemb) Then
       write(6,*) ' *** Warning:  KEonly works in OFembedding runs!'
      Else
       write(6,*) ' *** Warning:  Exc_nad set to NONE at this point ***'
      EndIf
      write(6,*)
      GoTo 1000
!
!>>>>>>>>>>>>> TWCO <<<<< activate Tw correlation <<<<<<<<<<
 4661 Continue
      Do_Tw =.true.
      GoTo 1000
!
!>>>>>>>>>>>>> ADDC << add correlation energy (CONStraint) <
 4662 Continue
      Do_Addc=.True.
      Line=Get_Ln(LuSpool)
      Call UpCase(Line)
      Line = adjustl(Line)
      ADDC_KSDFT=Line(1:80)
      GoTo 1000
!
!>>>>>>>>>>>>> SAVE << Spin-Averaged wavelets (CONStraint) <
 4663 Continue
      Do_SpinAV=.True.
      GoTo 1000
!
!>>>>>>>>>>>>> DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4700 Continue
       ivvloop=1
       Diis=.false.
       MiniDn=.false.
       Damping=.False.
      GoTo 1000
!>>>>>>>>>>>>> ZSPI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4800 Continue
      If (SpinSet) Then
         Call WarningMessage(2,'Multiple definition of SPIN/ZSPIn')
         Call Abend()
      End If
      If(OccSet) Then
      call WarningMessage(2,'Input Error!; Keywords OCCUpied and ZSPIn are mutually exclusive')
         Call Abend()
      End if
      Line=Get_Ln(LuSpool)
      Call Get_I(1,iArray,1)
      iAu_ab=iArray(1)
      If ((nD/=2).and.(iAu_ab.ne.0)) Then
         Call WarningMessage(2,'ZSPIn different from 0 requires UHF before it')
         Call Abend()
      End If
      SpinSet=.true.
      GoTo 1000
!>>>>>>>>>>>>> SPIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4850 Continue
      If (SpinSet) Then
         Call WarningMessage(2,'Multiple definition of SPIN/ZSPIn')
         Call Abend()
      End If
      If(OccSet) Then
      call WarningMessage(2,'Input Error!; Keywords OCCUpied and SPIN are mutually exclusive')
         Call Abend()
      End if
      Line=Get_Ln(LuSpool)
      Call Get_I(1,iArray,1)
      iAu_ab=iArray(1)-1
      If (iAu_ab.lt.0) Then
         Call WarningMessage(2,'SPIN must be a positive integer')
         Call Abend()
      End If
      If (iAu_ab/=0) Then
         nD = 2
         MiniDn = .False.
      End If
      If ((nD/=2).and.(iAu_ab.ne.0)) Then
         Call WarningMessage(2,'SPIN greater than 1 requires UHF before it')
         Call Abend()
      End If
      SpinSet=.true.
      GoTo 1000
!>>>>>>>>>>>>> EXFA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4900 Continue
      ExFac=Zero
      GoTo 1000

!>>>>>>>>>>>>> ONEG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 4901 Continue
      One_Grid=.True.
      GoTo 1000

!>>>>>>>>>>>>> ROTP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 5000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,RotLev)
      Call Get_F1(2,RotFac)
      Call Get_F1(3,RotMax)
      Write(6,'(a,ES15.3)') 'Fock matrix levelshift   ',RotLev
      Write(6,'(a,ES15.3)') 'Fock matrix scaling      ',RotFac
      Write(6,'(a,ES15.3)') 'Fock matrix max rotation ',RotMax
      GoTo 1000
!>>>>>>>>>>>>> HLGA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 5002 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,HLgap)
      Write(6,'(a,ES15.3)') 'Minimum HOMO-LUMO gap    ',HLgap
      DoHLgap=.true.
      QNRTh    = Zero
      GoTo 1000
!>>>>>>>>>>>>> FLIP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 5020 Continue
      Line=Get_Ln(LuSpool)
      Line(178:180)='0.1'
      Call Put_Ln(Line)
      Call Get_I(1,iArray,1)
      Call Get_F1(2,FlipThr)
      MaxFlip=iArray(1)
!     Write(6,*) 'MaxFlip:',MaxFlip
!     Write(6,*) 'FlipThr:',FlipThr
      Goto 1000
!>>>>>>>>>>>>> PMTI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Time (CPU *and* wall) subroutine pmat_scf (2-el Fock matrix)
 6000 Continue
      PmTime=.true.
      GoTo 1000
!>>>>>>>>>>>>> STAT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Print Statistic information
 6010 Continue
      iStatPRN=1
      GoTo 1000
!>>>>>>>>>>>>> ADDF <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Add the fragment atoms to the MOLDEN file
 6020 Continue
      AddFragments = .True.
      Goto 1000
!>>>>>>>>>>>>> ITPR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 7100 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,iterprlv)
      Goto 1000
!>>>>>>>>>>>>> PROP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 7200 Continue
      InVec=2
      One_Grid=.True.
      LstVec(1)=2
      LstVec(2)=-1
      OnlyProp=.true.
      Goto 1000
!>>>>>>>>>>>>> SKIP PROP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 7201 Continue
      NoProp=.true.
      Goto 1000
!>>>>>>>>>>>>> NOX  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Debug option: skip exchange in Fock matrix build
 8700 Continue
      NoExchange=.True.
      Goto 1000
!>>>>>>>>>>>>> PSDC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Debug option: check that full integral matrix is PSD
! input=1: diagonalization
! input=2: Cholesky decomposition
 8900 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,Mode,1)
      Goto 1000
!>>>>>>>>>>>>> USEX <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Debug option: use exact diagonal (1) or off-dagonal (2) blocks when
! checking PSD
 8901 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,Mode,1)
      Goto 1000
!>>>>>>>>>>>>> NEG2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Specify action when negative two-electron energies are
! encountered (stop, warn, continue).
 8903 Continue
      Line=Get_Ln(LuSpool)
      Call UpCase(Line)
      Neg2_Action=Line(1:4)
      If (Neg2_Action.ne.'STOP' .and. Neg2_Action.ne.'WARN' .and. Neg2_Action.ne.'CONT') Then
         Write(6,'(A,A)') 'Illegal input for NEG2 keyword: ',Line(1:4)
!        Call FindErrorLine()
         Call Quit_OnUserError()
      End If
      Goto 1000
!>>>>>>>>>>>>> MSYM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8904 Continue
      MSYMON = .True.
      GoTo 1000
!>>>>>>>>>>>>> ITDIIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8905 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,iTer2run)
      GoTo 1000
!>>>>>>>>>>>>> ITDIIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8906 Continue
      Line=Get_Ln(LuSpool)
      Call UpCase(Line)
      If (Index(Line,'TRUE').ne.0) Then
         FckAuf=.True.
      Else
         FckAuf=.False.
      End If
      GoTo 1000
!>>>>>>>>>>>>> DEPTH  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8907 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,kOptim_Max)
      If (kOptim_Max>MxOptm) Then
         Write (6,*) 'kOptim_Max>MxOptm'
         Write (6,*) 'kOptim_Max=',kOptim_Max
         Write (6,*) 'MxOptm=',MxOptm
         Write (6,*) 'Modify mxdm.f90 and recompile!'
         Call Abend()
      End If
      GoTo 1000
!>>>>>>>>>>>>> FALC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
30000 Continue
      Falcon = .True.
      GoTo 1000
!
!>>>>>>>>>>>>> END  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 9000 Continue
!
      If (iPrint.ge.3) iStatPRN=1
!
      Tot_El_charge=Tot_El_Charge+Tot_Ml_Charge
!
! xml tag method
!
      If(nD==1) Then
         If(KSDFT.eq.'SCF') Then
            Call xml_cDump('method','','',0,'rhf',1,1)
         Else
            Call xml_cDump('method','','',0,'rdft',1,1)
         End If
      Else
         If(KSDFT.eq.'SCF') Then
            Call xml_cDump('method','','',0,'uhf',1,1)
         Else
            Call xml_cDump('method','','',0,'udft',1,1)
         End If
      End If
!
! Even or odd number of electrons
!
      If(.not.SpinSet) Then
         nnn=Int(tot_nuc_charge-tot_charge+half)
         If((nnn/2)*2.ne.nnn) Then
            iAu_ab=1
         End If
      End If
!
! Check start orbital priority list
!
      If(.not.OccSet .and. .not.FermSet) Then
!        Write(6,*) 'rdinp: Checking OCCU/FERM'
         Call VecFind(OccSet,FermSet,CharSet,SpinSet)
         If(OccSet .and. .not.FermSet) Then
!           Write(6,*) 'Using OCCU'
            Aufb=.false.
            Teee=.false.
            Cho_Aufb=.false.
         Else If(FermSet .and. .not.OccSet) Then
!           Write(6,*) 'Using FERM'
            Aufb=.true.
            Teee=.true.
            Cho_Aufb=.true.
         Else
            call WarningMessage(2,'Internal (logical) error in rdinp_scf')
            Call Abend()
         End If
      End If
      If(SpinSet) Then
         nnn=Int(Tot_Nuc_Charge-Tot_Charge-iAu_ab+half)
         If((nnn/2)*2.ne.nnn) Then
         call WarningMessage(2, 'Input error!;zSpin inconsistent with number of electrons')
            Call Abend()
         End If
      End If
!
!---- Check if certain input parameters are not in conflict
!
      If (iCoCo.eq.1.and.iOccu.eq.0) Then
         call WarningMessage(2, 'Input error!; The OCCNumber option require that the OCCUpied option is specified!')
         Call Abend
      End If
!
      If (Max(nIter(0),nIter(1)).gt.MxIter) Then
         call WarningMessage(1, 'Input error!;Number of iterations specified in input exceeds allowed maximum!')
!         Write (6,*) 'nIter(0)=',nIter(0)
!         Write (6,*) 'nIter(1)=',nIter(1)
!         Write (6,*) 'MxIter=',MxIter
!         Write (6,*)
         nIter(0)=MxIter
         nIter(1)=MxIter

!         Write (6,*) 'Number of iteration reset to ',MxIter
!         Write (6,*)
      End If
!
      If (iOrbi.eq.1 .and. (InVec.ne.2 .and. InVec.ne.4)) Then
         call WarningMessage(2, 'Input error!; The ORBITAL option can only be used with input orbitals!;'//   &
                                ' Possible exclusive options are: LUMORB RESTART')
         Call Abend
      End If
!
      If (iFroz.eq.1 .and. InVec.eq.3) Then
       call WarningMessage(2, 'Input error!; The FROZEN option does not work with an density matrix input')
         Call Abend
      End If
!
      If (InVec.eq.1 .and. Aufb .and. .not.Chol) Then
         call WarningMessage(2, 'Input error!; Aufbau not compatible with NDDO option')
         Call Abend
      End If
!
      If (Invec.lt.-1 .or. InVec.gt.9) Then
         call WarningMessage(2, 'Input error!; inappropriate value for Invec')
         Call Abend
      End If
!
      If(nD==1 .and. UHF_HFC) Then
      call sysAbendMsg('rdinp','incorrect input','HFC keyword should be used with UHF')
      End If
!
!---- Print out warning informations (if any)
!
      If (iFroz.eq.1 .and. (InVec.eq.2 .or. InVec.eq.4)) Then
         call WarningMessage(1,'RdInp: Warning!;Option FROZEN is used together with input orbitals as'//  &
                               ' initial guess.;Freezing of orbitals is performed at MO level.;'//        &
                               'First orbitals in each symmetry will not be modified.')
      End If
!
      If (Aufb .AND.(iFroz.eq.1)) Then
         Do iSym = 1, nSym
            nAufb(1)=nAufb(1)+nFro(iSym)
          if(nD==2) nAufb(2)=nAufb(2)+nFro(iSym)
         End Do
         Call ICopy(nSym,[0],0,nFro,1)
         call WarningMessage(2, 'Input error!;Aufbau not allowed with frozen orbitals')
         Call Abend()
      End If
!
!---- Check parameters for semi direct SCF
!
!
!---- If semi-direct and I/O buffer value not specified set to
!     default value.
      If (nCore.eq.0.and.nDisc.ne.0) nCore=lDaRec*nSect*2*8/1024
      nCore=((nCore+7)/8)*8
!
!---- Adjust disk size to multiple of I/O buffer
      If (nDisc.ne.0) nDisc=(nDisc*1024+nCore-1)/1024
!
!     nDisc=Min(Int(DiskMx_MByte),nDisc)
      nDisc=Min(10*Allocdisk(),nDisc)
!
!---- Set up parameters that follow from others
!
      If (.Not.Diis .and. .Not.Damping) iDKeep = 0
!
      If (InVec.eq.3 .and. Max(nIter(0),nIter(1)).eq.0) Then
         iPrOrb = 0
         jVOut  = 0
      End If
!
!---- Check Cholesky vs. Aufbau
      If (Aufb .and. DoCholesky) Then
         Cho_Aufb=.true.
!         Write (6,*)
!         Write (6,*) ' ********** WARNING! *********'
!         Write(6,*)' Cholesky SCF runs much slower with AufBau !'
!         Write(6,*)' *** Do you really need AufBau in this case? ***'
!         Write (6,*)
      EndIf
!
!---- Check CONS vs. UHF+OCCU
      If (MxConstr.gt.0 .and. (nD-1+iOCCU).ne.2) Then
         call WarningMessage(2,'For CONStraints, keywords UHF and OCCUpied are compulsory!')
         Call Abend()
      EndIf
!
!---- Check CONS vs. ADDC
      If (MxConstr.eq.0 .and. Do_Addc  ) Then
         call WarningMessage(0,' In absence of CONStraints, ADDCorrelation is de-activated!')
         Do_Addc  =.false.
      EndIf
!
!---- Check CONS vs. SAVE
      If (MxConstr.eq.0 .and. Do_SpinAV) Then
         call WarningMessage(0,' In absence of CONStraints, SAVErage is de-activated!')
         Do_SpinAV=.false.
      EndIf
!
      If (Do_SpinAV) DThr=1.0d4 ! reset because it is not meaningful
      If (MxConstr.gt.0) InVec=6
!
!---- Check parameters of KS-DFT
!
      If (KSDFT.ne.'SCF') Then
         If (MiniDn) Then
            If (jPrint.ge.2) Then
            call WarningMessage(0,' Minimized-density-differences option turned off!')
            End If
            MiniDn = .False.
         End If
         If (Do_OFemb .and. dFMD.ne.Zero) Then
            call WarningMessage(0,' KSDFT/OFE requires DFMD=0 for correlation potential!')
            write(6,*) ' dFMD = ',dFMD
         EndIf
      Else
         If (Do_OFemb .and. dFMD.ne.One) Then
            call WarningMessage(0,' HF/OFE may require DFMD=1 for correlation potential!')
            write(6,*) ' dFMD = ',dFMD
         EndIf
      End If
!
      Call Put_iScalar('SCF mode',nD-1)
!
      LKon = ALGO.eq.4
!
      Method='RHF-SCF '
      If (nD==2) Method='UHF-SCF '
      If (kIvo.ne.0) Method='IVO-SCF '
      If (KSDFT.ne.'SCF') Method='KS-DFT  '
      Call Put_cArray('Relax Method',Method,8)
!
      ExFac=Get_ExFac(KSDFT)
      If (ExFac.eq.Zero .and. .not.TDen_UsrDef   &
                         .and. .not.Do_OFemb ) Then
         DDnOFF = .false. ! use always differential density
      EndIf
!     DDnOFF = .True.
!     MiniDn = .False.
!                                                                      *
!***********************************************************************
!                                                                      *
!     remove local copy of standard input
!
      Call Close_LuSpool(LuSpool)
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
!
!---- Error exit
  902 Continue
      call WarningMessage(2, 'Input error!;Error reading input file for OCCNO option')
      Call Abend()
  903 Continue
      call WarningMessage(2, 'Input error!;End of input file for OCCNO option')
      Call Abend()
!
      End SubRoutine RdInp_scf
