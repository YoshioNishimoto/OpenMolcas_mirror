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
!               1995, Martin Schuetz                                   *
!***********************************************************************
      SubRoutine WrInp_SCF(SIntTh)
!***********************************************************************
!                                                                      *
!     purpose: Write input                                             *
!                                                                      *
!     input: SIntTh: Threshold for Integral prescreening               *
!                                                                      *
!***********************************************************************
      Use Functionals, only: Print_Info
      Use KSDFT_Info, only: CoefR, CoefX
      Use InfSO, only: DltNth, QNRTh, IterSO_Max
      use InfSCF, only: Aufb, DDnoff, DelThr, DIIS, DIISTh, DoCholesky, DSCF, DThr, EThr, FThr, iAU_ab,        &
                        InVec, isHDF5, nD, jPrint, jVOut, kIVO,kOptim_Max, KSDFT, LKOn, lpaper, MiniDn, nCore,     &
                        nDIsc, nMem, NoExchange, nSym, nTit, One_Grid,PreSch, RFPert, rTemp, Scrmbl, StVec, Teee,    &
                        TemFac, Thize, Tot_Charge, Tot_El_Charge,Tot_Nuc_Charge, TStop, VTitle, Header, Title,       &
                        nFro, nAufb, nOcc, nOrb, nBas, nIter, nDel
      use ChoSCF, only: dmpk, Algo, ReOrd
      use Fock_util_global, only: Deco
      use RICD_Info, only: Do_DCCD
!
      Implicit None
      Real*8 SIntTh

!---- Define local variables
      Integer i, iCharge, iDoRI, iSym, iTit, mTmp
      Character(LEN=60) Fmt, FmtR, FmtI
      Character(LEN=72) Line
      Character(LEN=3) lIrrep(8)
      Logical NonEq

      If (jPrint.ge.2) Then
         Call CollapseOutput(1,'   Input section:')
         Write(6,'(3X,A)')     '   --------------'
         Write(6,*)
      End If
!
!---- Print out header of the integral file
      If (jPrint.ge.2) Then
         Write(6,'(6X,A)')  'Header of the integral files:'
         Write(Line,'(A72)') Header(1)
         Write(6,'(6X,A)') trim(Line)
         Write(Line,'(A72)') Header(2)
         Write(6,'(6X,A)') trim(Line)
         Write(6,*)
      End If
!
!---- Print out title (if any)
      If (nTit.gt.0) Then
         If (jPrint.ge.3) Then
             Call Banner(Title,nTit,lPaper-7)
             Write(6,*)
         Else If (jPrint.ge.2) Then
             Write (6,'(6X,A)') ' Title:'
             Do iTit = 1, nTit
                Write (6,'(8X,A)') trim(Title(iTit))
             End Do
         End If
      End If
!
!---- Print the coordinates of the system
      If (jPrint.ge.2) Call PrCoor
!
!---- Print out Print level
!     Write(*,'(1X,A,I3)') 'Print level=',iPrint
!     Write(*,*)
!
!
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         lIrrep(iSym) = adjustr(lIrrep(iSym))
      End Do
!
      If (jPrint.ge.2) Then
        Call CollapseOutput(0,'   Input section:')
        Write(6,*)
      End If
!---- Print out orbital informations
      Fmt='(6X,A,T35,8I4)'
      If (jPrint.ge.2) Then
         Call CollapseOutput(1,'   Orbital specifications:')
         Write(6,'(3X,A)')     '   -----------------------'
         Write(6,*)
         Write(6,Fmt)'Symmetry species',         (i,i=1,nSym)
         Write(6,'(6X,A,T35,8(1X,A))') '                ',         (lIrrep(i),i=1,nSym)
         Write(6,Fmt)'Frozen orbitals',          (nFro(i),i=1,nSym)
      End If
!
      If (Aufb) Then
       if(nD==1) Then
        If (nAufb(1).eq.-1) Then
           Tot_El_Charge=Tot_Charge-Tot_Nuc_Charge
!--------- Check that Tot_El_Charge is a multiple of two!
           mtmp=Int(-Tot_El_Charge+0.1D0)
           If (Mod(mtmp,2).ne.0) Then
              Write (6,*) 'WrInp: Error in definition of molecular charge!'
              Write (6,*) 'Current implementation only allows double occupations.'
              Write (6,*) 'Tot_Charge    :',Tot_Charge
              Write (6,*) 'Tot_El_Charge :',Tot_El_Charge
              Write (6,*) 'Tot_Nuc_Charge:',Tot_Nuc_Charge
              Call Abend()
           End If
           nAufb(1)=mtmp/2
        endif
       else
        If (nAufb(1).eq.-1) Then
           Tot_El_Charge=Tot_Charge-Tot_Nuc_Charge
!--------- Check that Tot_El_Charge is a multiple of two!
           mtmp=Int(-Tot_El_Charge+0.1D0)
! if ZSPIN is not set - make difference alpha-beta = 0 or 1
            nAufb(2)=(mtmp-iAu_ab)/2
            nAufb(1)=Int(-Tot_El_Charge-nAufb(2))
        endif
!          Write(6,*) ' CHARGE + UHF is un'
!           Call Abend()
        End If
        If (nD==1.and.jPrint.ge.2) then
           Write(6,Fmt)'Aufbau',                 nAufb(1)
        else if (jPrint.ge.3) Then
           Write(6,Fmt)'Aufbau alpha',                 nAufb(1)
           Write(6,Fmt)'Aufbau beta ',                 nAufb(2)
        endif
        If (Teee.and.jPrint.ge.2) Then
           Write (6,'(a,f6.3)') '      Start temperature =',RTemp
           Write (6,'(a,f6.3)') '      End temperature   =',TStop
           Write (6,'(a,f6.3)') '      Temperature Factor=',TemFac
        End If
      Else
        if(nD==1.and.jPrint.ge.2) then
        Write(6,Fmt)'Occupied orbitals',    (nOcc(i,1),i=1,nSym)
        Write(6,Fmt)'Secondary orbitals',   (nOrb(i)-nOcc(i,1),i=1,nSym)
        else if (jPrint.ge.2) Then
        Write(6,Fmt)'Occupied orbitals alpha',      (nOcc(i,1),i=1,nSym)
        Write(6,Fmt)'Occupied orbitals beta ',      (nOcc(i,2),i=1,nSym)
        Write(6,Fmt)'Secondary orbitals alpha', (nOrb(i)-nOcc(i,1),i=1,nSym)
        Write(6,Fmt)'Secondary orbitals beta',  (nOrb(i)-nOcc(i,2),i=1,nSym)

        endif
      End If
!
      Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
      Call put_dscalar('Total Charge    ',Tot_Charge)
      If (jPrint.ge.2) Then
         Write(6,Fmt)'Deleted orbitals         ',(nDel(i),i=1,nSym)
         Write(6,Fmt)'Total number of orbitals ',(nOrb(i),i=1,nSym)
         Write(6,Fmt)'Number of basis functions',(nBas(i),i=1,nSym)
         Call CollapseOutput(0,'   Orbital specifications:')
         Write(6,*)
         Write(6,'(6X,A,T45,F10.3)') 'Molecular charge',Tot_Charge
         Write(6,*)
      End If
!
!---- Print out reaction field specifications
!
      iCharge=Int(Tot_Charge)
      NonEq=.False.
      Call PrRF(.False.,NonEq,iCharge,jPrint)
      If ( RFpert.and.jPrint.ge.2 ) then
         Write(6,'(6X,A)')'Reaction field specifications:'
         Write(6,'(6X,A)')'------------------------------'
         Write(6,*)
         Write(6,'(6X,A)')'The Reaction field is added as a perturbation and has been determined in a previous calculation'
         Write(6,*)
      End If
!
!---- Print out grid information in case of DFT
!
      If (KSDFT.ne.'SCF') Then
         Call Put_dScalar('DFT exch coeff',CoefX)
         Call Put_dScalar('DFT corr coeff',CoefR)
         Call Put_dScalar('EThr',EThr)
         Call Funi_Print()
         If (jPrint.ge.2) Then
            If (One_Grid) Then
               Write (6,'(6X,A)') 'The same grid will be used for all iterations.'
            Else
               Write (6,'(6X,A)') 'A smaller intermediate grid will be used the first few iterations.'
            End If
            Write(6,*)
            Write(6,'(6X,A)') 'DFT functional specifications'
            Write(6,'(6X,A)') '-----------------------------'
            Call libxc_version()
            Call Print_Info()
            Write(6,*)
         End If
      End If
!
!---- Print out informations concerning Direct/Conventional scheme
      FmtI= '(6X,A,T50,I9)'
      FmtR= '(6X,A,T50,ES9.2)'
      If (jPrint.ge.2) Then
         Call CollapseOutput(1,'   Optimization specifications:')
         Write(6,'(3X,A)')     '   ----------------------------'
         Write(6,*)
      End If
      If (DSCF.and..NOT.Do_DCCD.and.jPrint.ge.2) Then
         If (nDisc.eq.0.and.nCore.eq.0) Then
            Write(6,'(6X,A)')'SCF Algorithm: Direct'
         Else If (nDisc*1024.ge.nCore) Then
            Write(6,'(6X,A)')'SCF Algorithm: Semi-direct'
!
!---------- The threshold to be used in the direct SCF procedure is
!           defined by the energy threshold.
            Write (6,FmtI) 'Max MByte of integrals on disk/process:',nDisc
            Write(6,FmtR) 'Threshold for saving integrals on disc', Thize
         Else
            Write(6,'(6X,A)')'SCF Algorithm: Semi-direct in-core'
!
!---------- The threshold to be used in the direct SCF procedure is
!           defined by the energy threshold.
            Write (6,FmtI) 'Max kByte of integrals in memory/process:',nCore
            Write(6,FmtR) 'Threshold for saving integrals in memory', Thize
         End If
         If (PreSch) Then
            Write(6,'(6X,A)')'Prescreening Scheme: Only Integral value'
         Else
            Write(6,'(6X,A)')'Prescreening Scheme: Integral*Density value'
         End If
      Else If (jPrint.ge.2) Then
         if(nD==1) then
            if(.not.DoCholesky)then
              Write(6,'(6X,A)')'SCF Algorithm: Conventional'
            else
              Call Get_iScalar('System BitSwitch',iDoRI)
              if (Iand(iDoRI,1024).Eq.1024) then
                 if (LKon) then
                    Write(6,'(6X,A)')'SCF Algorithm: LK-RI/DF'
                    Write(6,FmtR) 'LK screening threshold:',dmpk
                 else
                    Write(6,'(6X,A)')'SCF Algorithm: RI/DF'
                 endif
              else
                 if (LKon) then
                    Write(6,'(6X,A)')'SCF Algorithm: LK-Cholesky'
                    Write(6,FmtR) 'LK screening threshold:',dmpk
                 else
                    Write(6,'(6X,A,I1)') 'SCF Algorithm: Cholesky'
                 endif
              endif

             if(ALGO.eq.0)then
              if (Iand(iDoRI,1024).Eq.1024) then
               Write(6,'(6X,A)') 'Integral regeneration from RI vectors reordered on disk'
              else
               Write(6,'(6X,A)') 'Integral regeneration from Cholesky vectors reordered on disk'
              endif
             elseif(ALGO.eq.1)then
               Write(6,'(6X,A)') 'Density-based Cholesky. Default reorder: on the fly'
             elseif(ALGO.eq.2)then
               Write(6,'(6X,A)') 'MO-based-Exchange Cholesky. Default reorder: on the fly'
             elseif(ALGO.eq.3)then
               Write(6,'(6X,A)') 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
             elseif(ALGO.eq.4)then
               Write(6,'(6X,A)') 'Local-Exchange (LK) algorithm.'
             endif

             If (Do_DCCD) Write (6,'(6X,A)') ' - Corrected with exact 1-center two-electron integrals'
             If (ReOrd) Write (6,'(6X,A)')   ' - the Cholesky vectors are reordered'
             If (DeCo) Write (6,'(6X,A)')    ' - the density matrix is decomposed'
             If (DSCF.and.nDisc==0.and.nCore==0) Write (6,'(6X,A)') ' - SCF Algorithm: Direct'
            endif
         else
            if(.not.DoCholesky)then
              Write(6,'(6X,A)')'SCF Algorithm: Conventional USCF'
            else
              Call Get_iScalar('System BitSwitch',iDoRI)
              if (Iand(iDoRI,1024).Eq.1024) then
                 if (LKon) then
                    Write(6,'(6X,A)')'SCF Algorithm: LK-RI/DF USCF'
                    Write(6,FmtR) 'LK screening threshold:',dmpk
                 else
                    Write(6,'(6X,A)')'SCF Algorithm: RI/DF USCF'
                 endif
              else
                 if (LKon) then
                    Write(6,'(6X,A)') 'SCF Algorithm: LK-Cholesky USCF'
                    Write(6,FmtR) 'LK screening threshold:',dmpk
                 else
                    Write(6,'(6X,A)')'SCF Algorithm: Cholesky USCF'
                 endif
              endif

              if(ALGO.eq.0)then
               if (Iand(iDoRI,1024).Eq.1024) then
                Write(6,'(6X,A)') 'Integral regeneration from RI vectors reordered on disk'
               else
                Write(6,'(6X,A)') 'Integral regeneration from Cholesky vectors reordered on disk'
               endif
              elseif(ALGO.eq.1)then
                Write(6,'(6X,A)') 'Density-based Cholesky. Default reorder: on the fly'
              elseif(ALGO.eq.2)then
                Write(6,'(6X,A)') 'MO-based-Exchange Cholesky. Default reorder: on the fly'
              elseif(ALGO.eq.3)then
                Write(6,'(6X,A)') 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
              elseif(ALGO.eq.4)then
                Write(6,'(6X,A)') 'Local-Exchange (LK) algorithm.'
              endif

              If (Do_DCCD) Write (6,'(6X,A)') ' - Corrected with exact 1-center two-electron integrals'
              If (ReOrd) Write (6,'(6X,A)')   ' - the Cholesky vectors are reordered'
              If (DeCo) Write (6,'(6X,A)')    ' - the density matrix is decomposed'
            endif
         endif
      End If
      If (NoExchange) Then
         Write(6,'(6X,A)') 'NOTE: exchange contributions will not be computed'
      End If
!
      If (jPrint.ge.2) Then
!
!---- Print out informations concerning difference scheme used
      If (MiniDn) Then
         Write(6,'(6X,A)')'Minimized density differences are used'
      Else
        if(.not.DDnOFF)then
          Write(6,'(6X,A)')'D(i)-D(i-1) density differences are used'
        else
          Write(6,'(6X,A)')'The actual AO density is used'
        endif
      End If
      Write(6,FmtI)'Number of density matrices in core',nMem
!
!---- Print out number of iterations
      Write(6,FmtI) 'Maximum number of NDDO SCF iterations',nIter(0)
      Write(6,FmtI) 'Maximum number of HF SCF iterations',nIter(1)
!
!---- Print out thresholds for SCF
      Write(6,FmtR) 'Threshold for SCF energy change',EThr
      Write(6,FmtR) 'Threshold for density matrix',DThr
      Write(6,FmtR) 'Threshold for Fock matrix',FThr
      Write(6,FmtR) 'Threshold for linear dependence',DelThr
      If (Diis) Then
         Write(6,FmtR) 'Threshold at which DIIS is turned on', DiisTh
         Write(6,FmtI) 'Maximum depth in the DIIS procedure',kOptim_Max
         Write(6,FmtI) 'Maximum depth in the BFGS Hessian update', IterSO_Max
         Write(6,FmtR) 'Threshold at which QNR/C2DIIS is turned on', QNRTh
         Write(6,FmtR) 'Threshold for Norm(delta) (QNR/C2DIIS)', DltNTh
      End If
      If (DSCF) Then
         Write(6,FmtR) 'Threshold for contribution to Fock matrix', SIntTh
      End If
!
      Fmt = '(6x,A,A)'
!
!---- Print out IVO information (if any)
      If (kIvo.ne.0) Write(6,Fmt) 'Improved virtual orbitals.'
!
!---- Print out information about orbitals punched on the file
      If (jVOut.le.0) Then
         Write(6,Fmt) 'No vectors punched'
      Else If (jVOut.eq.1) Then
         If(nD==1) Then
            Write(6,Fmt) 'All non deleted orbitals punched on: SCFORB'
         Else
            Write(6,Fmt) 'All non deleted orbitals punched on: UHFORB'
         End If
      Else
         If(nD==1) Then
            Write(6,Fmt) 'All orbitals punched on: SCFORB'
         Else
            Write(6,Fmt) 'All orbitals punched on: UHFORB'
         End If
      End If
      Call CollapseOutput(0,'   Optimization specifications:')
      Write(6,*)
!
!---- Print out
      If (InVec.eq.0) Then
         Write(6,Fmt) 'Starting vectors from core diagonalization'
      Else If (InVec.eq.1) Then
         Write(6,Fmt) 'NDDO MOs are generated before actual HF SCF computation'
      Else If (InVec.eq.2) Then
         If (isHDF5) Then
            Write(6,Fmt) 'Input vectors read from HDF5 file'
         Else
            Write(6,Fmt) 'Input vectors read from INPORB'
            Write(6,Fmt) 'Orbital file label: ',trim(VTitle)
         End If
      Else If (InVec.eq.3) Then
         Write(6,Fmt) 'Input density matrix read from RUNFILE'
      Else If (InVec.eq.4) Then
         Write(6,Fmt) 'Restart...'
      Else If (InVec.eq.5) Then
         Write(6,Fmt) 'Input vectors from NDDO calculation'
      Else
         Write(6,Fmt) StVec
      End If
      If (Scrmbl) Write(6,Fmt) 'Start orbitals are scrambled in order to introduce symmetry breaking'
      Write(6,*)
!
      End If
!

      Call XFlush(6)
      End
