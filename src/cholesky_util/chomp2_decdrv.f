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
! Copyright (C) 2004,2008, Thomas Bondo Pedersen                       *
!***********************************************************************
      SubRoutine ChoMP2_DecDrv(irc,DelOrig,Diag,CD_Type)
!
!     Thomas Bondo Pedersen, Dec. 2004.
!     * Amplitude extension, Thomas Bondo Pedersen, Dec. 2007/Jan. 2008.
!
!     Purpose: decompose (ai|bj) integrals or
!              amplitudes [(ai|bj)/e(a)-e(i)+e(b)-e(j)]
!              for use in Cholesky MP2 program.
!
!     Arguments:
!     irc...... OUT: return code - if non-zero, decomposition failed.
!               Caller MUST check this!
!     DelOrig.. INP: flag for deleting files with original vectors after
!               decomposition completes.
!     Diag..... INP: integral diagonal (ai|ai) or
!               amplitude diagonal (ai|ai)/2[e(a)-e(i)]
!     CD_Type.. INP: string
!               'Integrals'  - integral decomposition
!               'Amplitudes' - amplitude decomposition
!
!     Other input such as orbital energies are read from mbpt2 include
!     files.
!
!
      use ChoMP2, only: OldVec
      use ChoMP2_dec, only: Incore, NowSym, iOption_MP2CD
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None

      Integer  irc
      Logical  DelOrig
      Real*8   Diag(*)
      Character(LEN=*) CD_Type

      External ChoMP2_Col, ChoMP2_Vec
      Integer :: IOPTION, ISYM, LERRSTAT, nBin, kOffD, nDim, iBin, iTyp,
     &           MxQual, LEFT, lB, nInc, lTot, iOpt, iAdr
      Real*8 :: THR, XMN, XMX, RMS

#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"

      Character(LEN=6), Parameter:: ThisNm = 'DecDrv'
      Character(LEN=13), Parameter:: SecNam = 'ChoMP2_DecDrv'

      Logical, Parameter:: Restart = .false.
      Logical Failed, ConventionalCD

      Integer, Parameter :: nOption = 2

      Character(LEN=18) Option

      Integer iClos(2)
      Integer MxCDVec(8)

      Real*8, Allocatable:: ErrStat(:), Bin(:), Qual(:), Buf(:)
      Integer, Allocatable:: iQual(:), iPivot(:)

!     Initializations.
!     ----------------

      irc = 0

#if defined (_DEBUGPRINT_)
      ChkDecoMP2 = .True.
#endif

      If (len(CD_Type) .lt. 1) Then ! input error
         iOption = 0
         Option = 'Unknown !?!?!     '
      Else
         If (CD_Type(1:1).eq.'i' .or. CD_Type(1:1).eq.'I') Then
            iOption = 1
            Option = '(ai|bj) integrals '
         Else If (CD_Type(1:1).eq.'a' .or. CD_Type(1:1).eq.'A') Then
            iOption = 2
            Option = 'MP2 amplitudes    '
         Else
            iOption = nOption + 1
            Option = 'Unknown !?!?!     '
         End If
      End If
      If (iOption.lt.1 .or. iOption.gt.nOption) Then
         irc = -98
         Write(6,*) SecNam,': illegal input option (argument CD_Type)'
         Return
      End If
      iOption_MP2CD = iOption  ! copy to include file chomp2_dec.fh
!-TBP:
! Frankie,
! I use the array MxCDVec(iSym) to decide whether the decomposition
! of a given symmetry block is conventional or "MaxVec":
! ConventionalCD = MxCDVec(iSym).lt.1 .or. MxCDVec(iSym).ge.nT1Am(iSym)
! You may want to calculate the MxCDVec values somewhere else and store
! them in an include-file (say, chomp2_dec.fh).
! For now, I simply define the array here:
      Do iSym = 1,nSym
         MxCDVec(iSym) = nT1Am(iSym)
      End Do

      lErrStat = 3
      Call mma_allocate(ErrStat,lErrStat,Label='ErrStat')
      If (Verbose) Then
         nBin = 18
         Call mma_allocate(Bin,nBin,Label='Bin')
      Else
         nBin  = 0
         Call mma_allocate(Bin,1,Label='Bin')
      End If

      Do iSym = 1,nSym
         nMP2Vec(iSym) = 0
         InCore(iSym)  = .false.
      End Do
      NowSym    = -999999

      If (DelOrig) Then
         iClos(1) = 3  ! signals close and delete original vectors
      Else
         iClos(1) = 2  ! signals close and keep original vectors
      End If
      iClos(2) = 2     ! signals close and keep result vectors

!     Print.
!     ------

      If (Verbose) Then
         Write(6,*)
         Call Cho_Head('Cholesky decomposition of '//Option,'=',80,6)
         Write(6,'(/,1X,A)') 'Configuration of decomposition:'
         Write(6,'(1X,A,1P,D15.6)') 'Threshold: ',ThrMP2
         Write(6,'(1X,A,1P,D15.6)') 'Span     : ',SpanMP2
         If (ChkDecoMP2) Then
            Write(6,'(1X,A)') 'Full decomposition check activated.'
         End If
      End If

!     Start symmetry loop.
!     --------------------

      kOffD = 1
      Do iSym = 1,nSym

         nDim = nT1am(iSym)
         If (nDim.gt.0 .and. NumCho(iSym).gt.0) Then

            ConventionalCD = MxCDVec(iSym).lt.1 .or.
     &                       MxCDVec(iSym).ge.nDim
            If (Verbose .and. nBin.gt.0) Then
               Bin(1) = 1.0D2
               Do iBin = 2,nBin
                  Bin(iBin) = Bin(iBin-1)*1.0D-1
               End Do
               If (ConventionalCD) Then
                  Write(6,'(//,1X,A,I2,A,I9)')
     &     '>>> Conventional Cholesky decomposition of symmetry block ',
     &            iSym,
     &            ', dimension: ',nDim
               Else
                  Write(6,'(//,1X,A,I2,A,I9)')
     &           '>>> MaxVec Cholesky decomposition of symmetry block ',
     &            iSym,
     &            ', dimension: ',nDim
               End If
               Write(6,'(/,1X,A)') 'Analysis of initial diagonal:'
               Call Cho_AnaSize(Diag(kOffD),nDim,Bin,nBin,6)
            End If

!           Open files.
!           -----------

            Do iTyp = 1,2
               Call ChoMP2_OpenF(1,iTyp,iSym)
            End Do

!           Setup decomposition.
!           --------------------

            NowSym = iSym

            If (MxQualMP2 .ne. MxQual_Def) Then ! user-defined
               MxQual = min(max(MxQualMP2,1),nDim)
            Else ! default
               If (nDim .gt. 10) Then
                  MxQual = max(min(nDim/10,MxQualMP2),1)
               Else
                  MxQual = max(min(nDim,MxQualMP2),1)
               End If
            End If
#if !defined (_I8_)
            lTstBuf = (nDim+MxQual)*MxQual
            lTstQua = nDim*(MxQual+1)
            Do While ((lTstBuf.lt.0.or.lTstQua.lt.0) .and. MxQual.gt.0)
               MxQual  = MxQual - 1
               lTstBuf = (nDim+MxQual)*MxQual
               lTstQua = nDim*(MxQual+1)
            End Do
            If (MxQual .lt. 1) Then
               Write(6,*) SecNam,': MxQual causes integer overflow!'
               Write(6,*) SecNam,': parameters:'
               Write(6,*) 'Symmetry block: ',iSym
               Write(6,*) 'Dimension     : ',nDim
               Write(6,*) 'MxQual        : ',MxQual
               irc = -99
               Go To 1 ! exit
            End If
#endif

            Call mma_allocate(Qual,nDim*(MxQual + 1),Label='Qual')
            Call mma_allocate(iQual,MxQual,Label='iQual')
            Call mma_allocate(iPivot,nDim,Label='iPivot')

            Call mma_maxDBLE(lB)
            lBuf = min((nDim+MxQual)*MxQual,lB)
            Left = lB - lBuf
            nInC = Left/nDim
            If (nInC .ge. NumCho(iSym)) Then
               InCore(iSym) = .true.
               lTot = nDim*NumCho(iSym)
               Call mma_allocate(OldVec,lTot,Label='OldVec')
               iOpt = 2
               iAdr = 1
               Call ddaFile(lUnit_F(iSym,1),iOpt,OldVec,lTot,iAdr)
            End If
            Call mma_maxDBLE(lBuf)
            Call mma_allocate(Buf,lBuf,Label='Buf')

!           Decompose this symmetry block.
!           ------------------------------

            Thr  = ThrMP2
            Span = SpanMP2
            If (ConventionalCD) Then
               Call ChoDec(ChoMP2_Col,ChoMP2_Vec,
     &                     Restart,Thr,Span,MxQual,
     &                     Diag(kOffD),Qual,Buf,
     &                     iPivot,iQual,
     &                     nDim,lBuf,
     &                     ErrStat,nMP2Vec(iSym),irc)
               If (irc .ne. 0) Then
                  Write(6,*) SecNam,': ChoDec returned ',irc,
     &                              '   Symmetry block: ',iSym
                  Go To 1 ! exit...
               End If
            Else
               Call ChoDec_MxVec(ChoMP2_Col,ChoMP2_Vec,MxCDVec(iSym),
     &                           Restart,Thr,Span,MxQual,
     &                           Diag(kOffD),Qual,Buf,
     &                           iPivot,iQual,
     &                           nDim,lBuf,
     &                           ErrStat,nMP2Vec(iSym),irc)
               If (irc .ne. 0) Then
                  Write(6,*) SecNam,': ChoDec_MxVec returned ',irc,
     &                              '   Symmetry block: ',iSym
                  Go To 1 ! exit...
               End If
            End If
            XMn = ErrStat(1)
            XMx = ErrStat(2)
            RMS = ErrStat(3)
            If (Verbose) Then
               Write(6,'(/,1X,A)')
     &         '- decomposition completed!'
               Write(6,'(1X,A,I9,A,I9,A)')
     &         'Number of vectors needed: ',nMP2Vec(iSym),
     &         ' (number of AO vectors: ',NumCho(iSym),')'
               If (.not. ConventionalCD) Then
                  Write(6,'(1X,A,I9)')
     &            'Max. number of vectors allowed: ',MxCDVec(iSym)
               End If
               Write(6,'(1X,A)')
     &         'Error statistics for diagonal [min,max,rms]:'
               Write(6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
            End If
            If (ConventionalCD) Then
               Failed = abs(Xmn).gt.Thr .or. abs(XMx).gt.thr .or.
     &                  RMS.gt.Thr
            Else
               Failed = abs(Xmn).gt.Thr
            End If
            If (Failed) Then
               If (.not. Verbose) Then
                  Write(6,'(1X,A)')
     &            'Error statistics for diagonal [min,max,rms]:'
                  Write(6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
               End If
               Write(6,'(A,A,A,A)')
     &         SecNam,': decomposition of ',Option,' failed!'
               irc = -9999
               Go To 1 ! exit
            End If

!           If requested, check decomposition.
!           ----------------------------------

            If (ChkDecoMP2) Then
               Write(6,*)
               Write(6,'(A,A,A)')
     &         SecNam,': Checking decomposition of ',Option
               Write(6,*) 'Symmetry block: ',iSym
               Write(6,*) 'Threshold, Span, MxQual: ',Thr,Span,MxQual
               Write(6,*) 'Error statistics for diagonal [min,max,rms]:'
               Write(6,*)  ErrStat(:)
               Call ChoMP2_DecChk(irc,iSym,Qual,nDim,MxQual,
     &                            Buf,lBuf,ErrStat)
               If (irc .ne. 0) Then
                  If (irc .eq. -123456) Then
                     Write(6,*)
     &                       ' -- Sorry, full decomposition check not ',
     &                       'yet implemented --'
                     irc = 0
                  Else
                     Write(6,*) SecNam,': ChoMP2_DecChk returned ',irc,
     &                                 '   Symmetry block: ',iSym
                     Call ChoMP2_Quit(SecNam,'decomposition failed!',
     &                                ' ')
                  End If
               Else
                  XMn = ErrStat(1)
                  XMx = ErrStat(2)
                  RMS = ErrStat(3)
                  Write(6,'(A,A,A)')
     &            'Error statistics for ',Option,' [min,max,rms]:'
                  Write(6,*) XMn,XMx,RMS
                  Failed = Failed .or. abs(Xmn).gt.Thr .or.
     &                     abs(XMx).gt.Thr .or. RMS.gt.Thr
                  If (ConventionalCD) Then
                     If (Failed) Then
                        Write(6,*) '==> DECOMPOSITION FAILURE <=='
                        irc = -9999
                        Go To 1 ! exit
                     Else
                        Write(6,*) '==> DECOMPOSITION SUCCESS <=='
                     End If
                  Else
                     If (Failed) Then
                        Write(6,*)
     &                  '==> DECOMPOSITION SUCCESS <== (by definition)'
                     Else
                        Write(6,*) '==> DECOMPOSITION SUCCESS <=='
                     End If
                  End If
                  Call xFlush(6)
               End If
            End If

!           Free memory.
!           ------------

            Call mma_deallocate(Buf)
            If (InCore(iSym)) Call mma_deallocate(OldVec)
            Call mma_deallocate(iPivot)
            Call mma_deallocate(iQual)
            Call mma_deallocate(Qual)

!           Close (possibly deleting original) files.
!           -----------------------------------------

            Do iTyp = 1,2
               Call ChoMP2_OpenF(iClos(iTyp),iTyp,iSym)
            End Do

!           Update pointer to diagonal block.
!           ---------------------------------

            kOffD = kOffD + nT1am(iSym)

         Else

            If (Verbose) Then
               Write(6,'(//,1X,A,I2,A)')
     &         '>>> Symmetry block',iSym,' is empty!'
            End If

         End If

      End Do

    1 Continue
      If (irc .ne. 0) Then ! make sure files are closed before exit
         Do iSym = 1,nSym
            Do iTyp = 1,2
               Call ChoMP2_OpenF(2,iTyp,iSym)
            End Do
         End Do
      End If

      Call mma_deallocate(Bin)
      Call mma_deallocate(ErrStat)
      End
