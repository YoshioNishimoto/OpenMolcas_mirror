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
* Copyright (C) 1990, Markus P. Fuelscher                              *
*               2013, Giovanni Li Manni                                *
*               2016, Andrew M. Sand                                   *
************************************************************************
      Subroutine MSCtl(CMO,FI,FA,Ref_Ener)
!CMO,F,FI
*
*     This routine is a modification of SGFCIN, adapted to a CASDFT
*     implementation in which the CI step of a CASDFT calculation is
*         not corrected by DFT. DFT will play a role only in the Orbital
*         optimization step.
*     Purpose:
*     Generate the Fock-matrix for the frozen and inactive orbitals.
*     Compute also the core energy. Finally, transform the Fock-matrix
*     into the basis of the active orbitals.
*
*     M.P. Fuelscher, Lund, July 1990
*     GLM, Minneapolis,   May 2013
*     AMS, Minneapolis,   Feb 2016
*
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      use OneDat, only: sNoNuc, sNoOri
      Use KSDFT_Info, only: do_pdftpot, ifav, ifiv
      Use hybridpdft, only: Do_Hybrid, E_NoHyb, Ratio_WF
      use mspdft, only: dogradmspd, do_rotate, iIntS, iDIDA, IP2MOt,
     &                  D1AOMS, D1SAOMS
      use mcpdft_output, only: debug, lf, iPrLoc

      Implicit Real*8 (A-H,O-Z)

      Dimension CMO(*), FI(*), FA(*), Ref_Ener(*)
*
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "pamint.fh"
#include "timers.fh"
#include "SysDef.fh"
#include "gugx.fh"
#include "wadr.fh"
#include "rasscf_lucia.fh"


      Character*8 Label
      Logical First, Dff, Do_DFT,Found
      Parameter ( Zero=0.0d0 , One=1.0d0 )
      integer iD1I,iD1Act,iD1ActAO,iD1Spin,iD1SpinAO,IAD19
      integer iJOB,dmDisk,iP2d
      integer itmp1,itmp2,itmp3,itmp4
      integer itmp5,itmp6,itmp7,itmpn,itmpk,itmpa
      integer ifocki,ifocka
      integer, dimension(30) :: IADR19
      integer LP,NQ,LQ,LPUVX
      integer LOEOTP
      integer jroot
      real(kind=wp), dimension(nroots) :: Energies
      integer  i_off1,i_off2,ifone
      integer isym
      integer LUGS
      real(kind=wp), allocatable, dimension(:) :: P2dt1, tmp_grd, scr1
      External IsFreeUnit

***********************************************************
C Local print level (if any)
***********************************************************
      IPRLEV=IPRLOC(3)

! Compute size of Q matrix...not sure if we really need it any more though
      NQ=0
      NIAIA=0
      do ISYM=1,NSYM
        NQ = MAX(NQ,NASH(ISYM)*NORB(ISYM))
        NIAIA = NIAIA+(NASH(ISYM)+NISH(ISYM))**2
      end do
      if(NQ < NIAIA) NQ=NIAIA


*TRS
      Call Get_iScalar('Relax CASSCF root',iRlxRoot)
***********************************************************
* Load the nuclear repulsion energy
***********************************************************
*TRS
*
      Call Get_dScalar('PotNuc',potNuc)

***********************************************************
* Load bare nuclei Hamiltonian
! This is h_pq but in the AO basis (so h_{mu, nu})
***********************************************************
      Call GetMem('Fcore','Allo','Real',iTmp1,nTot1)
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp1),iSyLbl)
      If ( iRc.ne.0 ) then
        Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
        Write(LF,*) 'Label = ',Label
        Write(LF,*) 'iRc = ',iRc
        Call Abend
      Endif
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' OneHam in AO basis in CASDFT_Terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=0
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(iTmp1+iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

      Call GetMem('Kincore','Allo','Real',iTmpk,nTot1)
c--reads kinetic energy integrals  Work(iTmpk)--(Label=Kinetic)----
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'Kinetic '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmpk),iSyLbl)
      If ( iRc.ne.0 ) then
        Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
        Write(LF,*) 'Label = ',Label
        Write(LF,*) 'iRc = ',iRc
        Call Abend
      Endif
      Call GetMem('NucElcore','Allo','Real',iTmpn,nTot1)
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'Attract '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmpn),iSyLbl)
      If ( iRc.ne.0 ) then
        Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
        Write(LF,*) 'Label = ',Label
        Write(LF,*) 'iRc = ',iRc
        Call Abend
      Endif

!Here we calculate the D1 Inactive matrix (AO).
      Call GetMem('D1Inact','Allo','Real',iD1I,NTOT2)
      Call Get_D1I_mcpdft(CMO,Work(iD1I))

      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'iD1inact'
        do i=1,ntot2
          write(6,*) Work(iD1i-1+i)
        end do
      END IF

      iJOB=0
      IAD19=0
      Call f_Inquire('JOBOLD',Found)
      If (.not.found) then
        Call f_Inquire('JOBIPH',Found)
        if (Found) JOBOLD=JOBIPH
      end if
      If (Found) iJOB=1
      If (iJOB.eq.1) Then
         if(JOBOLD.le.0) Then
           JOBOLD=20
           Call DaName(JOBOLD,'JOBOLD')
         end if
      end if
      IADR19(:)=0
      Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
      IADR15 = IADR19
      dmDisk = IADR19(3)

      Call GetMem('D1Active','Allo','Real',iD1Act,NACPAR)
      Call GetMem('D1ActiveAO','Allo','Real',iD1ActAO,NTOT2)
      Call GetMem('D1Spin','Allo','Real',iD1Spin,NACPAR)
      Call GetMem('D1SpinAO','Allo','Real',iD1SpinAO,NTOT2)

      Call GetMem('DtmpI','Allo','Real',iTmp3,nTot1)
      Call GetMem('DtmpA','Allo','Real',iTmp4,nTot1)
      Call GetMem('P2','Allo','Real',iP2d,NACPR2)


      Call GetMem('FockI','ALLO','Real',ifocki,ntot1)
      Call GetMem('FockA','ALLO','Real',ifocka,ntot1)

************************************************************************
* load back two-electron integrals (pu!vx)
************************************************************************
      lPUVX = 1

      If (nFint > 0) then
        iDisk = 0
        Call GetMem('PUVX','Allo','Real',lPUVX,nFint)
        Call DDaFile(LUINTM,2,Work(lPUVX),nFint,iDisk)
      End If

      IF(IPRLEV >= DEBUG) THEN
        write(6,*) 'PUVX integrals in msctl'
        call wrtmat(Work(lPUVX),1,nFInt,1,nFInt)
      END IF


!This iSA is used to control gradient calculations.  Analytic gradients
!(in ALASKA) will only run if iSA=1, and iSA will only be set to one if
!the on-top potentials are computed as part of this calculation.
      iSA = 99
      Call Put_iScalar('SA ready',iSA)


!We loop over the number of states.  For each state, we read the density
!matrices from the JOBIPH file.
      do jroot=1,lroots
        iIntS=jRoot

       !Load a fresh FockI and FockA
        Call dcopy_(ntot1,FI,1,Work(ifocki),1)
        Call dcopy_(ntot1,FA,1,Work(ifocka),1)
*
!Read in the density matrices for <jroot>.
        Call Fzero(Work(iD1Act),NACPAR)
        Call Fzero(Work(iD1ActAO),NTOT2)
        Call Fzero(Work(iD1Spin),NACPAR)
        Call Fzero(Work(iD1SpinAO),NTOT2)
        Call Fzero(Work(iTmp3),nTot1)
        Call Fzero(Work(iTmp4),nTot1)
        Call Fzero(Work(iP2d),NACPR2)

!Get the D1 Active matrix for this state.  These should probably be
!most easily read from the previous JOBIPH file.  Then, convert D1A from
!the MO to the AO basis.

        Call DDaFile(JOBOLD,2,Work(iD1Act),NACPAR,dmDisk)
        Call DDaFile(JOBOLD,2,Work(iD1Spin),NACPAR,dmDisk)
        Call DDaFile(JOBOLD,2,Work(iP2d),NACPR2,dmDisk)
        Call DDaFile(JOBOLD,0,Work(iP2d),NACPR2,dmDisk)

        if(DoGradPDFT .and. jroot == irlxroot) then
          call mma_allocate(P2dt1, nacpr2, "P2t")
          P2dt1 = 0
          call p2_contraction(work(id1act), p2dt1)
          call put_dArray('P2MOt',p2dt1, nacpr2)
          call mma_deallocate(P2dt1)
        else if(dogradmspd) then
          Call P2_contraction(Work(iD1Act),
     &                        Work(iP2MOt+(jroot-1)*NACPR2))
        end if

        Call Put_dArray('D1mo',Work(iD1Act),NACPAR)
        Call Put_dArray('P2mo',Work(iP2d),NACPR2)

!**********************************************************
! Generate total density
!**********************************************************

        If(NASH(1) /= NAC) Call DBLOCK_m(Work(iD1Act))
        Call cas_mo_to_ao(CMO,Work(iD1Act),Work(iD1ActAO))

        if((DoGradPDFT .and. jroot == irlxroot) .or. DoGradMSPD) then
          call mma_allocate(tmp_grd, ntot1, "dTmpA_g")
          call fold_pdft(nsym, nbas, work(iD1ActAO), tmp_grd)
          call put_darray('d1activeao', tmp_grd, ntot1)
          call mma_deallocate(tmp_grd)
        end if

        Call Fold(nSym,nBas,Work(iD1I),Work(iTmp3))
        Call Fold(nSym,nBas,Work(iD1ActAO),Work(iTmp4))

        if(DoGradMSPD) then
          Call Dcopy_(nTot1,Work(iTmp4),1,Work(iDIDA+(iIntS-1)*nTot1),1)
          if (iIntS == lRoots) then
            Call Dcopy_(ntot1,Work(iTmp3),1,Work(iDIDA+lRoots*nTot1),1)
          end if
        end if
        ! The full 1RDM (inactive/frozen + active) stored in iTmp3
        Call Daxpy_(nTot1,1.0D0,Work(iTmp4),1,Work(iTmp3),1)
!Maybe I can write all of these matrices to file, then modify stuff in
!the nq code to read in the needed density.  In other words, I need to
!replace the next call with something that supports multistates.
        Call Put_dArray('D1ao',Work(iTmp3),nTot1)
        IF(DoGradMSPD)THEN
          Call DCopy_(nTot1,Work(iTmp3),1,
     &                Work(D1AOMS+(jRoot-1)*nTot1),1)
        END IF

!Get the spin density matrix for open shell cases
***********************************************************
* Generate spin-density
***********************************************************
        IF(NASH(1) /= NAC) then
          CALL DBLOCK_m(Work(iD1Spin))
        end if

        call cas_mo_to_ao(CMO,Work(iD1Spin),Work(iD1SpinAO))
        Call GetMem('DtmpS','Allo','Real',iTmp7,nTot1)
        Call Fold(nSym,nBas,Work(iD1SpinAO),Work(iTmp7))
        Call Put_dArray('D1sao',Work(iTmp7),nTot1)

        IF(iSpin /= 1 .and. DoGradMSPD) THEN
          Call DCopy_(nTot1,Work(iTmp7),1,
     &                Work(D1SAOMS+(jRoot-1)*nTot1),1)
        END IF
        Call GetMem('DtmpS','Free','Real',iTmp7,nTot1)

***********************************************************
* Calculation of the CASDFT_energy
***********************************************************
!Perhaps ideally, we should rework how DrvXV (and its children) handles
!the AO to MO transformation on the grid.  It seems like perhaps we are
!doing redundant transformations by retransforming AOs (which may have
!been included in a previous batch) into MOs.
*
        Call GetMem('htmp','Allo','Real',iTmp5,nTot1)
        Call GetMem('gtmp','Allo','Real',iTmp6,nTot1)
        Call dCopy_(nTot1,[0.0d0],0,Work(iTmp5),1)
        Call dCopy_(nTot1,[0.0d0],0,Work(iTmp6),1)
*
        First=.True.
        Dff=.False.
        Do_DFT=.True.

        Call Put_iArray('nFro',nFro,nSym)
        Call Put_iArray('nAsh',nAsh,nSym)
        Call Put_iArray('nIsh',nIsh,nSym)

        call get_charge(iCharge)

        do_pdftPot=.false.
        if((DoGradPDFT.and.jroot == irlxroot).or.DoGradMSPD) then
          do_pdftPot=.true.

          CALL GETMEM('TE_POTG','ALLO','REAL',LTEOTPG,NFINT)
          CALL GETMEM('OE_POT','ALLO','REAL',LOEOTP,NTOT1)
          Call DCopy_(NTOT1,[0.0d0],0,Work(LOEOTP),1)
          Call DCopy_(NFINT,[0.0d0],0,work(LTEOTPG),1)

          Call Put_dArray('ONTOPO',work(LOEOTP),NTOT1)
          Call Put_dArray('ONTOPT',work(LTEOTPG),NFINT)

          Call DCopy_(NTOT1,[0.0d0],0,Work(LOEOTP),1)
          Call DCopy_(NFINT,[0.0d0],0,work(LTEOTPG),1)

          Call Put_dArray('FI_V',work(LOEOTP),NTOT1)
          Call Put_dArray('FA_V',work(LOEOTP),NTOT1)

          CALL GETMEM('OE_POT','FREE','REAL',LOEOTP,NTOT1)
          CALL GETMEM('TE_POTG','FREE','REAL',LTEOTPG,NFINT)
        end if !DoGradPdft

        Call DrvXV(Work(iTmp5),Work(iTmp6),Work(iTmp3),
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,
     &             KSDFT_TEMP,ExFac,iCharge,iSpin,
     &             Work(iD1I),Work(iD1ActAO),
     &             nTot1,DFTFOCK,Do_DFT)

        Call Daxpy_(nTot1,1.0d0,Work(iTmp5),1,Work(iTmp1),1)
        Call Daxpy_(nTot1,1.0d0,Work(iTmp6),1,Work(iFockI),1)

        Call GetMem('gtmp','Free','Real',iTmp6,nTot1)
        Call GetMem('htmp','Free','Real',iTmp5,nTot1)

!
***********************************************************
*     Compute energy contributions
***********************************************************
        iTmp2=0
        Call GetMem('DoneI','Allo','Real',iTmp2,nTot1)

        Call Fold(nSym,nBas,Work(iD1I),Work(iTmp2))

        Call GetMem('DoneA','Allo','Real',iTmpa,nTot1)
        Call Fold(nSym,nBas,Work(iD1ActAO),Work(iTmpa))

        Eone = dDot_(nTot1,Work(iTmp2),1,Work(iTmp1),1)
        Etwo = dDot_(nTot1,Work(iTmp2),1,Work(iFockI),1)

!**************Kinetic energy of inactive electrons********
        Ekin = dDot_(nTot1,Work(iTmp2),1,Work(iTmpk),1)

!*****Nuclear electron attraction for inactive electrons******
        Enuc = dDot_(nTot1,Work(iTmp2),1,Work(iTmpn),1)

c**************Kinetic energy of active electrons*********
        EactK = dDot_(nTot1,Work(iTmpk),1,Work(iTmpa),1)

        EactN = dDot_(nTot1,Work(iTmpn),1,Work(iTmpa),1)
        EMY  = PotNuc+Eone+0.5d0*Etwo

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,'(4X,A35,F18.8)')
     &     'Nuclear repulsion energy :',PotNuc
          Write(LF,'(4X,A35,F18.8)')
     &     'One-electron kinetic core energy:',Ekin
          Write(LF,'(4X,A35,F18.8)')
     &      'Nuc-elec attraction core energy:',Enuc
          Write(LF,'(4X,A35,F18.8)') 'One-electron core energy:',Eone
          Write(LF,'(4X,A35,F18.8)') 'Two-electron core energy:',Etwo
          Write(LF,'(4X,A35,F18.8)') 'Total core energy:',EMY
          Write(LF,'(4X,A35,F18.8)') 'Active Kinetic energy:',EactK
          Write(LF,'(4X,A35,F18.8)')
     &     'Active nuc-elec attraction energy:',EactN
        End If
***********************************************************
* Printing matrices
***********************************************************
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' FI matrix in CASDFT_Terms only 2-electron terms'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',Work(IFockI+iOff-1),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do
        End If

        Call DaXpY_(nTot1,One,Work(iTmp1),1,Work(iFockI),1)

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*)' Inactive Fock matrix in AO basis in CASDFT_terms'
          write(LF,*)'(it already contains OneHam and TwoEl contrib.)'
          Write(LF,*)' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',Work(iFockI+iOff-1),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do
        End If

************************************************************************
* update and transform the Fock matrices FI and FA ----> FMAT routine
************************************************************************
        if(iprlev.ge.debug) then
          write(lf,*) 'id1act before reading in'
           do i=1,nacpar
             write(lf,*) work(id1act-1+i)
           end do
          write(6,*) 'lpuvx before tractl'
          do i=1,nfint
            write(6,*) work(lpuvx-1+i)
          end do
          write(6,*) 'ltuvx after tractl'
          do i=1,nacpr2
            write(6,*) work(ltuvx-1+i)
          end do
          write(6,*) 'id1actao before tractl'
          do i=1,ntot2
            write(6,*) work(id1actao-1+i)
          end do
          write(6,*) 'ifocki before tractl'
          do i=1,ntot1
            write(6,*) work(ifocki-1+i)
          end do
        end if

        Call GetMem('FockI_save','ALLO','Real',ifocki_tmp,ntot1)
        call  dcopy_(ntot1,work(ifocki),1,work(ifocki_tmp),1)

        if (iprlev.ge.debug) then
          write(6,*) 'ifocki_tmp before tractl'
          do i=1,ntot1
            write(6,*) work(ifocki_tmp-1+i)
          end do

          write(6,*) 'ifocka before tractl'
          do i=1,ntot1
            write(6,*) work(ifocka-1+i)
          end do
        end if

        Call GetMem('ltuvx_tmp','ALLO','Real',ltuvx_tmp,nacpr2)
        Call GetMem('lpuvx_tmp','ALLO','Real',lpuvx_tmp,nfint)
        CALL DCOPY_(nacpr2,[Zero],0,WORK(ltuvx_tmp),1)
        CALL DCOPY_(nfint,[Zero],0,WORK(lpuvx_tmp),1)

        CALL TRACTL2(CMO,WORK(LPUVX_tmp),WORK(LTUVX_tmp),
     &                WORK(iD1I),WORK(ifocki_tmp),
     &                WORK(iD1ActAO),WORK(ifocka),
     &                IPR,lSquare,ExFac)

        Call GetMem('FockI_Save','Free','Real',ifocki_tmp,ntot1)
        Call GetMem('ltuvx_tmp','Free','Real',ltuvx_tmp,nacpr2)
        Call GetMem('lpuvx_tmp','Free','Real',lpuvx_tmp,nfint)

        if(iprlev.ge.debug) then
          write(6,*) 'lpuvx after tractl'
          do i=1,nfint
            write(6,*) work(lpuvx-1+i)
          end do
          write(6,*) 'ltuvx after tractl'
          do i=1,nacpr2
            write(6,*) work(ltuvx-1+i)
          end do
          write(6,*) 'id1actao after tractl'
          do i=1,ntot2
            write(6,*) work(id1actao-1+i)
          end do
          write(6,*) 'ifocki after tractl'
          do i=1,ntot1
            write(6,*) work(ifocki-1+i)
          end do
          write(lf,*) 'ifocka after tractl'
          do i=1,ntot1
            write(lf,*) work(ifocka-1+i)
          end do
        end  if

        ! Here we compute the inactive-active contribution (the so
        ! called interaction energy (saved into VIA).
        call mma_allocate(scr1, ntot1, "scr1")
        call Fold(nSym, nBas, Work(iD1ActAO), scr1)
        VIA = ddot_(ntot1, work(iFockI), 1, scr1, 1)
        call mma_deallocate(scr1)

        ECAS = EMY + VIA
        If ( iPrLev.ge.DEBUG ) then
          Write(LF,'(A,E20.10)') ' Total core energy:            ',EMY
          Write(LF,'(A,E20.10)') ' inactive-active interaction:  ',VIA
          Write(LF,'(A,E20.10)') ' CAS energy (core+interaction):',ECAS
        End If

! Compute the active-space DmatDmat (D_pqD_rs)
        IF(ISTORP(NSYM+1) > 0) THEN
          CALL GETMEM('ISTRP','ALLO','REAL',LP,ISTORP(NSYM+1))
          CALL DmatDmat_m(Work(iD1Act),WORK(LP))
        END IF

        if(iprlev >= debug) then
          write(lf,*) 'dmatdmat'
          do i=1,istorp(nsym+1)
            write(lf,*) Work(LP-1+i)
          end do
        end if
! Now we calculate the active-space Classical Coulomb Energy
        active_coulomb = 0.0D0
        call active_classical_coulomb(work(LP), work(LPUVX),
     &                                active_coulomb)
        ECAS = ECAS + active_coulomb
! Grab the E_ot energy which is computed somewhere above
        CASDFT_Funct = 0
        Call Get_dScalar('CASDFT energy',CASDFT_Funct)
        CASDFT_E = ECAS + CASDFT_Funct

! Of course, now we deal with a hybrid functional
        IF(Do_Hybrid) THEN
          E_NoHyb = CASDFT_E
          CASDFT_E=Ratio_WF*Ref_Ener(jRoot)+(1-Ratio_WF)*E_NoHyb
        END IF

        Call Print_MCPDFT_2(CASDFT_E,PotNuc,EMY,ECAS,CASDFT_Funct,
     &         jroot,Ref_Ener)

        Energies(jroot)=CASDFT_E
        IF(Do_Rotate) Then
*JB       replacing ref_ener with MC-PDFT energy for MS-PDFT use
          Ref_Ener(jroot)=CASDFT_E
        ELSE
          ener(jroot,1)=CASDFT_E
        END IF

!At this point, the energy calculation is done.  Now I need to build the
!fock matrix if this root corresponds to the relaxation root.

        if((DoGradPDFT .and. jroot == irlxroot) .or. dogradmspd) then
          ! This converts
          Call Fmat_m(CMO,Work(lPUVX),Work(iD1Act),
     &               Work(ifocki),Work(iFockA))
          ! This not just does
          !   a) calculate the generalized Fock matrix and stores it in
          !      lfock
          ! In theory, we don't need the lfock unless we are doing
          ! gradients! So we should try and remove that..
          ! This can, and should, be wrapped up with the fock_update
          !  procedure
          CALL GETMEM('FOCK','ALLO','REAL',LFOCK,NTOT4)
          CALL FOCK_m(WORK(LFOCK),Work(iFockI),Work(iFockA),
     &        Work(iD1Act),WORK(LP),WORK(LPUVX))
        end if

!***********************************************************************
*
*            BUILDING OF THE NEW FOCK MATRIX                           *
*
************************************************************************
        if(DoGradPDFT .and. jroot == irlxroot) then

          Write(LF,*) 'Loading potentials for analytic gradients...'
!MCLR requires two sets of things:
!1. The effective one-body Fock matrix and the effective two-body fock
!matrix.  These are used to generate the CI gradient inside of MCLR
!2. The effective generalized fock matrix.  This is used to calculate
!the orbital gradient inside of MCLR and is also used in the
!renormalization/effective lagrangian part of the final gradient
!evalutation.

!I think the plan should be to add on the missing pieces (to Fock_occ)
!which come from the one- and two-electron potentials.  These pieces are
!given by
! F_{xy} = \sum_{p} V_{py} D_{px} + \sum_{pqr} 2v_{pqry}d_{pqrx}.

!
!      write(6,*) 'NACPAR (input fock)',nacpar
!      write(6,*) 'ntot1 (# of V, fock_occ)',ntot1
!      write(6,*) 'nfint (# of v)',nfint
cPS         call xflush(6)

!I will read in the one- and two-electron potentials here
          Call GetMem('ONTOPT','ALLO','Real',ipTmpLTEOTP,nfint)
          Call GetMem('ONTOPO','ALLO','Real',ipTmpLOEOTP,ntot1)
          Call FZero(Work(iptmplteotp),Nfint)
          Call FZero(Work(iptmploeotp),ntot1)

          Call Get_dArray('ONTOPT',work(ipTmpLTEOTP),NFINT)
          Call Get_dArray('ONTOPO',work(ipTmpLOEOTP),NTOT1)
!
          If(IPRLEV.ge.DEBUG ) then
            write(6,*) 'One-electron potentials'
            do i=1,ntot1
              write(6,*) Work(iptmploeotp-1+i)
            end do
            write(6,*) 'Two-electron potentials'
            do i=1,nfint
              if (abs(work(lpuvx-1+i)).ge.1d-10)then
                write(6,*) Work(iptmplteotp-1+i),work(lpuvx-1+i)
              end if
            end do
          end if
!______________________________________________________
!Grab the active-active part of the FI+FA matrix (currently held in the
!FA matrix) and place it in an array of size NACPAR.  Add the oeotp to
!it.  Write to file.
          If ( IPRLEV.ge.DEBUG ) then
            write(6,*) "FA+FI to send to MCLR"
            do i=1,Ntot1
              write(6,*) Work(ifocka-1+i)
            end do
          end if

          Call GetMem('F_ONE','ALLO','Real',iFone,NTOT1)
          CALL DCOPY_(NTOT1,[0.0D0],0,WORK(iFone),1)

!I think I need to generate FI, which will contain both the one-electron
!potential contribution and the V_kkpu contribution.

          CALL GETMEM('FI_V','ALLO','REAL',ifiv,Ntot1)
          Call Get_dArray('FI_V',work(ifiv),NTOT1)
          Call daxpy_(ntot1,1.0d0,Work(ifiv),1,Work(iFocka),1)
          Call daxpy_(ntot1,1.0d0,Work(iptmploeotp),1,Work(iFocka),1)

          i_off1=0
          i_off2=0
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            do i=1,iBas
              do j=1,i
                Work(iFone+i_off1) = Work(ifone+i_off1) +
     &                               Work(ifocka+i_off2)
                i_off1 = i_off1 + 1
                i_off2 = i_off2 + 1
              end do
            end do
          end do
          If (IPRLEV >= DEBUG) then
            write(6,*) 'F1 to send'
            do i=1,NTOT1
              write(6,*) work(iFone-1+i)
            end do
          end if

      !Add the V_kktu contribution to Fone_tu?
!STILL MUST DO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This should be addressed in the upd_FI routine.

          !Write to file.
          LUTMP=87
          LUTMP=IsFreeUnit(LUTMP)
          Call Molcas_Open(LUTMP,'TmpFock')
          do i=1,ntot1
            write(LUTMP,*) Work(iFone+i-1)
          end do

          Call GetMem('ttTUVX','Allo','Real',ittTUVX,NACPR2)
          CALL DCOPY_(nacpr2,[0.0D0],0,WORK(ittTUVX),1)
          Call Get_TUVX(Work(ipTmpLTEOTP),Work(ittTUVX))

          !Unpack TUVX to size
          do i=1,nacpr2
            write(LUTMP,*) Work(ittTUVX+i-1)
          end do
          Close(LUTMP)
          Call GetMem('F_ONE','Free','Real',iFone,NTOT1)
          Call GetMem('ttTUVX','Free','Real',ittTUVX,NACPR2)

!____________________________________________________________
!This next part is to generate the MC-PDFT generalized fock operator.

!Zero out the matrices.  We will be adding the potential-containing
!terms as a correction to the Focc component already on the runfile.
          CALL DCOPY_(ntot1,[0.0D0],0,WORK(iFocka),1)
          CALL DCOPY_(ntot1,[0.0D0],0,WORK(iFocki),1)

!The corrections (from the potentials) to FI and FA are built in the NQ
!part of the code, for efficiency's sake.  It still needs to be
!debugged.
          CALL GETMEM('FA_V','ALLO','REAL',ifav,Ntot1)
          Call Get_dArray('FA_V',work(ifav),NTOT1)

          If(IPRLEV.ge.DEBUG) then
            write(lf,*) "extra terms to update FI"
            do i=1,ntot1
              write(lf,*) Work(ifiv-1+i)
            end do
            write(lf,*) "extra terms to update FA"
            do i=1,ntot1
              write(lf,*) Work(ifav-1+i)
            end do
          end if

          If(IPRLEV >= DEBUG ) then
            CALL GETMEM('FA_t','ALLO','REAL',ifat,Ntot1)
            Call dcopy_(ntot1,[0.0d0],0,work(ifat),1)
            Call DaXpY_(NTOT1,1.0D0,Work(ipTmpLOEOTP),1,Work(ifat),1)
            Call daxpy_(NTOT1,1.0D0,Work(ifiv),1,Work(ifat),1)
            Call daxpy_(NTOT1,1.0D0,Work(ifav),1,Work(ifat),1)
            write(lf,*) "Total F additions:"
            Call TriPrt(' ','(5G18.10)',Work(ifat),norb(1))
            CALL GETMEM('FA_t','free','REAL',ifat,Ntot1)
          end if

      !Add one e potential, too.
          Call DaXpY_(NTOT1,1.0D0,Work(ipTmpLOEOTP),1,Work(ifocki),1)

      !Add two e potentials
          Call daxpy_(NTOT1,1.0D0,Work(ifiv),1,Work(ifocki),1)
          Call daxpy_(NTOT1,1.0D0,Work(ifav),1,Work(ifocka),1)
          If(IPRLEV >= DEBUG) then
            write(6,*) "new FI"
            Call TriPrt(' ','(5G18.10)',Work(ifocki),norb(1))
            write(6,*) "new FA"
            Call TriPrt(' ','(5G18.10)',Work(ifocka),norb(1))
          end if

          CALL GETMEM('FI_V','Free','REAL',ifiv,Ntot1)
          CALL GETMEM('FA_V','Free','REAL',ifav,Ntot1)

          IF(ISTORP(NSYM+1) > 0) THEN
            CALL DCOPY_(ISTORP(NSYM+1),[0.0D0],0,WORK(LP),1)
            CALL PMAT_RASSCF_M(Work(iP2d),WORK(LP))
          END IF

!Must add to existing FOCK operator (occ/act). FOCK is not empty.
! delete NSXS (LBM) or BM
          CALL GETMEM('SXLQ','ALLO','REAL',LQ,NQ) ! q-matrix(1symmblock)
          CALL FOCK_update(WORK(LFOCK),Work(iFockI),
     &        Work(iFockA),Work(iD1Act),WORK(LP),
     &        WORK(LQ),WORK(ipTmpLTEOTP), cmo)

          Call Put_dArray('FockOcc',Work(ipFocc),ntot1)
          If (IPRLEV >= DEBUG) then
            write(lf,*) 'FOCC_OCC'
            call wrtmat(Work(ipFocc),1,ntot1,1,ntot1)
            write(lf,*) 'DONE WITH NEW FOCK OPERATOR'
          end if

          CALL GETMEM('SXLQ','Free','REAL',LQ,NQ) ! q-matrix(1symmblock)
          Call GetMem('ONTOPO','FREE','Real',ipTmpLOEOTP,ntot1)
          Call GetMem('ONTOPT','FREE','Real',ipTmpLTEOTP,nfint)

!Put information needed for geometry optimizations.
          iSA = 1 !need to do MCLR for gradient runs. (1 to run, 2 to
*skip)
       !MUST MODIFY THIS.  I need to check that the calculation is not
       !SA, and if it is, set iSA to -1.
          Call Put_iScalar('SA ready',iSA)
          Call Put_cArray('MCLR Root','****************',16)
          Call Put_iScalar('Relax CASSCF root',irlxroot)
        end if !DoGradPDFT

        if (dogradmspd) then
*        doing exactly the same thing as done in the previous chunck
*        starting from 'BUILDING OF THE NEW FOCK MATRIX'
*        Hopefully this code will be neater.
          call savefock_pdft(CMO,IFockI,IFockA,iD1Act,LFock,
     &                        LP,NQ,LPUVX,ip2d,jroot)
        end if

        if((DoGradPDFT .and. jroot == irlxroot) .or. dogradmspd) then
          CALL GETMEM('FOCK','Free','REAL',LFOCK,NTOT4)
        end if

        Call GetMem('DoneI','Free','Real',iTmp2,nTot1)
        Call GetMem('DoneA','Free','Real',iTmpa,nTot1)

        IF(ISTORP(NSYM+1) > 0) THEN
          CALL GETMEM('ISTRP','FREE','REAL',LP,ISTORP(NSYM+1))
        END IF
      end do !loop over roots

      Call Put_iScalar('Number of roots',nroots)
      Call Put_dArray('Last energies',Energies,nroots)
      Call Put_cArray('Relax Method','MCPDFT  ',8)
      Call Put_dScalar('Last energy',Energies(iRlxRoot))

      ! Put the densities into the runfile for MCLR
      if(DoGradPDFT .or. DoGradMSPD) then
        dmDisk = IADR19(3)
        do jroot=1,irlxroot-1
          Call DDaFile(JOBOLD,0,Work(iD1Act),NACPAR,dmDisk)
          Call DDaFile(JOBOLD,0,Work(iD1Spin),NACPAR,dmDisk)
          Call DDaFile(JOBOLD,0,Work(iP2d),NACPR2,dmDisk)
          Call DDaFile(JOBOLD,0,Work(iP2d),NACPR2,dmDisk)
        end do

        Call DDaFile(JOBOLD,2,Work(iD1Act),NACPAR,dmDisk)
        Call DDaFile(JOBOLD,2,Work(iD1Spin),NACPAR,dmDisk)
        Call DDaFile(JOBOLD,2,Work(iP2d),NACPR2,dmDisk)
        Call DDaFile(JOBOLD,0,Work(iP2d),NACPR2,dmDisk)

        Call Put_dArray('D1mo',Work(iD1Act),NACPAR)
        Call Put_dArray('P2mo',Work(iP2d),NACPR2)

        If(NASH(1) /= NAC) then
          Call DBLOCK_m(Work(iD1Act))
          CALL DBLOCK_m(Work(iD1Spin))
        end if

        Call cas_mo_to_ao(CMO,Work(iD1Act),Work(iD1ActAO))
        Call Fold(nSym,nBas,Work(iD1I),Work(iTmp3))
        Call Fold(nSym,nBas,Work(iD1ActAO),Work(iTmp4))
        Call Daxpy_(nTot1,1.0D0,Work(iTmp4),1,Work(iTmp3),1)
        Call Put_dArray('D1ao',Work(iTmp3),nTot1)

!Get the spin density matrix for open shell cases
        Call cas_mo_to_ao(CMO,Work(iD1Spin),
     &                     Work(iD1SpinAO))
        Call GetMem('DtmpS','Allo','Real',iTmp7,nTot1)
        Call Fold(nSym,nBas,Work(iD1SpinAO),Work(iTmp7))
        Call Put_dArray('D1Sao',Work(iTmp7),nTot1)
        Call GetMem('DtmpS','Free','Real',iTmp7,nTot1)
      end if

      if(doGSOR) then
        LUGS=25
        LUGS=IsFreeUnit(LUGS)
        IAD19=0
        Call DaName(LUGS,'JOBGS')
        Call IDaFile(LUGS,2,IADR19,15,IAD19)
        Call DDAFile(LUGS,1,Energies,lroots,IADR19(6))
        Call DaClos(LUGS)
      end if

      if (nFint > 0) then
        Call GetMem('PUVX','FREE','Real',lPUVX,nFint)
      end if

!Free up all the memory we can here, eh?
      Call GetMem('DtmpA','Free','Real',iTmp4,nTot1)
      Call GetMem('DtmpI','Free','Real',iTmp3,nTot1)

      Call GetMem('D1Active','free','Real',iD1Act,NACPAR)
      Call GetMem('D1ActiveAO','free','Real',iD1ActAO,NTOT2)
      Call GetMem('D1Spin','free','Real',iD1Spin,NACPAR)
      Call GetMem('D1SpinAO','free','Real',iD1SpinAO,NTOT2)
      Call GetMem('Fcore','Free','Real',iTmp1,nTot1)
      Call GetMem('FockI','FREE','Real',ifocki,ntot1)
      Call GetMem('FockA','FREE','Real',ifocka,ntot1)
      Call GetMem('P2','Free','Real',iP2d,NACPR2)
      Call GetMem('D1Inact','Free','Real',iD1i,NTOT2)
      Call GetMem('Kincore','free','Real',iTmpk,nTot1)
      Call GetMem('NucElcore','free','Real',iTmpn,nTot1)
      Return
      END

      Subroutine P2_contraction(D1MO,P2MO)
      use definitions, only: wp
      implicit none

      real(kind=wp), dimension(*), intent(in) :: d1mo
      real(kind=wp), dimension(*), intent(out) :: p2mo

#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      integer :: i, j, k, l, ij, kl, ijkl, lmax
      real(kind=wp) :: fact

      ijkl=0
      do i=1,nac
        do j=1,i
          ij = iTrii(i,j)
          do k=1,i
            if(i == k) then
              lmax = j
            else
              lmax = k
            end if
            do l=1,lmax
              kl = iTrii(k,l)
              ijkl = ijkl + 1
              fact=1.0d0
              if(k == l) fact=0.5d0
              p2MO(ijkl) = fact*D1MO(ij)*D1MO(kl)
            end do
          end do
        end do
      end do
      contains
        integer function iTrii(i,j)
          integer, intent(in) :: i, j
          itrii = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
        end function
      end subroutine


      Subroutine Fold_pdft(nSym,nBas,A,B)

      Implicit Real*8 (A-H,O-Z)

      Dimension nBas(*) , A(*) , B(*)

      iOff1 = 0
      iOff2 = 0
      Do iSym = 1, nSym
        mBas = nBas(iSym)
        Do iBas= 1, mBas
          Do jBas = 1 , iBas-1
            B(iOff2+jBas) =   A(iOff1+jBas)
          End Do
          B(iOff2+iBas) =  A(iOff1+iBas)
          iOff1 = iOff1 + mBas
          iOff2 = iOff2 + iBas
        End Do
      End Do

      Return
      end
************ columbus interface ****************************************
