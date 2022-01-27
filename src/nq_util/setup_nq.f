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
* Copyright (C) 1999, Roland Lindh                                     *
************************************************************************
      Subroutine Setup_NQ(Maps2p,nShell,nSym,nNQ,Do_Grad,On_Top,nD,
     &                    Pck_Old,PMode_old,R_Min,nR_Min)
************************************************************************
*                                                                      *
* Object: to set up information for calculation of integrals via a     *
*         numerical quadrature.                                        *
* Warning: The exponents of each shell are reordered diffuse to compact*
*                                                                      *
*     Author: Roland Lindh,                                            *
*             Dept of Chemical Physics,                                *
*             University of Lund, Sweden                               *
*             August 1999                                              *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iOper
      use nq_Grid, only: nGridMax, Coor, Pax, Fact, Tmp, nR_Eff, SOs
      use nq_Grid, only: Angular, Mem
      use nq_structure, only: NQ_Data
      use Grid_On_Disk
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "status.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
#include "print.fh"
      Real*8 XYZ(3), C(3)
      Logical EQ, Do_Grad, On_Top, PMode_Old
      Real*8 Alpha(2),rm(2), R_Min(0:nR_Min)
      Integer Maps2p(nShell,0:nSym-1)
      Integer iDCRR(0:7)
      Dimension Dummy(1)
      Real*8, Allocatable:: TempC(:,:), ZA(:), Crd(:,:), dOdx(:,:,:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement Functions
*
      nElem(i)=(i+1)*(i+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      Call ICopy(nShell*nSym,[-99999999],0,Maps2p,1)
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
c     Write(6,*) '********** Setup_NQ ***********'
      ntotgp=0
*                                                                      *
************************************************************************
*                                                                      *
*-----Check if NQ environment has been activated
*
      If (NQ_Status.ne.Active.and.NQ_Status.ne.Inactive) Then
         Call WarningMessage(2,'Setup_NQ: NQ_Status not initialized')
         Call Quit_OnUserError()
      End If
      If (NQ_Status.eq.Active) Return
      NQ_Status=Active
*                                                                      *
************************************************************************
*                                                                      *
*---- Get the coordinates to the centers of all Voronoi polyhedra
*
*     Note that this will be all centers with valence basis sets on
*     them. Hence this will also include any pseudo centers!
*
      Call mma_allocate(TempC,3,nShell*nSym,Label='TempC')
      nAtoms = 0
      If (nShell.gt.nskal_iSD) Then
         Write (6,*) 'nShell.gt.nSkal_iSD'
         Write (6,*) 'nShell=',nShell
         Write (6,*) 'nSkal_iSD=',nSkal_iSD
         Call AbEnd()
      End If
      Do iShell = 1, nShell
         iCnttp=iSD(13,iShell)
         iCnt  =iSD(14,iShell)
         XYZ(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
         Call Process_Coor(XYZ,TempC,nAtoms,nSym,iOper)
      End Do
      Call mma_allocate(Coor,3,nAtoms,Label='Coor')
      call dcopy_(3*nAtoms,TempC,1,Coor,1)
      Call mma_deallocate(TempC)
*                                                                      *
************************************************************************
*                                                                      *
*---- Get the symmetry unique coordinates
*
      nNQ=nAtoms
      Allocate(NQ_data(1:nNQ))
      Do iNQ = 1, nNQ
         Call mma_allocate(NQ_data(iNQ)%Coor,3,
     &              Label='NQ_data(iNQ)%Coor')
         call dcopy_(3,Coor(1:3,iNQ),1,NQ_data(iNQ)%Coor,1)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*-----Pick up the requested accuracy.
*
*     Dr=-Log10(Thr)
*                                                                      *
************************************************************************
*                                                                      *
*-----Loop over each unique center, and find the highest and lowest    *
*     Gaussians exponents associated with this center. The later       *
*     information will be used to design the radial grid associated    *
*     with this center.                                                *
*                                                                      *
      iAngMax=0
      NbrMxBas=0
      Do iShell=1,nShell
         iAng  =iSD(1,iShell)
         nCntrc=iSD(3,iShell) !Get the # of contracted functions
                              !for iShell
         mExp  =iSD(5,iShell) ! Get the number of exponents of ishell
         NbrMxBas=Max(NbrMxbas,nCntrc)
         iAngMax=Max(iAngMax,iAng)
      End Do
*
*-----Loop over the shells
*
      nMaxExp=0
      nAOMax=0
      Do iShell = 1, nShell
*
*------- Get the Atom number
         iANr=dbsc(iSD(13,iShell))%AtmNr
*
         iShll=iSD(0,iShell)   ! Get the angular momentum of ishell
         iAng=iSD(1,iShell)   ! Get the angular momentum of ishell
         iCmp=iSD(2,iShell)   ! Get the # of angular components
         nCntrc=iSD(3,iShell) !Get the # of contracted functions
                              !for iShell
         mExp=iSD(5,iShell)   ! Get the number of exponents of ishell
         nMaxExp=Max(nMaxExp,mExp)
         nAOMax=Max(nAOMax,iCmp*nCntrc)
*
************************************************************************
*                                                                      *
*     Order the exponents diffuse to compact for the active shell      *
*                                                                      *
************************************************************************
         Call OrdExpD2C(mExp,Shells(iShll)%Exp,nCntrc,
     &                       Shells(iShll)%pCff)
*
*-----Get the extreme exponents for the active shell.
         A_low =Shells(iShll)%Exp(1)
         A_high=Shells(iShll)%Exp(mExp)
*
         iCnttp=iSD(13,iShell)
         iCnt  =iSD(14,iShell)
         C(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
         Do iIrrep = 0, nIrrep-1
            Call OA(iOper(iIrrep),C,XYZ)
         Do iNQ=1,nNQ
*
            If (EQ(NQ_data(iNQ)%Coor,XYZ)) Then
*
               NQ_Data(iNQ)%Atom_Nr=iANR
*
*------------- Assign the BS radius to the center
               NQ_Data(iNQ)%R_RS=Bragg_Slater(iANr)
*
*------------- What is the maximum angular momentum for the active center ?
               NQ_Data(iNQ)%l_Max=Max(NQ_Data(iNQ)%l_max,iAng)
*
*------------- Get the extreme exponents for the atom
               NQ_Data(iNQ)%A_high=Max(NQ_Data(iNQ)%A_high,A_High)
               NQ_Data(iNQ)%A_low =Min(NQ_Data(iNQ)%A_low ,A_low)
*
               Maps2p(iShell,iIrrep)=iNQ
               Go To 100
            End If
         End Do                 ! iNQ
         Call WarningMessage(2,
     &              'Didn''t find a center associated with the shell!')
         Call Abend()
 100     Continue
         End Do                 ! iIrrep
*
      End Do                    !iShell
*
************************************************************************
*                                                                      *
*               END OF THE LOOP OVER THE SHELLS.                       *
*     Now we have the number of unique center and their associated     *
*     exponents and maximum angular momentum. And for each shell the   *
*     exponents are ordered diffuse to compact.                        *
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over all the atoms to create the radial quadrature          *
*                                                                      *
************************************************************************
*
*-----Allocate memory to store the number of effective radial points for
*     each center and the radius of this center.
      Call mma_Allocate(nR_Eff,nNQ,Label='nR_Eff')
*
      iNQ_MBC=0
      iReset=0
      Threshold_tmp=Zero
      nR_tmp=0
      Do iNQ=1,nNQ
*--------Get the extreme exponents for the atom
         Alpha(1)=NQ_Data(iNQ)%A_low
         Alpha(2)=NQ_Data(iNQ)%A_high
*
*--------Get the coordinates of the atom
         call dcopy_(3,NQ_Data(iNQ)%Coor,1,XYZ,1)
*
*        For a special center we can increase the accuracy.
*
         If (MBC.ne.' ') Then
            Do iS = 1, nShell
               iCnttp=iSD(13,iS)
               iCnt  =iSD(14,iS)
               C(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
               If ( EQ(NQ_Data(iNQ)%Coor,C) ) Then
                  mdci=iSD(10,iS)
                  If (dc(mdci)%LblCnt.eq.MBC) Then
                     nR_tmp=nR
                     nR=INT(DBLE(nR)*2.0D0)
                     Threshold_tmp=Threshold
                     Threshold=Threshold*1.0D-6
*
                     iReset=1
                     iNQ_MBC=iNQ
                     Go To 1771
                  End If
               End If
            End Do
         End If
 1771    Continue
*
*        Max angular momentum for the atom -> rm(1)
*        Max Relative Error -> rm(2)
         rm(1)=DBLE(NQ_Data(iNQ)%l_Max)
         rm(2)=Threshold
*
         Call GenVoronoi(XYZ,nR_Eff,nNQ,Alpha,rm,iNQ)
*
         If (iReset.eq.1) Then
            nR=nR_tmp
            Threshold=Threshold_tmp
            iReset=0
         End If
*
      End Do  ! iNQ
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the principal axis system and optionally to
*     compute derivatives of the principal axis. Needed in order to
*     compute the gradient of the rotationally invariant DFT energy.
*
      Call mma_allocate(Pax,3,3,Label='Pax')
      Call mma_allocate(dOdx,3,3,nNQ,3,Label='dOdx')
      dOdx(:,:,:,:)=Zero
      Call mma_Allocate(ZA,nNQ,Label='ZA')
      Call mma_Allocate(Crd,3,nNQ)
*
*     Collect coordinates and charges of the nuclei
*
      Do iNQ = 1, nNQ
         ZA(iNQ)=DBLE(NQ_Data(iNQ)%Atom_Nr)
         call dcopy_(3,NQ_data(iNQ)%Coor,1,Crd(:,iNQ),1)
      End Do
*
      Call RotGrd(Crd,ZA,Pax,dOdx,Dummy,nNQ,Do_Grad,.False.)
*
*     Distribute derivative of the principle axis system
*
      If (Do_Grad) Then
      Do iNQ = 1, nNQ
         Call mma_allocate(NQ_Data(iNQ)%dOdx,3,3,3,Label='dOdx')
         Do iCar = 1, 3
           call dcopy_(9,dOdx(:,:,iNQ,iCar),1,
     &                   NQ_Data(iNQ)%dOdx(:,:,iCar),1)
         End Do
      End Do
      End If
*
      Call mma_deallocate(dOdX)
      Call mma_deallocate(Crd)
      Call mma_deallocate(ZA)
*                                                                      *
************************************************************************
*                                                                      *
      If (Rotational_Invariance.eq.Off) Then
         Call FZero(Pax,9)
         call dcopy_(3,[One],0,Pax,4)
         Do iNQ = 1, nNQ
            If (.Not.Allocated(NQ_Data(iNQ)%dOdx))
     &         Call mma_allocate(NQ_Data(iNQ)%dOdx,3,3,3,Label='dOdx')
               NQ_Data(iNQ)%dOdx(:,:,:)=Zero
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the angular grid
*
      Call Angular_grid()
*
      Crowding_tmp=Zero
      Do iNQ = 1, nNQ
*
*        Assign the angular grid to be used with each radial grid point
*
         Call mma_allocate(NQ_Data(iNQ)%Angular,nR_Eff(iNQ),
     &                     Label='Angular')
         NQ_Data(iNQ)%Angular(:)=nAngularGrids
*
*        Prune the angular grid
*
         If (Angular_Prunning.eq.On) Then
*
*
*---------- Find the R_min values of each angular shell
*
            lAng=NQ_Data(iNQ)%l_max
            Do iAng = 0, lAng
               R_Min(iAng)=Zero
               ValExp=-One
               iSet=0
               Do iShell=1,nShell
                  iShll =iSD(0,iShell)
                  iAng_ =iSD(1,iShell)
                  NrExp =iSD(5,iShell)
*                 Write (6,*) 'iAng_,iAng=',iAng_,iAng
                  If (iAng_.eq.iAng.and.NrExp.ge.1) Then
                     Do iSym = 0, nSym-1
                        iNQ_=Maps2p(iShell,iSym)
*                       Write (6,*) 'iNQ_,iNQ=',iNQ_,iNQ
                        If (iNQ_.eq.iNQ) Then
                           ValExp=Shells(iShll)%Exp(NrExp)
                           iSet=1
                        End If
                     End Do
                  End If
               End Do
               If (ValExp.lt.Zero.and.iSet.eq.1) Then
                  Call WarningMessage(2,'ValExp.lt.Zero')
                  Call Abend()
               End If
               If (iSet.eq.1) Then
                  R_Min(iAng)=Eval_RMin(ValExp,iAng,Threshold)
                  If (iAng.eq.0) R_Min(iAng)=Zero
               End If
            End Do
*
            R_BS = NQ_Data(iNQ)%R_RS
*
            If (iNQ.eq.iNQ_MBC) Then
               Crowding_tmp=Crowding
               Crowding=One + (Crowding-One)*0.25D0
               iReset=1
            End If
*
            Call Angular_Prune(NQ_Data(iNQ)%R_Quad,nR_Eff(iNQ),
     &                         NQ_Data(iNQ)%Angular,Crowding,
     &                         Fade,R_BS,L_Quad,R_Min,lAng,
     &                         nAngularGrids)
*
            If (iReset.eq.1) Then
               Crowding=Crowding_tmp
               iReset=0
            End If
*
         End if
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,'(A)') ' =================================='
      Write (6,'(A)') ' =        Grid information        ='
      Write (6,'(A)') ' =================================='
      Write (6,'(A)') ' Legend             '
      Write (6,'(A)') ' ----------------------------------'
      Write (6,'(A)') ' ANr: element number'
      Write (6,'(A)') ' nR : number of radial grid points'
      Write (6,'(A)') ' iNQ: grid index'
      Write (6,'(A)') ' ----------------------------------'
      Write (6,*)
      Write (6,'(A)') ' iNQ ANr  nR'
      Do iNQ=1,nNQ
         iANr=NQ_Data(iNQ)%Atom_Nr
         kR=nR_Eff(iNQ)
         Write (6,'(3I4)') iNQ, iANr, kR
      End Do
      Write (6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Determine the spatial extension of the molecular system
*
!     Box_Size=Four           ! Angstrom
      Box_Size=Two            ! Angstrom
!     Box_Size=1.0d0/Two      ! Angstrom
      Block_size=Box_Size
      x_min=1.0D99
      y_min=1.0D99
      z_min=1.0D99
      x_max=-1.0D99
      y_max=-1.0D99
      z_max=-1.0D99
      Do iAt = 1, nAtoms
         x_min=Min(x_min,Coor(1,iAt))
         y_min=Min(y_min,Coor(2,iAt))
         z_min=Min(z_min,Coor(3,iAt))
         x_max=Max(x_max,Coor(1,iAt))
         y_max=Max(y_max,Coor(2,iAt))
         z_max=Max(z_max,Coor(3,iAt))
      End Do
*
*---- Add half an box size around the whole molecule
*
      x_min=x_min-Box_Size/Two
      y_min=y_min-Box_Size/Two
      z_min=z_min-Box_Size/Two
      x_max=x_max+Box_Size/Two
      y_max=y_max+Box_Size/Two
      z_max=z_max+Box_Size/Two
*
*---- At least one finite box. Adjust to an even number of boxes.
*
      nx=Int((x_max-x_min+Box_Size)/Box_Size)
      nx=2*((nx+1)/2)
      ny=Int((y_max-y_min+Box_Size)/Box_Size)
      ny=2*((ny+1)/2)
      nz=Int((z_max-z_min+Box_Size)/Box_Size)
      nz=2*((nz+1)/2)
*
*---- Adjust extremal values to fit exactly with the
*     box size.
*
      dx=(DBLE(nx)*Box_Size-(x_max-x_min))/Two
      dy=(DBLE(ny)*Box_Size-(y_max-y_min))/Two
      dz=(DBLE(nz)*Box_Size-(z_max-z_min))/Two
*
      x_min=x_min-dx
      y_min=y_min-dy
      z_min=z_min-dz
      x_max=x_max+dx
      y_max=y_max+dy
      z_max=z_max+dz
*
*---- Add the infinite edge boxes
*
      nx=nx+2
      ny=ny+2
      nz=nz+2
#ifdef _DEBUGPRINT_
      Write (6,*) 'x_min=',x_min,dx
      Write (6,*) 'y_min=',y_min,dy
      Write (6,*) 'z_min=',z_min,dz
      Write (6,*) 'x_max=',x_max
      Write (6,*) 'y_max=',y_max
      Write (6,*) 'z_max=',z_max
      Write (6,*) 'nx,ny,nz=',nx,ny,nz
      Write (6,*) 'Total number of blocks=',nx*ny*nz
#endif
      number_of_subblocks=nx*ny*nz
*                                                                      *
************************************************************************
*                                                                      *
*     nFOrd: the order of the functional. nFOrd-1 is the number of times
*            the basis functions has to be differentiated to compute the
*            energy contribution.
*     mTmp: seems to be redundant?
*     mRad: number of different radial functions associated with a
*           basis function. This number depends on the type of
*           functional and the number of times the basis function has
*           to be differentiated in order to produce the values of the
*           parameters which the functional depends on (rho, grad rho,
*           and nabla rho).
*     nScr: Used to assemble integrals, (rho, grad rho, nabla rho)
*           (1,4, or 5). This is not needed in gradient calculations!
*     mAO: number of elements a basis function generates upon
*          differentiation (1,4,10,20, etc.)
*
      If (Functional_type.eq.LDA_type) Then
         nFOrd=1
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         if(do_grad) mAO=4!AMS - GRADIENTS?
         If (.Not.Do_Grad) Then
            mTmp=1
            mRad=nFOrd
            nScr=nD*nAOMax
         Else
            mTmp=7
            mRad=nFOrd+1
            nScr=0
         End If
*
      Else If (Functional_type.eq.GGA_type) Then
         nFOrd=2
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         if(do_grad) mAO=10
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*4*nAOMax
         Else
            mTmp=28
            mRad=nFOrd+1
            nScr=0
         End If
      Else If (Functional_type.eq.CASDFT_type) Then
*        I need to discuss this with Sergey!
*        nFOrd=3 !?
*        mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         nFOrd=2
         mAO=10
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*4*nAOMax
         Else
            mTmp=28
            mRad=nFOrd+1
            nScr=0
         End If
*
      Else If (Functional_type.eq.meta_GGA_type1) Then
         nFOrd=2
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*4*nAOMax
         Else
            mTmp=28
            mRad=NFOrd+1
            nScr=0
         End If
*
      Else If (Functional_type.eq.meta_GGA_type2) Then
         nFOrd=3
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*5*nAOMax
         Else
            mTmp=28
            mRad=NFOrd+1
            nScr=0
         End If
*
*     Else If (Functional_type.eq.Other_type) Then
      Else
         mTmp=0 ! Dummy initialize
         mRad=0 ! Dummy initialize
         mAO=0  ! Dummy initialize
         nScr=0 ! Dummy initialize
         Call WarningMessage(2,'Functional_type.eq.Other_type')
         Call Abend()
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     Allocate scratch for processing AO's on the grid
!
      nTmp=nGridMax*Max(mTmp,nScr)
      Call mma_allocate(Tmp,nTmp,Label='Tmp')
*
      nMem=0
      nSO =0
      lSO =0
      lAngular=0
      Do ish = 1, nShell
         iAng  = iSD( 1,iSh)
         iCmp  = iSD( 2,iSh)
         iBas  = iSD( 3,iSh)
         iPrim = iSD( 5,iSh)
*
         nxyz    = nGridMax*3*(iAng+mRad)
         nDrv     = mRad - 1
         nForm    = 0
         Do iDrv  = 0, nDrv
            nForm = nForm + nElem(iDrv)
         End Do
         nTerm    = 2**nDrv
         nAngular = 5*nForm*nTerm
         nRad    = iPrim*nGridMax*mRad
         nRadial = iBas*nGridMax*mRad
         If (On_Top) Then
            mdci  = iSD(10,iSh)
            kAO=iCmp*iBas*nGridMax
            nSO=kAO*nSym/dc(mdci)%nStab*mAO
         End If
         nMem=Max(nMem,nxyz+nRad+nRadial)
         lSO=Max(lSO,nSO)
         lAngular=Max(lAngular,nAngular)
      End Do
*
      Call mma_allocate(SOs,2*lSO,Label='SOs')
      Call mma_allocate(Angular,lAngular,Label='Angular')
      Call mma_allocate(Mem,nMem,Label='Mem')
*                                                                      *
************************************************************************
*                                                                      *
*     Access the file with Grid points and weights.
*
*---- Open the file.
      Lu_Grid=88
      Call DaName_MF_WA(Lu_Grid,'NQGRID')
*
      If (iGrid_Set.eq.Not_Specified) iGrid_Set=Final
*
*---- Read the status flag.
      iDisk_Grid=0
      Call iDaFile(Lu_Grid,2,G_S,5,iDisk_Grid)
*
      Grid_Status=G_S(iGrid_Set)
      If (Old_Functional_Type.ne.Functional_Type) Then
         G_S(Final)=Regenerate
         G_S(Intermediate)=Regenerate
         Grid_Status=Regenerate
      End If
      iDisk_Grid=iDisk_Set(iGrid_Set)
*
*---- Allocate memory for the master TOC.
      Call mma_Allocate(GridInfo,2,number_of_subblocks,
     &                  Label='GridInfo')
*
*---- Retrieve the TOC or regenerate it.
*
*     The table contains two data items per subblock.
*     1) disk address and 2) number of batches.
*
      If (Grid_Status.eq.Regenerate) Then
C        Write (6,*) 'Grid_Status.eq.Regenerate'
         Grid_Status=Regenerate
         GridInfo(:,:)=0
         Call iDaFile(Lu_Grid,1,GridInfo,
     &                2*number_of_subblocks,iDisk_Grid)
         Old_Functional_Type=Functional_Type
      Else If (Grid_Status.eq.Use_Old) Then
C        Write (6,*) 'Grid_Status.eq.Use_Old'
         Call iDaFile(Lu_Grid,2,GridInfo,
     &                2*number_of_subblocks,iDisk_Grid)
      Else
         Call WarningMessage(2,'Illegal Grid Status!')
         Call Abend()
      End If
*
      Call ParmPkR8(Pck_Old,PMode_old)
      Call IniPkR8(T_X,.True.)
*                                                                      *
************************************************************************
*                                                                      *
*     Setup some symmetry stuff outside the loop
*
      ndc = 0
      Do iSh = 1, nShell
         ndc = Max(ndc,iSD(10,iSh))
      End Do
      Call mma_allocate(Fact,ndc,ndc,Label='Fact')
      Do mdci = 1, ndc
         nDegi=nIrrep/dc(mdci)%nStab
         Do mdcj = 1, ndc
            nDegj=nIrrep/dc(mdcj)%nStab
*
            Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,
     &                     dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
*
            iuv = dc(mdci)%nStab*dc(mdcj)%nStab
            If (MolWgh.eq.1) Then
               Fct = DBLE(nIrrep) / DBLE(LmbdR)
            Else If (MolWgh.eq.0) Then
               Fct = DBLE(iuv) / DBLE(nIrrep * LmbdR)
            Else
               Fct = Sqrt(DBLE(iuv))/ DBLE(LmbdR)
            End If
            Fct=Fct*DBLE(nDCRR)/DBLE(nDegi*nDegj)
*
*---------- Save: Fact
*
            Fact(mdci,mdcj) = Fct
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
