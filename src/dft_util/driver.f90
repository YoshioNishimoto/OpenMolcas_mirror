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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
Subroutine Driver(KSDFA,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nD,DFTFOCK)

use libxc_parameters
use OFembed, only: KEOnly, dFMD, Do_Core
use libxc,   only: Only_exc
use nq_Grid, only: l_casdft
use nq_pdft, only: lft
Implicit None
#include "real.fh"
#include "nq_info.fh"
Character*(*) KSDFA
Logical Do_Grad
Integer :: nGrad, nh1, nD
Real*8 :: Func, Grad(nGrad)
Logical Do_MO, Do_TwoEl
Real*8 :: D_DS(nh1,nD), F_DFT(nh1,nD)
Real*8 :: A=Zero, B=Zero, C=Zero
Character*4 DFTFOCK
logical :: LDTF=.False., NDSD=.False.
character(LEN=12) :: FLabel=''

abstract interface
   Subroutine DFT_FUNCTIONAL(mGrid,nD)
      Integer mGrid, nD
   end subroutine
end interface

!***********************************************************************
!     Define external functions not defined in LibXC. These are either
!     accessed through the procedure pointer sub or External_sub.

procedure(DFT_FUNCTIONAL) :: Overlap, NucAtt, ndsd_ts
!***********************************************************************
procedure(DFT_FUNCTIONAL), pointer :: sub => null()
!     Sometime we need an external routine which covers something which
!     Libxc doesn't support.
procedure(DFT_FUNCTIONAL), pointer :: External_sub => null()
!                                                                      *
!***********************************************************************
! Global variable for MCPDFT functionals                               *
FLabel=KSDFA ! The user could be passing an explicit string! Hence, the local copy.

!
!     Set some flags and clean up the label to be just the label of the
!     underlaying DFT functional.
!
l_casdft = FLabel(1:2).eq.'T:' .or. FLabel(1:3).eq.'FT:'

lft      = FLabel(1:3).eq.'FT:'

If (l_casdft) Then
   If (lft) Then
      FLabel(:)=FLabel(4:)//'   '
   Else
      FLabel(:)=FLabel(3:)//'  '
   End If
   Do_MO=.true.
   Do_TwoEl=.true.
   If (.NOT.Do_PDFTPOT .and. .Not.DO_Grad) Only_exc=.True.
End If

If (FLabel(1:5)=='LDTF/')Then
   LDTF=.true.
   FLabel(:)=FLabel(6:)//'      '
End If
If (FLabel(1:5)=='NDSD/')Then
   NDSD=.true.
   FLabel(:)=FLabel(6:)//'      '
End If
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!      Default is to use the libxc interface
!      Coefficient for the individual contibutions are defaulted to 1.0D0

Sub => libxc_functionals     ! Default
Coeffs(:)=1.0D0              ! Default
!                                                                      *
!***********************************************************************
!                                                                      *
    Select Case(FLabel)
!                                                                      *
!***********************************************************************
!                                                                      *
!      LSDA LDA SVWN                                                   *
!                                                                      *
      Case('LSDA ','LDA ','SVWN ')
         Functional_type=LDA_type

!----    Slater exchange
!
!----    Vosko-Wilk-Nusair correlation functional III

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN_RPA]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_LDA_X,XC_LDA_C_VWN_RPA]
               End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_LDA_X,XC_LDA_C_VWN_RPA]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      LSDA5 LDA5 SVWN5                                                *
!                                                                      *
       Case('LSDA5','LDA5','SVWN5')
         Functional_type=LDA_type

!----    Slater exchange
!
!----    Vosko-Wilk-Nusair correlation functional V

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_LDA_X,XC_LDA_C_VWN]
              End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_LDA_X,XC_LDA_C_VWN]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFB                                                              *
!                                                                      *
       Case('HFB')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_B88]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFO                                                              *
!                                                                      *
       Case('HFO')
         Functional_type=GGA_type

!----    Slater exchange + OPTx exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_OPTX]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFG                                                              *
!                                                                      *
       Case('HFG')
         Functional_type=GGA_type

!----    Slater exchange + G96 exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_G96]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFB86                                                            *
!                                                                      *
       Case('HFB86')
         Functional_type=GGA_type

!----    Slater exchange + B86 exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_B86]
!                                                                      *
!***********************************************************************
!                                                                      *
!      HFS                                                             *
!                                                                      *
       Case('HFS')
         Functional_type=LDA_type

!----    Slater exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_LDA_X]
!                                                                      *
!***********************************************************************
!                                                                      *
!      XALPHA                                                          *
!                                                                      *
       Case('XALPHA')
         Functional_type=LDA_type

!----    Slater exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_LDA_C_XALPHA]
!        func_id(1:nFuncs)=[XC_LDA_X]
!        Coeffs(1)=0.70D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     Overlap                                                          *
!                                                                      *
      Case('Overlap')
         Functional_type=LDA_type
         Sub => Overlap
!                                                                      *
!***********************************************************************
!                                                                      *
!     NucAtt                                                           *
!                                                                      *
      Case('NucAtt')
         Functional_type=LDA_type
         Sub => NucAtt
!                                                                      *
!***********************************************************************
!                                                                      *
!     BWIG                                                             *
!                                                                      *
      Case('BWIG')
         Functional_type=GGA_type

!----    Slater exchange + B88 exchange
!
!----    Wigner correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B88,XC_LDA_C_OW_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     BLYP                                                             *
!                                                                      *
      Case('BLYP')
         Functional_type=GGA_type

!----    Slater exchange + B88 exchange
!
!----    Lee-Yang-Parr correlation

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_LYP]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_GGA_X_B88,XC_GGA_C_LYP]
               End If
            Else If (NDSD) Then
               If (KEOnly) Then
                  Sub => ndsd_ts
               Else
                  Only_exc=.True.
                  nFuncs=2
                  func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
                  External_Sub => ndsd_ts
               End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     OLYP                                                             *
!                                                                      *
      Case ('OLYP')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_OPTX,XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     KT3                                                              *
!                                                                      *
      Case('KT3')
         Functional_type=GGA_type

!----    Slater exchange + + OPT exchange + Keal-Tozer exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_XC_KT3]

!                                                                      *
!***********************************************************************
!                                                                      *
!     KT2                                                              *
!                                                                      *
      Case('KT2')
         Functional_type=GGA_type

!----    Slater exchange + Keal-Tozer exchange
!
!----    Vosko-Wilk-Nusair correlation functional III

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_XC_KT2]
!                                                                      *
!***********************************************************************
!                                                                      *
!     GLYP                                                             *
!                                                                      *
      Case('GLYP')
         Functional_type=GGA_type

!----    Slater exchange + Gill 96 exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_G96,XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     B86LYP                                                           *
!                                                                      *
      Case('B86LYP')
         Functional_type=GGA_type

!----    Slater exchange + Becke 86 exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B86,XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     BPBE                                                             *
!                                                                      *
      Case('BPBE')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     OreLYP
!                                                                      *
      Case('ORELYP')
         Functional_type=GGA_type

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_OPTX,XC_GGA_C_TM_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     OPBE                                                             *
!                                                                      *
      Case('OPBE')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_OPTX,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     GPBE                                                             *
!                                                                      *
      Case('GPBE')
         Functional_type=GGA_type

!----    Slater exchange + Gill 96 exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_G96,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     B86PBE                                                           *
!                                                                      *
      Case('B86PBE')
         Functional_type=GGA_type

!----    Slater exchange + Becke 86 exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B86,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     TLYP                                                             *
!                                                                      *
      Case('TLYP')
         Functional_type=GGA_type

!----    Lee-Yang-Parr correlation

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     B3LYP                                                            *
!                                                                      *
      Case('B3LYP ')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Perdew-Burk-Ernzerhof correlation
!
!----    Vosko-Wilk-Nusair correlation functional III
!
!----    Lee-Yang-Parr correlation

         nFuncs=1
         func_id(1:nFuncs)=[XC_HYB_GGA_XC_B3LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     O3LYP                                                            *
!                                                                      *
      Case('O3LYP ')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    Perdew-Burk-Ernzerhof correlation
!
!----    Vosko-Wilk-Nusair correlation functional III
!
!----    Lee-Yang-Parr correlation

         nFuncs=1
         func_id(1:nFuncs)=[XC_HYB_GGA_XC_O3LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     B2PLYP                                                           *
!                                                                      *
      Case('B2PLYP')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
         Coeffs(1)=0.470D0
         Coeffs(2)=0.73D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     O2PLYP                                                           *
!                                                                      *
      Case('O2PLYP')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    OPT exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=3
         func_id(1:nFuncs)=[XC_LDA_X,XC_GGA_X_OPTX,XC_GGA_C_LYP]
         Coeffs(1)=(0.525755D0-1.051510d0*(0.715845D0/1.431690d0))
         Coeffs(2)=(0.715845D0/1.431690d0)
         Coeffs(3)=0.75D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     B3LYP5                                                           *
!                                                                      *
      Case('B3LYP5')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Perdew-Burk-Ernzerhof correlation
!
!----    Vosko-Wilk-Nusair correlation functional V
!
!----    Lee-Yang-Parr correlation

         nFuncs=1
         func_id(1:nFuncs)=[XC_HYB_GGA_XC_B3LYP5]
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBE                                                              *
!                                                                      *
      Case('PBE ')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_PBE]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_GGA_X_PBE,XC_GGA_C_PBE]
               End If
            Else If (NDSD) Then
               If (KEOnly) Then
                  Sub => ndsd_ts
               Else
                  Only_exc=.True.
                  nFuncs=2
                  func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
                  External_Sub => ndsd_ts
               End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     revPBE                                                           *
!                                                                      *
      Case('REVPBE')
         Functional_type=GGA_type

!----    Revised Perdew-Burk-Ernzerhof exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE_R,XC_GGA_C_PBE]

!                                                                      *
!***********************************************************************
!                                                                      *
!     SSBSW                                                              *
!                                                                      *
      Case('SSBSW')
         Functional_type=GGA_type

!----    Swarts-Solo-Bickelhaupt sw exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_SSB_SW,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     SSBD                                                             *
!                                                                      *
      Case('SSBD')
         Functional_type=GGA_type

!----    Swarts-Solo-Bickelhaupt-D exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_SSB_D,XC_GGA_C_PBE]
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
!     S12H                                                             *
!                                                                      *
      Case('S12H')
         Functional_type=GGA_type

!----    Swarts hybrid exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_GGA_X_S12H,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     S12G                                                             *
!                                                                      *
      Case('S12G')
         Functional_type=GGA_type

!----    Swarts GGA exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_S12G,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBEsol                                                           *
!                                                                      *
      Case('PBESOL')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof SOL exchange
!
!----    Perdew-Burk-Ernzerhof SOL correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE_SOL,XC_GGA_C_PBE_SOL]
!                                                                      *
!***********************************************************************
!                                                                      *
!     RGE2                                                             *
!                                                                      *
      Case('RGE2')
         Functional_type=GGA_type

!----    Ruzsinzsky-Csonka-Scusceria exchange
!
!----    Perdew-Burk-Ernzerhof SOL correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_RGE2,XC_GGA_C_PBE_SOL]
!                                                                      *
!***********************************************************************
!                                                                      *
!     PTCA                                                             *
!                                                                      *
      Case('PTCA')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof exchange
!
!----    Tognotti-Cortona-Adamo correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_TCA]
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBE0                                                             *
!                                                                      *
      Case('PBE0')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
         Coeffs(1)=0.75D0
!        Coeffs(2)=1,00D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-L                                                            *
!                                                                      *
      Case('M06L')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 L exchange
!
!----    Minnesota 2006 L correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_MGGA_X_M06_L,XC_MGGA_C_M06_L]
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06                                                              *
!                                                                      *
      Case('M06 ')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 exchange
!
!----    Minnesota 2006 correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_MGGA_X_M06,XC_MGGA_C_M06]
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-2X                                                           *
!                                                                      *
      Case('M062X')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 2X exchange
!
!----    Minnesota 2006 2X correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_MGGA_X_M06_2X,XC_MGGA_C_M06_2X]
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-HF                                                           *
!                                                                      *
      Case('M06HF')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 HF exchange
!
!----    Minnesota 2006 HF correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_MGGA_X_M06_HF,XC_MGGA_C_M06_HF]
!                                                                      *
!***********************************************************************
!                                                                      *
!     BR3P86, DOI:10.1016/0009-2614(95)01143-2                         *
!                                                                      *
      Case('BR3P86')
         Functional_type=meta_GGA_type2

         A = 0.22D0
         B = 0.67D0
         C = 0.85D0

         nFuncs=4
         func_id(1:nFuncs)=[XC_LDA_X,XC_MGGA_X_BR89_1,XC_LDA_C_VWN_RPA,XC_GGA_C_PBE]
         Coeffs(1)=One-A-B
         Coeffs(2)=B
         Coeffs(3)=One
         Coeffs(4)=C
!                                                                      *
!***********************************************************************
!                                                                      *
!      Kinetic only  (Thomas-Fermi)                                    *
!                                                                      *
       Case('TF_only')
         Functional_type=LDA_type

         nFuncs=1
         func_id(1:nFuncs)=[XC_LDA_K_TF]
!                                                                      *
!***********************************************************************
!                                                                      *
!      HUNTER  (von Weizsacker KE, no calc of potential)               *
!                                                                      *
       Case('HUNTER')
         Functional_type=GGA_type

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_K_TFVW]
         Only_exc=.True.
!                                                                      *
!***********************************************************************
!                                                                      *
      Case default
         Call WarningMessage(2,' Driver: Undefined functional type!')
         Write (6,*) '         Functional=',FLabel(1:LEN(FLabel))
         Call Quit_OnUserError()
       End Select
!                                                                      *
!***********************************************************************
!                                                                      *
       If (Functional_type/=LDA_type.and.Functional_type/=GGA_type.and.l_CasDFT) Then
          Write (6,*) ' MC-PDFT combined with invalid functional class'
          Call Abend()
       End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     Now let's do some integration!
!     If the libxc interface is used do the proper initialization and closure.

      If (Associated(Sub,libxc_functionals)) Call Initiate_libxc_functionals(nD)

      Call DrvNQ(Sub,F_DFT,nD,Func,D_DS,nh1,nD,                         &
     &           Do_Grad,Grad,nGrad,Do_MO,Do_TwoEl,DFTFOCK)

      If (Associated(Sub,libxc_functionals)) Call Remove_libxc_functionals()

      If (Associated(External_Sub)) Call DrvNQ(External_Sub,F_DFT,nD,Func,D_DS,nh1,nD,         &
     &           Do_Grad,Grad,nGrad,Do_MO,Do_TwoEl,DFTFOCK)

      Sub          => Null()
      External_Sub => Null()
      Only_exc=.False.
      LDTF=.False.
      NDSD=.False.
!                                                                      *
!***********************************************************************
!                                                                      *
      End Subroutine Driver
