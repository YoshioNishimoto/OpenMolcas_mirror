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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine put array double data to the runfile.                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_dArray
!
!> @brief
!>   Add/update array data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put array data of type
!> ``Real*8`` into the runfile. The data items are
!> identified by the \p label. Below is a list of the
!> data items that are recognized. The labels are
!> case insensitive and significant to 16 characters.
!> (May change to 24 characters.)
!>
!> @warning
!> Naming convention is under development
!> and labels may change to the next version.
!>
!> For development purposes you can use an unsupported
!> label. Whenever such a field is accessed a warning
!> message is printed in the output, to remind the
!> developer to update this routine.
!>
!> List of known labels:
!>
!> - '``Analytic Hessian``'         Analytic Hessian.
!> - '``Center of Charge``'         Nuclear center of charge.
!> - '``Center of Mass``'           Nuclear center of mass.
!> - '``CMO_ab``'
!> - '``D1ao``'                     One particle density matrix, AO basis.
!> - '``D1ao_ab``'
!> - '``D1aoVar``'                  Generalized one particle density matrix, AO basis.
!> - '``D1av``'                     Average one particle density matrix, AO basis.
!> - '``D1mo``'                     One particle density matrix, MO basis.
!> - '``D1sao``'                    One particle spin density matrix, AO basis.
!> - '``D1activeao``'               One particle density matrix, AO basis, active orbitals
!> - '``D2av``'                     Average two particle density matrix for the active space, AO basis.
!> - '``dExcdRa``'                  The potential of the exchange-correlation functional.
!> - '``DLAO``'
!> - '``DLMO``'
!> - '``Effective nuclear Charge``' Effective nuclear charge for each unique atom.
!> - '``FockO_ab``'
!> - '``FockOcc``'                  Generalized Fock matrix, AO basis.
!> - '``GeoNew``'                   Next guess for Cartesian coordinates for the unique atoms.
!> - '``GeoNewPC``'                 Next guess for Cartesian coordinates for the unique point charges.
!> - '``GRAD``'                     Gradient with respect to nuclear displacements, for all unique atoms.
!> - '``Hess``'
!> - '``HF-forces``'                Hellmann--Feynman forces.
!> - '``Last orbitals``'            Last set of orbital computed.
!> - '``LCMO``'
!> - '``MEP-Coor``'                 List of nuclear coordinates along a minimum energy path, for unique atoms.
!> - '``MEP-Energies``'             List of energies along a minimum energy path.
!> - '``MEP-Grad``'                 List of nuclear gradients along a minimum energy path.
!> - '``MP2 restart``'              Information for restarting direct MP2 calculations.
!> - '``Mulliken Charge``'          Mulliken population charges for each unique atom.
!> - '``NEMO TPC``'
!> - '``Nuclear charge``'           Actual nuclear charge for each unique atom.
!> - '``OrbE``'                     SCF orbital energies.
!> - '``OrbE_ab``'
!> - '``P2MO``'
!> - '``PCM Charges``'              Charges for each tessera in the PCM model.
!> - '``PCM Info``'                 Misc. information needed for the PCM model.
!> - '``PLMO``'
!> - '``RASSCF orbitals``'          Last orbitals generated by the RASSCF module.
!> - '``Reaction field``'           Misc. information needed for the Kirkwood model.
!> - '``SCFInfoR``'                 Misc. information needed by the SCF module.
!> - '``SCF orbitals``'             Last orbitals generated by the SCF module.
!> - '``Slapaf Info 2``'            Misc. information needed by the SLAPAF module.
!> - '``Unique Coordinates``'       Cartesian coordinates for the symmetry unique atoms.
!> - '``Last energies``'            Energies for all roots in the last calculation.
!> - '``Dipole moment``'            The last computed dipole moment.
!> - '``MkNemo.vDisp``'             The displacements matrix as specified in the mknemo module.
!> - '``MkNemo.tqCluster``'         The transformation matrix for clusters as specified in the mknemo module.
!> - '``MkNemo.Energies``'          The energies of super-system and clusters as specified in the mknemo module.
!> - '``RASSCF OrbE``'              RASSCF orbital energies.
!> - '``GRD1``'                     MR-CISD gradient state1.
!> - '``GRD2``'                     MR-CISD gradient state2.
!> - '``NADC``'                     MR-CISD NADC vector state1/state2.
!> - '``MR-CISD energy``'           MR-CISD energies state1,state2.
!> - '``NOSEHOOVER``'               Extra-degrees of fredom needed to generated canonical ensemble.
!> - '``T-Matrix``'                 T-Matrix associated with PCO.
!> - '``rInt0``'                    Stored constrained values.
!> - '``Weights``'                  Weights used for alignment and hypersphere constraint.
!> - '``MEP-Lengths``'              Lengths of the MEP steps.
!> - '``MEP-Curvatures``'           Curvatures of the MEP steps.
!> - '``D1ao-``'                    Antisymmetric transition density matrix, in AO
!> - '``P2MOT``'                    "Fake" two-body density needed for MC-PDFT gradient calculations
!>
!> @param[in] Label Name of field
!> @param[in] Data  Data to put on runfile
!> @param[in] nData Length of array
!***********************************************************************
      Subroutine Put_dArray(Label,Data,nData)
      Implicit None
#include "pg_da_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
      Character*(*) Label
      Integer       nData
      Real*8        Data(nData)
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
      Character*16 RecLab(nTocDA)
      Integer      RecIdx(nTocDA)
      Integer      RecLen(nTocDA)
      Save         RecLab
      Save         RecIdx
      Save         RecLen
!
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      nTmp
      Integer      item
      Integer      iTmp
      Integer      i
!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
      Call ffRun('dArray labels',nTmp,iTmp)
      If(nTmp.eq.0) Then
         Do i=1,nTocDA
            RecLab(i)=' '
            RecIdx(i)=sNotUsed
            RecLen(i)=0
         End Do
!
!        Observe that label is at most 16 characters!
!
!                     1234567890123456
         RecLab(  1)='Analytic Hessian'
         RecLab(  2)='Center of Charge'
         RecLab(  3)='Center of Mass  '
         RecLab(  4)='CMO_ab          '
         RecLab(  5)='D1ao            '
         RecLab(  6)='D1ao_ab         '
         RecLab(  7)='D1aoVar         '
         RecLab(  8)='D1av            '
         RecLab(  9)='D1mo            '
         RecLab( 10)='D1sao           '
         RecLab( 11)='D2av            '
         RecLab( 12)='dExcdRa         '
         RecLab( 13)='DLAO            '
         RecLab( 14)='DLMO            '
         RecLab( 15)='Effective nuclea' !r charge
         RecLab( 16)='FockO_ab        '
         RecLab( 17)='FockOcc         '
         RecLab( 18)='GeoNew          '
         RecLab( 19)='GeoNewPC        '
         RecLab( 20)='GRAD            '
         RecLab( 21)='Hess            '
         RecLab( 22)='HF-forces       '
         RecLab( 23)='Last orbitals   '
         RecLab( 24)='LCMO            '
         RecLab( 25)='MEP-Coor        '
         RecLab( 26)='MEP-Energies    '
         RecLab( 27)='MEP-Grad        '
         RecLab( 28)='MP2 restart     '
         RecLab( 29)='Mulliken Charge '
         RecLab( 30)='NEMO TPC        '
         RecLab( 31)='Nuclear charge  '
         RecLab( 32)='OrbE            '
         RecLab( 33)='OrbE_ab         '
         RecLab( 34)='P2MO            '
         RecLab( 35)='PCM Charges     '
         RecLab( 36)='PCM Info        '
         RecLab( 37)='PLMO            '
         RecLab( 38)='RASSCF orbitals '
         RecLab( 39)='Reaction field  '
         RecLab( 40)='SCFInfoR        '
         RecLab( 41)='SCF orbitals    '
         RecLab( 42)='Slapaf Info 2   '
         RecLab( 43)='Unique Coordinat' !es
         RecLab( 44)='Vxc_ref         '
! mess started here :)
         RecLab( 45)='PotNuc00        '
         RecLab( 46)='h1_raw          '
         RecLab( 47)='h1    XX        '
         RecLab( 48)='HEFF            '
         RecLab( 49)='PotNucXX        '
         RecLab( 50)='Quad_r          '
         RecLab( 51)='RCTFLD          '
         RecLab( 52)='RFrInfo         '
         RecLab( 53)='                '
         RecLab( 54)='DKH_Info        '
         RecLab( 55)='Real_Info       '
         RecLab( 56)='Last orbitals_ab'
         RecLab( 57)='SCFInfoI_ab     '
         RecLab( 58)='SCFInfoR_ab     '
         RecLab( 59)='Transverse      '
         RecLab( 60)='SM              '
         RecLab( 61)='LP_Coor         '
         RecLab( 62)='LP_Q            '
         RecLab( 63)='DFT_TwoEl       '
         RecLab( 64)='Unit Cell Vector'
         RecLab( 65)='SCF orbitals_ab '
         RecLab( 66)='Guessorb        '
         RecLab( 67)='Guessorb energie' !s
         RecLab( 68)='Last energies   '
         RecLab( 69)='LoProp Dens 0   '
         RecLab( 70)='LoProp Dens 1   '
         RecLab( 71)='LoProp Dens 2   '
         RecLab( 72)='LoProp Dens 3   '
         RecLab( 73)='LoProp Dens 4   '
         RecLab( 74)='LoProp Dens 5   '
         RecLab( 75)='LoProp Dens 6   '
         RecLab( 76)='LoProp Integrals'
         RecLab( 77)='MpProp Orb Ener '
         RecLab( 78)='LoProp H0       '
         RecLab( 79)='Dipole moment   '
         RecLab( 80)='RICD_Info       '
         RecLab( 81)='BMtrx           '
         RecLab( 82)='CList           '
         RecLab( 83)='DList           '
         RecLab( 84)='                ' ! Unused
         RecLab( 85)='MkNemo.vDisp    '
         RecLab( 86)='MkNemo.tqCluster'
         RecLab( 87)='MkNemo.Energies '
         RecLab( 88)='MMHessian       '
         RecLab( 89)='Bfn Coordinates '
         RecLab( 90)='Pseudo Coordinat' !es
         RecLab( 91)='Pseudo Charge   '
         RecLab( 92)='RASSCF OrbE     '
         RecLab( 93)='Ref_Geom        '
         RecLab( 94)='LoProp Charge   '
         RecLab( 95)='Initial Coordina' !tes
         RecLab( 96)='Grad State1     '
         RecLab( 97)='Grad State2     '
         RecLab( 98)='NADC            '
         RecLab( 99)='MR-CISD energy  '
         RecLab(100)='Saddle          '
         RecLab(101)='Reaction Vector '
         RecLab(102)='IRC-Coor        '
         RecLab(103)='IRC-Energies    '
         RecLab(104)='IRC-Grad        '
         RecLab(105)='MM Grad         '
         RecLab(106)='Velocities      '
         RecLab(107)='FC-Matrix       '
         RecLab(108)='umass           '
         RecLab(109)='ESO_SINGLE      '
         RecLab(110)='UMATR_SINGLE    '
         RecLab(111)='UMATI_SINGLE    '
         RecLab(112)='ANGM_SINGLE     '
         RecLab(113)='TanVec          '
         RecLab(114)='Nuc Potential   '
         RecLab(115)='RF CASSCF Vector'
         RecLab(116)='Cholesky BkmThr '
         RecLab(117)='NOSEHOOVER      '
         RecLab(118)='T-Matrix        '
         RecLab(119)='rInt0           '
         RecLab(120)='Weights         '
         RecLab(121)='MEP-Lengths     '
         RecLab(122)='MEP-Curvatures  '
         RecLab(123)='Hss_X           '
         RecLab(124)='Hss_Q           '
         RecLab(125)='KtB             '
         RecLab(126)='BMxOld          '
         RecLab(127)='TROld           '
         RecLab(128)='qInt            '
         RecLab(129)='dqInt           '
         RecLab(130)='Fragment_Fock   '
         RecLab(131)='RAmatrixV       '
         RecLab(132)='IAmatrixV       '
         RecLab(133)='AllCIP          '
         RecLab(134)='AllCIPP         '
         RecLab(135)='VenergyP        '
         RecLab(136)='K               '
         RecLab(137)='MMO Coords      '
         RecLab(138)='MMO Grad        '
         RecLab(139)='Hss_upd         '
         RecLab(140)='TR              '
         RecLab(141)='D1ao-           '
         RecLab(142)='ESFS_SINGLE     '
         RecLab(143)='LA Fact         '
         RecLab(144)='primitives      '
         RecLab(145)='Isotopes        ' ! Atom masses in a.u.
         RecLab(146)='P2AO            '
         RecLab(147)='State Overlaps  '
         RecLab(148)='EFP_Coors       ' ! EFP fragment coordinates
         RecLab(149)='DIP1_SINGLE     '
         RecLab(150)='P2MOT           '
         RecLab(151)='ONTOPO          '
         RecLab(152)='ONTOPT          '
         RecLab(153)='OE_OT           '
         RecLab(154)='TEG_OT          '
         RecLab(155)='FI_V            '
         RecLab(156)='FA_V            '
         RecLab(157)='FOCK_PDFT       '
         RecLab(158)='AMFI_SINGLE     '
         RecLab(159)='HAMSOR_SINGLE   '
         RecLab(160)='HAMSOI_SINGLE   '
         RecLab(161)='Last Dipole Mome' !nts
         RecLab(162)='Un_cen Effective' ! Charge
         RecLab(163)='Un_cen Coordinat' !es
         RecLab(164)='ANGMR_NSS       '
         RecLab(165)='ANGMI_NSS       '
         RecLab(166)='EDIPR_NSS       '
         RecLab(167)='EDIPI_NSS       '
         RecLab(168)='SPINR_NSS       '
         RecLab(169)='SPINI_NSS       '
         RecLab(170)='Proj_Coord      '
         RecLab(171)='d1activeao      '
         RecLab(172)='Keep_Coord      '
         RecLab(173)='PCMSph          '
         RecLab(174)='PCMTess         '
         RecLab(175)='Vert            '
         RecLab(176)='Centr           '
         RecLab(177)='SSph            '
         RecLab(178)='PCMDM           '
         RecLab(179)='EF_Centers      '
         RecLab(180)='OAM_Center      '
         RecLab(181)='OMQ_Center      '
         RecLab(182)='DMS_Centers     '
         RecLab(183)='Wel_Info        '
         RecLab(184)='AMP_Center      '
         RecLab(185)='RP_Centers      '
         RecLab(186)='XF              '
         RecLab(187)='rDmp            '
         RecLab(188)='rDmp:A          '
         RecLab(189)='rDmp:S          '
         RecLab(190)='D1saoVar        '
         Reclab(191)='ESFS_SINGLEAU   '
         Reclab(192)='ESO_LOW         '
         Reclab(193)='SFS_HAM         '
         Reclab(194)='SFS_OVLP        '
!Arrays 195-207 for computing MS-PDFT gradient
         RecLab(195)='FocMS           '
         RecLab(196)='MSPDFTD5        '
         RecLab(197)='MSPDFTD6        '
         RecLab(198)='TwoEIntegral    '
         RecLab(199)='D1MOt           '
         RecLab(200)='D1INTER         '
         RecLab(201)='P2INTER         '
         RecLab(202)='D1AO_MS         '
         RecLab(203)='D1SAO_MS        '
         RecLab(204)='MS_FINAL_ROT    '
         RecLab(205)='F1MS            '
         RecLab(206)='F2MS            '
         RecLab(207)='FxyMS           '
         RecLab(208)='SH_Ovlp_Save    '
!                     1234567890123456
!
!        If you go beyond 256: update pg_da_info.fh and this line!
         Call cWrRun('dArray labels',RecLab,16*nTocDA)
         Call iWrRun('dArray indices',RecIdx,nTocDA)
         Call iWrRun('dArray lengths',RecLen,nTocDA)
      Else
         Call cRdRun('dArray labels',RecLab,16*nTocDA)
         Call iRdRun('dArray indices',RecIdx,nTocDA)
         Call iRdRun('dArray lengths',RecLen,nTocDA)
      End If
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocDA
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
!
! Do we create a new temporary field?
!
      If(item.eq.-1) Then
         Do i=1,nTocDA
            If(RecLab(i).eq.' ') item=i
         End Do
         If(item.ne.-1) Then
            RecLab(item)=Label
            RecIdx(item)=sSpecialField
            Call cWrRun('dArray labels',RecLab,16*nTocDA)
            Call iWrRun('dArray indices',RecIdx,nTocDA)
         End If
      End If
!
! Is this a temporary field?
!
      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, writing temporary dArray field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
#ifndef _DEVEL_
            Call AbEnd()
#endif
         End If
      End If
!----------------------------------------------------------------------*
! Write data to disk.                                                  *
!----------------------------------------------------------------------*
      If(item.eq.-1) Then
         Call SysAbendMsg('put_dArray','Could not locate',Label)
      End If
      Call dWrRun(RecLab(item),Data,ndata)
!     Write(6,*) Data
      If(RecIdx(item).eq.0) Then
         RecIdx(item)=sRegularField
         Call iWrRun('dArray indices',RecIdx,nTocDA)
      End If
      If(RecLen(item).ne.nData) Then
         RecLen(item)=nData
         Call iWrRun('dArray lengths',RecLen,nTocDA)
      End If
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
      Return
      End
