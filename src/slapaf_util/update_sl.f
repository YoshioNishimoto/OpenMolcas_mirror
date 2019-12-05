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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine Update_sl(iter,MaxItr,NmIter,iInt,nFix,nInter,qInt,
     &                     Shift,
     &                     Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,
     &                     Energy,UpMeth,ed,Line_Search,Step_Trunc,
     &                     nLambda,iRow_c,nsAtom,AtomLbl,nSym,iOper,
     &                     mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                     rLambda,ipCx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                     iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                     nWndw,Mode,ipMF,
     &                     iOptH,HUpMet,kIter,GNrm_Threshold,IRC,
     &                     dMass,HrmFrq_Show,CnstWght,Curvilinear,
     &                     Redundant,Degen,ThrEne,ThrGrd)
************************************************************************
*                                                                      *
*     Object: to update coordinates                                    *
*                                                                      *
*    Input:                                                            *
*      iter           : iteration counter                              *
*      MaxItr         : max number of iteration                        *
*      NmIter         : number of iteration in numerical approach      *
*      iInt           : number of internal coordinates to vary         *
*      nFix           : number of frozen internal coordinates          *
*      nInter         : total number of internal coordinates           *
*      qInt(*,iter )  : the internal coordinates                       *
*      Grad(*,iter )  : the gradient in the internal coordinates       *
*      iOptC          : option flag for update methods                 *
*      Beta           : damping factor                                 *
*      Lbl            : character labels for internal coordinates      *
*      nLbl           : length of Lbl                                  *
*      GNrm           : the norm of the gradient in each iteration     *
*      Energy         : the energy of each iteration                   *
*      Line_Search    : logical flag for line search                   *
*      nLambda        : number of contraints                           *
*      iRow_c         : number of lines on the UDC file                *
*      nsAtom         : number of symmetry unique atoms                *
*      AtomLbl        : character string with atom labels              *
*      nSym           : number of irreps                               *
*      iOper          : integer representations of symmetry operators  *
*      mxdc           : max number of nsAtom                           *
*      jStab          : integer list of stabilizers                    *
*      nStab          : number of stabilizers                          *
*      BMx            : the so-called Wilson B matrix                  *
*      Smmtrc         : logical flag for symmetry properties           *
*      nDimBC         : dimension of redundant coordinates(?)          *
*      rLambda        : vector for Lagrange multipliers                *
*      ipCx           : pointer to cartesian coordinates               *
*      iNeg           : Hessian index                                  *
*      Labels         : character string of primitive int. coord.      *
*      nLabels        : length of Labels                               *
*      CnstWght       : constraints weight                             *
*                                                                      *
*    OutPut:                                                           *
*      Shift(*,iter ) : the shift of the internal coordinates          *
*      qInt(*,iter +1): the internal coordinates to be used in the     *
*                       next iteration                                 *
*      UpMeth         : character label with update method abrivation  *
*      ed             : estimated energy change to the next point      *
*      Step_Trunc     : character label to denote truncation type      *
*      GrdMax         : largest gradient component                     *
*      StpMax         : largest step component                         *
*      GrdLbl         : character label of component with GrdMax       *
*      StpLbl         : character label of component with StpLbl       *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh                                             *
*             2000                                                     *
************************************************************************
      Use NewH_mod
      Use AI, only: Kriging, miAI, meAI, nspAI
      Implicit Real*8 (a-h,o-z)
      External Restriction_Step, Restriction_Dispersion
      Real*8 Restriction_Step, Restriction_Dispersion
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "Molcas.fh"
      Real*8 qInt(nInter,MaxItr), Shift(nInter,MaxItr),
     &       Grad(nInter,MaxItr), GNrm(MaxItr), Energy(MaxItr),
     &       BMx(3*nsAtom,3*nsAtom), rLambda(nLambda,MaxItr),
     &       dMass(nsAtom), Degen(3*nsAtom), dEner
      Integer iOper(0:nSym-1), jStab(0:7,nsAtom), nStab(nsAtom),
     &        iNeg(2)
      Logical Line_Search, Smmtrc(3*nsAtom),
     &        FindTS, TSC, HrmFrq_Show, Curvilinear, Redundant
      Character Lbl(nLbl)*8, GrdLbl*8, StpLbl*8, Step_Trunc,
     &          Labels(nLabels)*8, AtomLbl(nsAtom)*(LENIN), UpMeth*6,
     &          HUpMet*6
*
      Logical Kriging_Hessian, Not_Converged, Single_l_value
      Real*8, Allocatable:: Array_l(:)
      Real*8, Allocatable:: Energy_s(:), qInt_s(:,:), Grad_s(:,:)
*
      iRout=153
      iPrint=nPrint(iRout)
      Call QEnter('Update')
*
      If (iPrint.ge.99) Then
         Call RecPrt('Update: qInt',' ',qInt,nInter,Iter)
         Call RecPrt('Update: Shift',' ',Shift,nInter,Iter-1)
         Call RecPrt('Update: GNrm',' ',GNrm,Iter,1)
      End If
*define UNIT_MM
#ifdef UNIT_MM
*
*---- If redundant Cartesians and no symmetry, use unit matrix for MM atoms
*
      If (Redundant.and.(.not.Curvilinear).and.
     &    (3*nsAtom.eq.nInter)) Then
         Call mma_allocate(UpdMask,nInter,label="UpdMask")
         Call MMCount(nsAtom,nAtMM,ipIsMM)
         Do iAtom=1,nsAtom
            If (iWork(ipIsMM+iAtom-1).eq.1) Then
               Do i=1,3
                 UpdMask((iAtom-1)*3+i)=1
               End Do
            Else
               Do i=1,3
                 UpdMask((iAtom-1)*3+i)=0
               End Do
            End If
         End Do
         Call GetMem('IsMM for atoms','Free','Inte',ipIsMM,nsAtom)
      End If
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Redundant)
#endif
*
      iOpt_RS=0
      Kriging_Hessian =.FALSE.
      If (iter.eq.NmIter.and.NmIter.ne.1) Then
*
*------- On the first iteration after a numerical evaluation of the
*        Hessian we like the step to be relative to the initial
*        structure.
*
         If (iPrint.ge.99) Write(6,*)'UpDate_SL: first iteration'
         iter_=1
         Call GetMem('t_Shift','Allo','Real',iptmp1,nInter*iter_)
         Call GetMem('t_qInt ','Allo','Real',iptmp2,nInter*(iter_+1))
*
         call dcopy_(nInter,qInt,1,Work(iptmp2),1)
*
         Call Update_sl_(iter_,iInt,nFix,nInter,Work(iptmp2),
     &                   Work(iptmp1),Grad,iOptC,Beta,Beta_Disp,
     &                   Lbl,GNrm,Energy,UpMeth,ed,Line_Search,
     &                   Step_Trunc,nLambda,iRow_c,nsAtom,AtomLbl,nSym,
     &                   iOper,mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                   rLambda,ipCx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                   iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                   nWndw,Mode,ipMF,
     &                   iOptH,HUpMet,kIter,GNrm_Threshold,IRC,dMass,
     &                   HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                   Kriging_Hessian,qBeta,Restriction_step,
     &                   iOpt_RS)
*
*------- Move new coordinates to the correct position and compute the
*        corresponding shift.
*
         iOffq=iptmp2+nInter
         call dcopy_(nInter,Work(iOffq),1,qInt(1,iter+1),1)
         call dcopy_(nInter,Work(iOffq),1,Shift(1,iter),1)
         Call DaXpY_(nInter,-One,qInt(1,iter),1,Shift(1,iter),1)
*
         Call GetMem('t_qInt ','Free','Real',iptmp2,nInter*(iter_+1))
         Call GetMem('t_Shift','Free','Real',iptmp1,nInter*iter_)
*
      Else
*        ------- AI loop begin here
*#define _DEBUG_
         If (Kriging .AND. iter.ge.nspAI) then
            Kriging_Hessian =.TRUE.
            iOpt_RS=1   ! Activate restricted variance.
            iterAI=iter
            dEner=meAI
            nRaw=Min(iter,nWndw/2)
            iFirst = iter - nRaw + 1
            iterK=0
            dqdq=Zero
            qBeta=Beta
#ifdef _DEBUG_
            Write (6,*) 'iFirst,nRaw=',iFirst,nRaw
            Call RecPrt('qInt(0)',  ' ',qInt(1,iFirst),nInter,nRaw)
*
            Call RecPrt('Energy(0)',' ',Energy(iFirst),1,nRaw)
*
            Call RecPrt('Grad(0)',  ' ',Grad(1,iFirst),nInter,nRaw)
*
            Call RecPrt('Shift',  ' ',Shift(1,iFirst),nInter,nRaw)
#endif
*
*           Pass the data point to the GEK routine. Remember to
*           change the sign of the gradients.
*
*           Note that we could have some kind of sorting here if we
*           like!
*
*define _UNSORTED_
#ifdef _UNSORTED_
            Write (6,*) 'Unsorted'
            Call DScal_(nInter*nRaw,-One,Grad(1,iFirst),1)
            Call Start_Kriging(nRaw,nInter,
     &                            qInt(1,iFirst),
     &                            Grad(1,iFirst),
     &                            Energy(iFirst))
            Call DScal_(nInter*nRaw,-One,Grad(1,iFirst),1)
#else
*
*           Sort the data so that it some in an order of the points
*           closest to the last point. Make sure that the reference
*           point is the last point such that is always will define
*           the reference point. This is improtant when the bias is
*           defined relative to the last points energy.
*
*           This code will have to be cleanup up later.
*
            ipCx_Ref=ipCx + (iter-1)*(3*nsAtom)
            Call mma_Allocate(Energy_s,nRaw,Label="Energy_s")
            Call mma_Allocate(qInt_s,nInter,nRaw,Label="qInt_s")
            Call mma_Allocate(Grad_s,nInter,nRaw,Label="Grad_s")
*
            Call DCopy_(nInter,qInt(1,iter),1,qInt_s(1,nRaw),1)
            Call DCopy_(nInter,Grad(1,iter),1,Grad_s(1,nRaw),1)
            Energy_s(nRaw)=Energy(iter)
*
*           Pick up the coordinates in descending order with the ones
*           that are the closest to the current structure.
*
            iSt=Max(1,iter-nWndw+1)
            Thr_low = 0.0D0
            Thr_high= 99.0D0
            Do iRaw = nRaw-1, 1, -1
*
               kter=-1
               Do jter = iSt, iter-1
*
*                 Compute the distance in Cartesian coordinates.
*
                  Distance=Zero
                  Do ix = 1, 3*nsAtom
                     Distance= Distance +
     &                         Degen(ix) *
     &                        (Work(ipCx_ref+ix-1) -
     &                         Work(ipCx + (jter-1)*3*nsAtom+ix-1))**2
                  End Do
                  Distance = sqrt(Distance)
*                 Write (*,*) 'jter,Distance=',jter,Distance
*
                  If (Distance.gt.Thr_low .and.
     &                Distance.lt.Thr_high) Then
                     kter=jter
                     Thr_high=Distance
                  End If
               End Do
               If (kter.eq.-1) Then
                  Write (6,*) 'kter not set!'
               Else
                  Call DCopy_(nInter,qInt(1,kter),1,qInt_s(1,iRaw),1)
                  Call DCopy_(nInter,Grad(1,kter),1,Grad_s(1,iRaw),1)
                  Energy_s(iRaw)=Energy(kter)
                  Thr_low=Thr_high
                  Thr_high= 99.0D0
               End If
            End Do
#ifdef _DEBUG_
            Call RecPrt('qInt_s(s)',  ' ',qInt_s,nInter,nRaw)
            Call RecPrt('Energy_s(s)',' ',Energy_s,1,nRaw)
            Call RecPrt('Grad_s(s)',  ' ',Grad_s,nInter,nRaw)
#endif
*
            Call DScal_(nInter*nRaw,-One,Grad_s,1)
            Call Start_Kriging(nRaw,nInter,
     &                            qInt_s,
     &                            Grad_s,
     &                            Energy_s)
*
            Call mma_deAllocate(Energy_s)
            Call mma_deAllocate(qInt_s)
            Call mma_deAllocate(Grad_s)
#endif
*
*           Update the l value dynamically. Here we compare the actual,
*           ab inito value with the GEK prediction of the gradient.
*
            Call Get_dScalar('Value_l',Value_l)
*
            If (iter.gt.nspAI) Then
               iOld=iFirst+nRaw-2
               xxx=DDot_(nInter,Grad(1,iOld),1,Grad(1,iOld),1)
               iNew=iFirst+nRaw-1
               yyy=DDot_(nInter,Grad(1,iNew),1,Grad(1,iNew),1)
*              Write (*,*) 'yyy,xxx=',yyy,xxx
               If (yyy.gt.xxx) Then
                  Value_l=Value_l * 0.95D0
               Else
                  Value_l=Value_l * 1.05D0
               End If
*
*              Update the restricted variance threshold.
*
*              Beta_Disp=Max(Abs(Energy(iNew)-Energy(iOld)),
*    &                       1.0D-6)
            End If
*           Write (6,*) 'Modified l value:',Value_l
*
*           Single_l_value=.True.
            Single_l_value=.False.
            If (Single_l_value) Then
               Call set_l_kriging([Value_l],1)
            Else
               Call mma_Allocate(Array_l,nInter,Label='Array_l')
*              Call DCopy_(nInter,[1.0D0],0,Array_l,1)
               Call Set_l_Array(Array_l,nInter)
*              Call RecPrt('l values',' ',Array_l,nInter,1)
*              Write (6,*) 'Value_l=',Value_l
               Call DScal_(nInter,Value_l,Array_l,1)
*              Call RecPrt('l values',' ',Array_l,nInter,1)
               Call set_l_kriging(Array_l,nInter)
               Call mma_DeAllocate(Array_l)
            End If
            Call Put_dScalar('Value_l',Value_l)
*
*
*           Start the Kriging loop.
*
            Not_Converged = .True.
            do while (Not_Converged)
               kIter_=iterAI
#ifdef _DEBUG_
               Write (6,*)
               Write (6,*) 'Do iterAI: ',iterAI
#endif
               Call Update_sl_(iterAI,iInt,nFix,nInter,
     &                qInt,Shift,Grad,iOptC,Beta,Beta_Disp,
     &                Lbl,GNrm,Energy,
     &                UpMeth,ed,Line_Search,Step_Trunc,nLambda,
     &                iRow_c,nsAtom,AtomLbl,nSym,iOper,mxdc,jStab,
     &                nStab,BMx,Smmtrc,nDimBC,rLambda,ipCx,
     &                GrdMax,StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &                Labels,nLabels,FindTS,TSC,nRowH,
     &                nWndw/2,Mode,ipMF,
     &                iOptH,HUpMet,kIter_,GNrm_Threshold,IRC,dMass,
     &                HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                Kriging_Hessian,qBeta,Restriction_dispersion,
     &                iOpt_RS)
*
*              Change lable of updating method if kriging points has
*              been used.
*
               If (iterK.gt.0) Then
                  UpMeth='GPR   '
                  Write (UpMeth(4:6),'(I3)') iterK
               Else
                  UpMeth='RV-RFO'
               End If
*
*              Compute the step length from the last ab inito point
*              to the most recent krining point.
*
               dqdq=Zero
*              Write (6,*) 'iterAI+1=',iterAI+1
*              Write (6,*) 'iFirst+nRaw-1=',iFirst+nRaw-1
               Do iInter=1,nInter
*                 Write (6,*) qInt(iInter,iterAI+1)
*    &                 ,qInt(iInter,iFirst+nRaw-1)
                  dqdq=(qInt(iInter,iterAI+1)
     &                 -qInt(iInter,iFirst+nRaw-1))**2
               End Do
               If (iterK.eq.0.and.Step_trunc.eq.'*') dqdq=qBeta**2
*              Write (6,*) 'dqdq=',dqdq
*              Write (6,*) 'qBeta**2=',qBeta**2
#ifdef _TEST1_
               If (iterK.gt.0.and.dqdq.gt.qBeta**2) Then
*
                  sh2=Zero
                  D2 =Zero
                  shD=Zero
                  Do iInter=1,nInter
                     sh2=sh2 + Shift(iInter,iterAI)**2
                     D2 =D2 + (qInt(iInter,iterAI)-
     &                         qInt(iInter,iFirst+nRaw-1))**2
                     shD=shD + Shift(iInter,iterAI)
     &                       *(qInt(iInter,iterAI)-
     &                         qInt(iInter,iFirst+nRaw-1))
                  End Do
*                 Write (6,*) 'Sh2,shD,D2=',sh2,shD,D2
                  D2 = D2 - qBeta**2
                  shD=Two*shD
                  D2=D2/sh2
                  shD=shD/sh2
*                 Write (6,*) 'shD,D2=',shD,D2
*
*                 compute the scaling factor
*
                  alpha=-shD/Two + Sqrt((shD/Two)**2 - D2)
*                 Write (6,*) 'alpha=',Alpha
*
*                 rescale the step
*
                  Do iInter=1,nInter
                     Shift(iInter,iterAI) = alpha*Shift(iInter,iterAI)
                     qInt(iInter,iterAI+1) = qInt(iInter,iterAI)
     &                                      + Shift(iInter,iterAI)
                  End Do
*
                  dqdq=qBeta**2
                  Step_trunc='*'
                  Write (6,*) ' Step has been scaled'
               End If
#endif

*              Compute the energy and gradient according to the
*              surrogate model for the new coordinates.
*
               Call Energy_Kriging(qInt(1,iterAI+1),Energy(iterAI+1),
     &                             nInter)
               Call Dispersion_Kriging(qInt(1,iterAI+1),Dummy,nInter)
               E_Disp = Dummy
*              Write (6,*) 'E_Disp=',E_disp
               Call Gradient_Kriging(qInt(1,iterAI+1),Grad(1,iterAI+1),
     &                               nInter)
               Call DScal_(nInter,-One,Grad(1,iterAI+1),1)
*
               iterK  = iterK  + 1
               iterAI = iterAI + 1
               dEner = Energy(iterAI) - Energy(iterAI-1)
*
#ifdef _DEBUG_
               Call RecPrt('qInt(x):',' ',qInt,nInter,iterAI)
               Call RecPrt('Ener(x):',' ',Energy,1,iterAI)
               Call RecPrt('Grad(x):',' ',Grad,nInter,iterAI)
#endif
*              Write (6,*) 'ThrEne,ThrGrd',ThrEne,ThrGrd
*              Write (6,*) 'dEner=',dEner
               If (ThrEne.gt.Zero) Then
                  Not_Converged = Abs(dEner).ge.ThrEne
                  Not_Converged = Not_Converged .and. iterK.lt.miAI
                  Not_Converged = Not_Converged .and. dqdq.lt.qBeta**2
               Else
                  FAbs=Sqrt(DDot_(nInter,Grad(1,iterAI),1,
     &                                  Grad(1,iterAI),1)/DBLE(nInter))
                  RMS =Sqrt(DDot_(nInter,Shift(1,iterAI-1),1,
     &                                Shift(1,iterAI-1),1)/DBLE(nInter))
                  GrdMx=Zero
                  RMSMx=Zero
                  Do iInter = 1, nInter
                     GrdMx=Max(GrdMx,Abs(Grad(iInter,iterAI)))
                     RMSMx=Max(RMSMx,Abs(Shift(iInter,iterAI-1)))
                  End Do
*
*                 Write (6,*) 'FAbs,GrdMx,RMS,RMSMx:',
*    &                         FAbs,GrdMx,RMS,RMSMx
                  Not_Converged = FAbs.gt.ThrGrd
*                 Write (6,*) FAbs.gt.ThrGrd
                  Not_Converged = Not_Converged .or.
     &                            GrdMx.gt.ThrGrd*OneHalf
*                 Write (6,*)     GrdMx.gt.ThrGrd*OneHalf
                  Not_Converged = Not_Converged .or.
     &                            RMS.gt.ThrGrd*Four
*                 Write (6,*)     RMS.gt.ThrGrd*Four
                  Not_Converged = Not_Converged .or.
     &                            RMSMx.gt.ThrGrd*Six
*                 Write (6,*)     RMSMx.gt.ThrGrd*Six
                  Not_Converged = Not_Converged .and. iterK.le.miAI
*                 Write (6,*)     iterK,miAI
*                 Write (6,*)     iterK.le.miAI
*                 Write (6,*)     'Not_Converged=',Not_Converged
                  If (Step_trunc.eq.'*') Not_Converged=.False.
               End If
            End Do  ! Do While
            If (Step_trunc.eq.'*') Then
               Call Get_dScalar('Value_l',Value_l)
               Value_l=Value_l * 0.95D0
*              Write (6,*) ' Set l value to:',Value_l
               Call Put_dScalar('Value_l',Value_l)
            End If
*
*           Save the optimized kriging coordinates as the coordinates
*           for the next macro iteration.
*
            Call DCopy_(nInter,qInt(1,iterAI),1,qInt(1,iter+1),1)
            Call DCopy_(nInter,qInt(1,iter+1),1,Shift(1,iter),1)
            Call DaXpY_(nInter,-One,qInt(1,iter),1,Shift(1,iter),1)
            Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter,
     &                  Grad(1,iter),Shift(1,iter),Lbl)
#ifdef _DEBUG_
            Call RecPrt('qInt(3):',' ',qInt,nInter,iter+1)
            Call RecPrt('Shift(3):',' ',Shift,nInter,iter)
#endif
*
*           write(6,*) 'finished do iter',iterAI
*           De allocating memory used by Kriging
            Call Finish_Kriging()
*        ------- AI loop ends here
         Else
            Call Update_sl_(iter,iInt,nFix,nInter,qInt,Shift,
     &                   Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,Energy,
     &                   UpMeth,ed,Line_Search,Step_Trunc,nLambda,
     &                   iRow_c,nsAtom,AtomLbl,nSym,iOper,mxdc,jStab,
     &                   nStab,BMx,Smmtrc,nDimBC,rLambda,ipCx,
     &                   GrdMax,StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &                   Labels,nLabels,FindTS,TSC,nRowH,
     &                   nWndw,Mode,ipMF,
     &                   iOptH,HUpMet,kIter,GNrm_Threshold,IRC,dMass,
     &                   HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                   Kriging_Hessian,qBeta,Restriction_step,
     &                   iOpt_RS)
         End If
      End If
*
*-----Write out the shift in internal coordinates basis.
*
      If (iPrint.ge.99) Then
         Call RecPrt(
     &      'Shifts in internal coordinate basis / au or rad',
     &      ' ',Shift,nInter,Iter)
         Call RecPrt(
     &      'qInt in internal coordinate basis / au or rad',
     &      ' ',qInt,nInter,Iter+1)
      End If
*
*---- Remove unneeded fields from the runfile
      Call Put_dArray('BMxOld',Work(ip_Dummy),0)
      Call Put_dArray('TROld',Work(ip_Dummy),0)
      If (Allocated(UpdMask)) Call mma_deallocate(UpdMask)
*
      Call QExit('Update')
      Return
      End
      Subroutine Update_sl_(kIter,iInt,nFix,nInter,qInt,Shift,
     &                     Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,
     &                     Energy,UpMeth,ed,Line_Search,Step_Trunc,
     &                     nLambda,iRow_c,nsAtom,AtomLbl,nSym,iOper,
     &                     mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                     rLambda,ipCx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                     iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                     nWndw,Mode,ipMF,
     &                     iOptH,HUpMet,mIter,GNrm_Threshold,IRC,
     &                     dMass,HrmFrq_Show,CnstWght,Curvilinear,
     &                     Degen,Kriging_Hessian,qBeta,Restriction,
     &                     iOpt_RS)
************************************************************************
*     Object: to update coordinates                                    *
*                                                                      *
*    Input:                                                            *
*      kIter          : iteration counter                              *
*      iInt           : number of internal coordinates to vary         *
*      nFix           : number of frozen internal coordinates          *
*      nInter         : total number of internal coordinates           *
*      qInt(*,kIter)  : the internal coordinates                       *
*      Shift(*,kIter) : the shift of the internal coordinates          *
*      Grad(*,kIter)  : the gradient in the internal coordinates       *
*      iOptC          : option flag for update methods                 *
*      Beta           : damping factor step length                     *
*      Beta_Disp      : damping factor variance                        *
*      Lbl            : character labels for internal coordinates      *
*      nLbl           : length of Lbl                                  *
*      GNrm           : the norm of the gradient in each iteration     *
*      Energy         : the energy of each iteration                   *
*      Line_Search    : logical flag for line search                   *
*      nLambda        : number of contraints                           *
*      iRow_c         : number of lines on the UDC file                *
*      nsAtom         : number of symmetry unique atoms                *
*      AtomLbl        : character string with atom labels              *
*      nSym           : number of irreps                               *
*      iOper          : integer representations of symmetry operators  *
*      mxdc           : max number of nsAtom                           *
*      jStab          : integer list of stabilizers                    *
*      nStab          : number of stabilizers                          *
*      BMx            : the so-called Wilson B matrix                  *
*      Smmtrc         : logical flag for symmetry properties           *
*      nDimBC         : dimension of redundant coordinates(?)          *
*      rLambda        : vector for Lagrange multipliers                *
*      ipCx           : pointer to cartesian coordinates               *
*      iNeg           : Hessian index                                  *
*      Labels         : character string of primitive int. coord.      *
*      nLabels        : length of Labels                               *
*      CnstWght       : constraints weight                             *
*                                                                      *
*    OutPut:                                                           *
*      qInt(*,kIter+1): the internal coordinates to be used in the     *
*                       next iteration                                 *
*      UpMeth         : character label with update method abrivation  *
*      ed             : estimated energy change to the next point      *
*      Step_Trunc     : character label to denote truncation type      *
*      GrdMax         : largest gradient component                     *
*      StpMax         : largest step component                         *
*      GrdLbl         : character label of component with GrdMax       *
*      StpLbl         : character label of component with StpLbl       *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh                                             *
*             2000                                                     *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      External Restriction
      Real*8 Restriction
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8 qInt(nInter,kIter+1), Shift(nInter,kIter),
     &       Grad(nInter,kIter), GNrm(kIter), Energy(kIter),
     &       dMass(nsAtom), BMx(3*nsAtom,3*nsAtom),
     &       rLambda(nLambda,kIter+1), Degen(3*nsAtom)
      Integer iOper(0:nSym-1), jStab(0:7,nsAtom), nStab(nsAtom),
     &        iNeg(2)
      Logical Line_Search, Smmtrc(3*nsAtom),
     &        FindTS, TSC, HrmFrq_Show,Found,
     &        Curvilinear, Kriging_Hessian
      Character Lbl(nLbl)*8, GrdLbl*8, StpLbl*8, Step_Trunc,
     &          Labels(nLabels)*8, AtomLbl(nsAtom)*(LENIN), UpMeth*6,
     &          HUpMet*6, File1*8, File2*8
      Real*8, Allocatable:: Hessian(:,:)
      iRout=153
      iPrint=nPrint(iRout)
      Lu=6
      If (iPrint.ge.99) Then
         Write (Lu,*)'Update_:iOpt_RS,Beta,Beta_Disp=',
     &                        iOpt_RS,Beta,Beta_Dispa
         Call RecPrt('Update_: qInt',' ',qInt,nInter,kIter)
         Call RecPrt('Update_: Shift',' ',Shift,nInter,kIter-1)
         Call RecPrt('Update_: GNrm',' ',GNrm,kIter,1)
      End If
*
      GrdMax=Zero
      StpMax=Zero
      qBeta=Beta
*
      mInter=nInter
      nA = (Max(mInter,kIter)+1)**2
      nAF= (Max(mInter,kIter)+1)**2
      Call GetMem(' A ','Allo','Real',ipA,nA)
      Call GetMem(' dg','Allo','Real',ipdg,mInter)
      Call GetMem(' Pivot','Allo','Inte',iPvt,kIter+1)
      Call GetMem(' Index','Allo','Inte',iP,kIter)
      Call GetMem('ErrVec','Allo','Real',ipErr,mInter*(kIter+1))
      Call GetMem('EMtrx ','Allo','Real',ipEMx,(kIter+1)**2)
      Call GetMem('RHS   ','Allo','Real',ipRHS,kIter+1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Generate the Hessian in internal coordinates and then update it
*     according to some Hessian update method (BFGS, MSP, etc.)
*
*
      Call mma_Allocate(Hessian,nInter,nInter,Label='Hessian')
      If (Kriging_Hessian) Then
*
*        Temporary code until we have the 2nd derivatives from the
*        kriging code.
*
         Call DCopy_(nInter**2,[Zero],0,Hessian,1)
         Call DCopy_(nInter,[1.0D-2],0,Hessian,nInter+1)
         Call Hessian_Kriging(qInt(1,kIter),Hessian,nInter)
         iNeg(1)=0
         iNeg(2)=0
         HUpMet=' GPR'
      Else
         Call Mk_Hss_Q()
         Call Get_dArray('Hss_Q',Hessian,nInter**2)
*
*        Perform the Hessian update
*
         If (iPrint.ge.6) Then
            Write (Lu,*)
            Write (Lu,*)
            Write (Lu,*) ' *** Updating the molecular Hessian ***'
            Write (Lu,*)
         End If
         iRout=154
         jPrint=nPrint(iRout)
         Call Update_H(nWndw,Hessian,nInter,
     &                 mIter,iOptC,Mode,ipMF,
     &                 Shift(1,kIter-mIter+1),Grad(1,kIter-mIter+1),
     &                 iNeg,iOptH,HUpMet,nRowH,jPrint,GNrm(kIter),
     &                 GNrm_Threshold,nsAtom,IRC,.True.)
      End If
*
*     Save the number of internal coordinates on the runfile.
*
      Call Put_iScalar('No of Internal coordinates',mInter)
*
*     Save the force constant matrix on the runfile.
*
      Call Put_dArray('Hess',Hessian,mInter**2)
*
*     Optional frequency analysis
*
      iDo_dDipM=0
      If (HrmFrq_Show) Call GF_on_the_fly(iDo_DipM)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- If this is a constrained optimization with the goal of finding a
*     transition state and we have one negative eigenvalue of the
*     Hessian then shift over to Mode Following RF to find the TS.
*
*     If TSConstraints were given, merge them with global constraints if
*     they exist, unless we are in TS regime.
*     If TS constraints were not given, remove global constraints when in
*     TS regime.
*
      Call qpg_darray('TanVec',Found,nRP)
      If (FindTS) Then
         File1='UDC'
         File2='TSC'
         If (.not.TSC) File2=''
         If (iNeg(1).ge.1 ) Then
            If ((GNrm(kIter).le.GNrm_Threshold).or.Found) Then
*              Change to MFRF optimization.
               Mask=1+2+4+8+16+32+64+256+512+1024+2048+8192
               iOptC=iAnd(iOptC,Mask)
               Call Put_lScalar('TS Search',.True.)
               If (TSC) Then
                  File2=''
               Else
                  File1=''
               End If
            End If
         End If
         Call Merge_Constraints(File1,File2,'UDC',nLambda,iRow_c)
         Call Fix_UDC(iRow_c,nLambda,AtomLbl,nsAtom,nStab,.False.)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If ( nLambda.eq.0) Then
         M=3*nsAtom
         N=nInter
         NRHS=1
         Call Allocate_Work(ip_Tmp,M)
*
*        Iterate to avoid too large displacement in cartesians.
*        If the maximum displacement is more than 2*Beta, reduce the step.
*
*        Initial setup to ensure fCart=1.0 at first iteration.
         rInter=Beta
         fCart=Ten
         rCart=fCart*rInter
         nLoop=0
         Do While (rCart.ge.Two*Beta)
            nLoop=nLoop+1
            If (nLoop.gt.100) Exit
            If (rCart.gt.rInter) Then
              fCart=fCart*rInter/rCart
            Else
              fCart=fCart*0.9D0
            End If
            mInter=nInter+nLambda
            nA = (Max(mInter,kIter)+1)**2
            nAF= (Max(mInter,kIter)+1)**2
*                                                                      *
************************************************************************
*                                                                      *
            gBeta=One
            xBeta=One
            gg_last=Beta
            dxdx_last=Beta
            Sf=Sqrt(Two)
            kStart=Max(1,kIter-4)
            Thr=1.0D-6
            Do iIter = kStart, kIter
*
               If (iIter.ne.kIter) Then
                  dxdx=
     &             Sqrt(DDot_(mInter,Shift(1,iIter),1,Shift(1,iIter),1))
*
                  If (dxdx.lt.0.75D0*dxdx_last.and.
     &                dxdx.lt.(Beta-Thr)) Then
*                    Increase trust radius
C                    xBeta=xBeta*Sf
                     xBeta=Min(Two,xBeta*Sf)
                  Else If (dxdx.gt.1.25D0*dxdx_last.or.
     &                     dxdx.ge.(Beta+Thr)) Then
*                    Reduce trust radius
                     xBeta=Max(One/Five,xBeta/Sf)
                  End If
                  dxdx_last=dxdx
               End If
*
               gg=Sqrt(DDot_(nInter-nLambda,Grad(1,iIter),1,
     &                                     Grad(1,iIter),1))
               If (gg.lt.0.75D0*gg_last.and.gg.lt.(Beta-Thr)) Then
*                 Increase trust radius
C                 gBeta=gBeta*Sf
                  gBeta=Min(Two,gBeta*Sf)
               Else If (gg.gt.1.25D0*gg_last.or.gg.ge.(Beta+Thr)) Then
*                 Reduce trust radius
                  gBeta=Max(One/Five,gBeta/Sf)
               End If
               gg_last=gg
*
            End Do
            tBeta= Max(Beta*Min(xBeta,gBeta),Beta/Ten)
C           Write (*,*) 'tBeta=',tBeta
*                                                                      *
************************************************************************
*                                                                      *
*----------... Compute updated geometry in Internal coordinates
*
*           Select restriction if step or variance.
            If (iOpt_RS.eq.0) Then
               qBeta=fCart*tBeta
               Thr_RS=1.0D-7
            Else
               qBeta=Beta_Disp
               Thr_RS=1.0D-5
            End If
            Call Newq(qInt,mInter,kIter,Shift,Hessian,Grad,
     &                Work(ipErr),Work(ipEMx),Work(ipRHS),iWork(iPvt),
     &                Work(ipdg),Work(ipA),nA,
     &                ed,iOptC,qBeta,nFix,iWork(ip),UpMeth,
     &                Energy,Line_Search,Step_Trunc,
     &                Restriction,Thr_RS)
            Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,mInter,
     &                  Grad(1,kIter),Shift(1,kIter),Lbl)
*
*---------- Set shift vector to zero for frozen internal coordinates.
*
           If (nFix.gt.0)
     &        call dcopy_(nFix,[Zero],0,Shift(iInt+1,kIter),1)
*
*           Rough conversion to Cartesians
*
            Call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Shift(1,kIter),Work(ip_Tmp))
            rInter=Sqrt(dDot_(N,Shift(1,kIter),1,Shift(1,kIter),1))
            rCart=Zero
            Do i=ip_Tmp,ip_Tmp+nsAtom-1
               rCart=Max(rCart,
     &                   Sqrt(dDot_(3,Work(i),nsAtom,Work(i),nsAtom)))
            End Do
         End Do
         Call Free_Work(ip_Tmp)
         Call Put_dScalar('Max error',Zero)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*        Optimization with constraints using Lagrangian technique.
*
*        Allocate memory for the new arrays
*
         Call Getmem('ipR',   'Allo','Real',ipR,nLambda*kIter)
         Call Getmem('ipdRdq','Allo','Real',ipdRdq,nInter*nLambda*kIter)
         Call Getmem('ipT',   'Allo','Real',ipT,nInter*nInter)
         Call FZero(Work(ipT),nInter**2)
         Call Getmem('ipdy',  'Allo','Real',ipdy,nLambda)
         Call GetMem('dx',    'Allo','Real',ipdx,(nInter-nLambda)*kIter)
         Call FZero(Work(ipdx),(nInter-nLambda)*kIter)
         Call GetMem('dEdq_', 'Allo','Real',ipdEdq_,nInter*kIter)
         Call GetMem('du',    'Allo','Real',ipdu,nInter)
         Call GetMem('x','Allo','Real',ipx,(nInter-nLambda)*(kIter+1))
         Call GetMem('dEdx','Allo','Real',ipdEdx,(nInter-nLambda)*kIter)
         Call GetMem('W   ','Allo','Real',ipW   ,(nInter-nLambda)**2)
         Call GetMem('Energy','Allo','Real',ipEnergy,kIter)
*
         call dcopy_(kIter,Energy,1,Work(ipEnergy),1)
*
         If (mInter.gt.nLbl) Then
            Call WarningMessage(2,'Update_Sl: mInter.gt.nLbl')
            Write (6,*) 'mInter=',mInter
            Write (6,*) 'nLbl=',nLbl
            Call Abend()
         End If
*
         nBVec=iRow_c-nLambda-1
         If (nBVec.gt.nLabels) Then
            Call WarningMessage(2,'Update_Sl: nBVec.gt.nLabels')
            Write (6,*) 'nBVec=',nBVec
            Write (6,*) 'nLabels=',nLabels
            Call Abend()
         End If
         Call GetMem('dBMx', 'Allo','Real',ipdBMx,nLambda*(3*nsAtom)**2)
         Call GetMem('BMx',  'Allo','Real',ipBMx,3*nsAtom*nLambda)
*
         Call GetMem('BVec', 'Allo','Real',ipBVec,3*nsAtom*nBVec)
         Call GetMem('Value','Allo','Real',ipValue,nBVec)
         Call GetMem('Value0','Allo','Real',ipValue0,nBVec)
         Call GetMem('cInt', 'Allo','Real',ipcInt,nLambda)
         Call GetMem('cInt0','Allo','Real',ipcInt0,nLambda)
         Call GetMem('Mult', 'Allo','Real',ipMult,nBVec**2)
         Call GetMem('iFlip','Allo','Inte',ip_iFlip,nBVec)
         Call GetMem('dBVec','Allo','Real',ipdBVec,nBVec*(3*nsAtom)**2)
*
*        Compute the constraints
*
         n1=3*nsAtom
         n2=n1**2
         ip_drdq=ipdrdq
         Do lIter = 1, kIter
            ipCoor_l=ipCx + (lIter-1)*n1
            Call DefInt2(Work(ipBVec),Work(ipdBVec),nBVec,Labels,
     &                   Work(ipBMx),nLambda,nsAtom,iRow_c,
     &                   Work(ipValue),Work(ipcInt),Work(ipcInt0),
     &                   Lbl(nInter+1),AtomLbl,Work(ipCoor_l),
     &                   (lIter.eq.kIter),nSym,iOper,jStab,nStab,mxdc,
     &                   Work(ipMult),Smmtrc,nDimBC,Work(ipdBMx),
     &                   Work(ipValue0),lIter,iWork(ip_iFlip),dMass)
*
*           Assemble r
*
            iOff = ipr + (lIter-1)*nLambda
            Do j = 0, nLambda-1
               Ci=Work(ipcInt+j)-Work(ipcInt0+j)
               Work(iOff)=Ci
               iOff = iOff + 1
            End Do
*
*           Assemble dr/dq: Solve  B dr/dq = dr/dx
*
            ip_drdq=ipdrdq+(lIter-1)*nInter*nLambda
            Call FZero(Work(ip_drdq),nInter*nLambda)
*           Call RecPrt('BMx',' ',BMx,3*nsAtom,nInter)
*           Call RecPrt('Work(ipBMx)',' ',Work(ipBMx),3*nsAtom,nLambda)
*
            M=3*nsAtom
            N=nInter
            NRHS=nLambda
*
*           Temporary fix of the dC/dx vector which always
*           is propted up with the full degeneracy factor.
*
            If (.NOT.Curvilinear) Then
               Do iLambda=1,nLambda
                  Do i = 1, 3*nsAtom
                     ij = (iLambda-1)*3*nsAtom + i -1
                     Work(ij+ipBMx)=Work(ij+ipBMx)/Degen(i)
                  End Do
               End Do
            End If
            LudRdX=30
            Call DaName(LudRdX,'dRdX')
            iAd=0
            Call iDaFile(LudRdX,1,[nLambda],1,iAd)
            Call iDaFile(LudRdX,1,[3*nsAtom],1,iAd)
            Call dDaFile(LudRdX,1,Work(ipBMx),nLambda*3*nsAtom,iAd)
            Call DaClos(LudRdX)
            Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Work(ipBMx),Work(ip_drdq))
*           Call RecPrt('drdq',' ',Work(ip_drdq),nInter,nLambda)
*
         End Do     ! lIter
         Call GetMem('dBVec','Free','Real',ipdBVec,nBVec*(3*nsAtom)**2)
         Call GetMem('iFlip','Free','Inte',ip_iFlip,nBVec)
         Call GetMem('Mult', 'Free','Real',ipMult,nBVec**2)
         Call GetMem('cInt0','Free','Real',ipcInt0,nLambda)
         Call GetMem('cInt', 'Free','Real',ipcInt,nLambda)
         Call GetMem('Value0','Free','Real',ipValue0,nBVec)
         Call GetMem('Value','Free','Real',ipValue,nBVec)
         Call GetMem('BVec', 'Free','Real',ipBVec,3*nsAtom*nBVec)
*
         If (iPrint.ge.99) Then
            Call RecPrt('r',' ',Work(ipr),nLambda,kIter)
            Call RecPrt('drdq',' ',Work(ipdrdq),nInter*nLambda,kIter)
            Call RecPrt('dC/dx',' ',Work(ipBMx),n1,nLambda)
            Call RecPrt('dQ/dx',' ',BMx,n1,nInter)
            Do i = 1, nLambda
               jpx=ipdBMx+(i-1)*n2
               Write (6,*) ' iLambda=',i
               Call RecPrt('d2C/dx2',' ',Work(jpx),n1,n1)
            End Do
            If (Curvilinear) Call dBPrint(nInter,nDimBC)
         End If
         Call GetMem('BMx',  'Free','Real',ipBMx,3*nsAtom*nLambda)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Assemble d^2C/dQ^2
*
*        The expression is written as
*
*        dQ/dx * d^2C/dQ^2 * dQ/dx = d^2C/dx^2 - Sum (d^2Q/dx^2 * dC/dQ)
*
*        Solve first for y in
*                dQ/dx y = d^2C/dx^2 - Sum (d^2Q/dx^2 x dC/dQ)
*
#ifdef _NEW_CODE_
*        Compute Sum (d^2Q/dx^2 * dC/dQ)
*
         Call GetMem('d2Qdx2_dCdQ','Allo','Real',ip_QC,
     &                nDimBC**2*nLambda)
         If (Curvilinear)
     &      Call dBMult(Work(ip_drdq),Work(ip_QC),nInter,nDimBC,nLambda)
         Call RecPrt('Work(ip_drdq)',' ',Work(ip_drdq),n1,1)
         Call RecPrt('Work(ip_QC)',' ',Work(ip_QC),nDimBC**2,nLambda)
*
*        Subtract the term d^2C/dx^2
*
         Do k = 1, nLambda
            j = 0
            Do jx = 1, n1
               If (Smmtrc(jx)) Then
                  j = j + 1
*
                  i=0
                  Do ix = 1, n1
                     If (Smmtrc(ix)) Then
                        i = i + 1
                        ijk = (k-1)*n2 + (jx-1)*n1 + ix-1 + ipdBMx
                        ijk1= (k-1)*nDimBC**2 + (j-1)*nDimBC
     &                      + i-1 + ip_QC
                        Work(ijk) = Work(ijk) - Work(ijk1)
                     End If
                  End Do
               End If
            End Do
         End Do
         Call Free_Work(ip_QC)
#endif
*
         Call GetMem('Scr2','Allo','Real',ipScr2,nInter*n1*nLambda)
         Call GetMem('Scr1','Allo','Real',ipScr1,nInter*n1*nLambda)
*        Call RecPrt('d^2C/dx^2',' ',Work(ipdBMx),n2,nLambda)
*
*        Temporary fix of the d^2C/dx^2 vector which always
*        is propted up with the full degeneracy factor.
*
         If (.NOT.Curvilinear) Then
            Do k = 1, nLambda
               Do j = 1, n1
                  Do i = 1, n1
                     ijk = (k-1)*n2 + (j-1)*n1 + i-1 + ipdBMx
                     Work(ijk)=Work(ijk)/(Degen(i)*Degen(j))
                  End Do
               End Do
            End Do
         End If
*
         M=3*nsAtom
         N=nInter
         NRHS=n1*nLambda
*        Call RecPrt('dQ/dx',' ',BMx,n1,nInter)
         Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,
     &                  Work(ipdBMx),Work(ipScr1))
         Call GetMem('dBMx', 'Free','Real',ipdBMx,nLambda*(3*nsAtom)**2)
         Call TRNSPS(nInter,n1*nLambda,Work(ipScr1),Work(ipScr2))
*        Call RecPrt('Scr1',' ',Work(ipScr2),n1*nLambda,nInter)
*
*        Followed by solve for x in B x = y^T for which x = d^2r/dq^2
*
         M=3*nsAtom
         N=nInter
         NRHS=nLambda*nInter
         Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,
     &                  Work(ipScr2),Work(ipScr1))
         Call Free_Work(ipScr2)
         Call GetMem('d2L','Allo','Real',ipd2L,nInter*nInter*nLambda)
         Call TRNSPS(nInter*nLambda,nInter,Work(ipScr1),Work(ipd2L))
         Call Free_Work(ipScr1)
*        Call RecPrt('ipd2L',' ',Work(ipd2L),nInter**2,nLambda)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Compute updated geometry in Internal coordinates
*
         Mode_=0
         M=3*nsAtom
         N=nInter
         NRHS=1
         Call Allocate_Work(ip_Tmp,M)
*
*        Iterate to avoid too large displacement in cartesians.
*        If the maximum displacement is more than 2*Beta, reduce the step.
*
*        Initial setup to ensure fCart=1.0 at first iteration.
         rInter=Beta
         fCart=Ten
         rCart=fCart*rInter
         nLoop=0
         Do While (rCart.ge.Two*Beta)
            nLoop=nLoop+1
            If (nLoop.gt.100) Exit
            If (rCart.gt.rInter) Then
              fCart=fCart*rInter/rCart
            Else
              fCart=fCart*0.9D0
            End If
            qBeta=fCart*Beta
            Call Con_Opt(Work(ipr),Work(ipdrdq),Work(ipT),Grad,
     &                rLambda,qInt,Shift,Work(ipdy),Work(ipdx),
     &                Work(ipdEdq_),Work(ipdu),Work(ipx),Work(ipdEdx),
     &                Work(ipW),GNrm(kIter),
     &                nWndw,Hessian,nInter,kIter,
     &                iOptC,Mode_,ipMF,iOptH,HUpMet,jPrint,
     &                Work(ipEnergy),nLambda,mIter,nRowH,
     &                Work(ipErr),Work(ipEMx),Work(ipRHS),iWork(iPvt),
     &                Work(ipdg),Work(ipA),nA,ed,qBeta,Beta_Disp,nFix,
     &                iWork(iP),UpMeth,Line_Search,Step_Trunc,Lbl,
     &                GrdLbl,StpLbl,GrdMax,StpMax,Work(ipd2L),nsAtom,
     &                IRC,CnstWght,Restriction,iOpt_RS,Thr_RS)
*
*           Rough conversion to Cartesians
*
            Call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Shift(1,kIter),Work(ip_Tmp))
            rInter=Sqrt(dDot_(N,Shift(1,kIter),1,Shift(1,kIter),1))
            rCart=Zero
            Do i=ip_Tmp,ip_Tmp+nsAtom-1
               rCart=Max(rCart,
     &                   Sqrt(dDot_(3,Work(i),nsAtom,Work(i),nsAtom)))
            End Do
         End Do
         Call Free_Work(ip_Tmp)
*
         If (iPrint.ge.99) Then
            Write (Lu,*)
            Write (Lu,*) '********************************************'
            Write (Lu,*) '* Lagrange multipliers for the constraints *'
            Write (Lu,*) '********************************************'
            Write (Lu,'(1X,A,2X,ES13.6)')
     &          (Lbl(nInter+iInt),-One*rLambda(iInt,mIter),
     &           iInt=1,nLambda)
            Write (Lu,*)
         End If
*
         Call Free_Work(ipd2L   )
         Call Free_Work(ipEnergy)
         Call Free_Work(ipW     )
         Call Free_Work(ipdEdx  )
         Call Free_Work(ipx     )
         Call Free_Work(ipdu    )
         Call Free_Work(ipdEdq_ )
         Call Free_Work(ipdx    )
         Call Free_Work(ipdy    )
         Call Free_Work(ipT     )
         Call Free_Work(ipdrdq  )
         Call Free_Work(ipr     )
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Deallocate(Hessian)
      Call GetMem('RHS   ','Free','Real',ipRHS,kIter+1)
      Call GetMem('EMtrx ','Free','Real',ipEMx,(kIter+1)**2)
      Call GetMem('ErrVec','Free','Real',ipErr,mInter*(kIter+1))
      Call GetMem(' Index','Free','Inte',iP,kIter)
      Call GetMem(' Pivot','Free','Inte',iPvt,kIter+1)
      Call GetMem(' dg','Free','Real',ipdg,mInter)
      Call GetMem(' A ','Free','Real',ipA,nA)
*
      Return
      End
