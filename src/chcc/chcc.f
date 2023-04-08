!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
        subroutine chcc (ireturn)
!
!       docasny drajver reorder procesu
!
        use Para_Info, only: MyRank, nProcs
#ifdef _MOLCAS_MPP_
        use Para_Info, only: Is_Real_Par
#endif
        implicit none
#include "chcc1.fh"
#include "parcc.fh"
#include "chcc_files.fh"
#include "chcc_reord.fh"
!#include "SysDef.fh"
#include "WrkSpc.fh"
#include "chcc_casy.fh"
        integer ireturn
!
        integer NvGrp,NvSGrp,NchBlk
        integer LunAux
        integer wrksize
        integer maxspace,iOff
!
        integer maxdim
        integer NChHere
!        jalove
        integer Jal1,Jal2
!
        real*8 e1new,e2new,e1old,e2old,e2os,escf
        integer idum,iter
!mp
!mp
!
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
         write (6,*) ' Parallel run'
        else
         write (6,*) ' Serial run'
        end if
!mp!         write (6,*) ' MyRank, Nprocs',MyRank, Nprocs
#else
         write (6,*) ' Serial run'
!mp!         write (6,*) ' MyRank, Nprocs',MyRank, Nprocs
#endif
!mp
!        vynuluj hodiny
        Call CWTime(TCpu,TWall)
        TCpu0=TCpu
        TWall0=TWall
        TCpu_l=TCpu
        TCpu_l=TCpu
        TWall_l=TWall
        TWall_l=TWall
!mp
!
!@@     real*8, allocatable :: wrk(:)
!        real*8 wrk(1:40000000)
!        real*8 wrk(1:26593281)
!        wrksize=26593281
!        wrksize=40000000
!
!mp! ##########################################################
!0      info o cholesky vektoroch

        call frankie_drv_fake (NChHere)
        write (6,'(A,i9,A,i4)') ' Number of Cholesky vectors ',         &
     & NChHere,' on node ',myRank
!
#ifdef _MOLCAS_MPP_
!
        do jal1=0,Nprocs-1
          NChLoc(jal1)=0
        end do
!
        NChLoc(MyRank)=NChHere
        call gaigop (NChLoc(0),NProcs,'+')
!
        jal2=0
        do jal1=0,NProcs-1
          jal2=jal2+NChLoc(jal1)
!mp!          write (6,*) ' NChLoc',jal1,NChLoc(jal1)
        end do

        nc = jal2
#else
        nc = NChHere
#endif
!mp! ##########################################################
!
!       Get the maximum available memory
!
        Call GetMem('CCSD','Max','Real',idum,maxspace)
        maxspace=maxspace-8
        write (6,'(A,i13,A,f9.1,A,f5.1,A)') ' Max Size              : ',&
     & maxspace,' in r*8 Words,',                                       &
     & 1.0d0*maxspace*8/(1024*1024),' Mb,',                             &
     & 1.0d0*maxspace*8/(1024*1024*1024),' Gb'
!
!1      Nacitanie vstupu (Docasne) + time Delay
!
        call IniReord (NvGrp,NvSGrp,NchBlk,LunAux,wrksize)
!mp!        call TimeDelay (LunAux) ! temporarily disabled (MP)
!
!       Decide on automatic segmentation generation
        if (NvGrp.eq.0) then
           call autoSegmentation(Nprocs, maxspace,                      &
     & Jal1, Jal2,                                                      &
     & NvGrp, NvSGrp, NchBlk, wrksize, maxdim)
        else
           call checkMem(NvGrp, NvSGrp, NchBlk,                         &
     & Jal1, Jal2, wrksize, maxdim)

           if (wrksize.gt.maxspace) then
              write (6,*) ' Not Enough Memory! '//                      &
     & 'Increase large and/or small segmentation',                      &
     & (1.0d0*wrksize)/(1.0d0*maxspace)
             call abend
           end if
        endif
!
!         write Large segmentation to RunFile
        Call Put_iScalar('CHCCLarge',NvGrp)
!
! ---------------------------------------------------------
!
!        Transformation of Local L(ao) -> L(mo)
!        CD1tmp file is created
!
        call frankie_drv (NChHere)
        if (printkey.ge.10) then
        write (6,*) ' After Frenkie',myRank,NChHere
        end if
!
! ---------------------------------------------------------
!
        Call GetMem('CCSD','Allo','Real',iOff,wrksize)
        write (6,'(A,i13,A,f9.1,A,f5.1,A)') ' Real Allocated Memory : ',&
     & wrksize,' in r*8 Words,',                                        &
     & 1.0d0*wrksize*8/(1024*1024),' Mb,',                              &
     & 1.0d0*wrksize*8/(1024*1024*1024),' Gb'
        call mv0zero (wrksize,wrksize,Work(iOff))
!
!3      Priprava integralov
!
        call Reord_chcc (Work(iOff),wrksize,                            &
     &              NvGrp,NvSGrp,NchBlk,LunAux)
        if (generkey.eq.1) then
!mp!          if (printkey.ge.10) then ! uvidime ...
          write (6,*) ' Generation of integrals (Reord_chcc) done'
          write (6,*)
!mp!          end if
        else
          write (6,*) ' Generation of integrals (Reord_chcc) skipped,'//&
     &          ' only basic'
          write (6,*)
        end if
!mp
        Call CWTime(TCpu,TWall)
!
        if (printkey.gt.1) then
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',                  &
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',                  &
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',                      &
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',                      &
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',                      &
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
!mp
!
!4        Restart, ak treba
!
        if (restkey.eq.1) then
!        read T1o,niter,E1old,E2old (T2 are in T2files)
          call GetRest (Work(iOff),wrksize,LunAux,iter,E1old,E2old)
          goto 1
        else
!Bug    v originale chyba
!        set T1o=0
           call VanishT1 (Work(iOff),wrksize)
        end if
!
!@@
!        write (6,*) ' Pred Chck'
!        Call MakeChckData (Work(iOff),wrksize,LunAux)
!        Call SaveChckData (LunAux)
!        Call GetChckData (LunAux)
!        write (6,*) ' Chck  done'
!@@
!
!
!5        iteracny cyklus
!
        e1old=0.0d0
        e2old=0.0d0
        iter=1
!
!mp
!
        write (6,*)
        write (6,*) '------------------------'
        write (6,*) 'Starting CCSD iterations'
        write (6,*) '------------------------'
        write (6,*)

        write (6,*) '                  CCSD Energy  ',                  &
     & '    Difference'
        write (6,*)
!mp!
        call xflush(6)
!mp!
!
        Call CWTime(TCpu,TWall)
!
        if (printkey.gt.1) then
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',                  &
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',                  &
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',                      &
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',                      &
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',                      &
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
!mp

1        continue
!
           call o3v3ctl (Work(iOff),wrksize,NvGrp,LunAux)
        if (printkey.gt.1) then
         write (6,*) ' o3v3 done'
        end if
!mp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',                  &
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',                  &
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',                      &
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',                      &
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',                      &
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
!mp
        call o2v4ctl (Work(iOff),wrksize,NvGrp,NvSGrp,LunAux)
        if (printkey.gt.1) then
        write (6,*) ' o2v4 done'
        end if
!mp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',                  &
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',                  &
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',                      &
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',                      &
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',                      &
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
!mp
           call summary (Work(iOff),wrksize,NvGrp,LunAux,maxdim,        &
     &                e1new,e2new,e2os)
        if (printkey.gt.1) then
        write (6,*) ' summary done'
        end if
!mp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',                  &
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',                  &
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',                      &
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',                      &
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',                      &
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
!mp
        call SaveRest (Work(iOff),wrksize,LunAux,(iter+1),E1new,E2new)
!
!mp!        write (6,91) ' Iteration :',iter,e1new,e2new,
!mp!     c                             e1old+e2old-e1new-e2new
!mp!91        format (a12,1x,i3,1x,3(f15.12,1x))

        if (iter.eq.1) then
           write (6,91) ' Iteration :',iter,      e2new
        else
           write (6,93) ' Iteration :',iter,      e2new,                &
     &                             e1old+e2old-e1new-e2new
        end if

        call xflush(6)

91        format (a12,1x,i3,1x,f15.12,1x)
93        format (a12,1x,i3,1x,2(f15.12,1x))
!
        if ((abs(e1old+e2old-e1new-e2new).gt.conv).and.                 &
     &      (iter.lt.maxiter)) then
          e1old=e1new
          e2old=e2new
          iter=iter+1
          goto 1
        end if
!
        write (6,*)
        write (6,*)  ' Final CCSD energy decomposition'
        write (6,92) ' E1 CCSD energy :',e1new
        write (6,92) ' E2 CCSD energy :',e2new
        write (6,92) ' E2 CCSD ss     :',e2new-e2os
        write (6,92) ' E2 CCSD os     :',e2os
        write (6,*)
92        format (a17,1x,f15.12)
!mp! for MOLCAS verify
        Call Get_dScalar('SCF energy',escf)
        Call Add_Info('CHCCene', [e2new], 1, 6)
        Call Add_Info('E_CHCC', [e2new], 1, 6)
        Call Add_Info('E_HYPE', [e2new+escf], 1, 6)
!    for NUMERICAL_GRADIENTS
        Call Put_cArray('Relax Method','CHCC    ',8)
        Call Store_Energies(1,e2new+escf,1)
!mp!
!
!@@        deallocate (wrk)
        Call GetMem('CCSD','Free','Real',iOff,wrksize)
!
        ireturn=0
!
        return
        end
