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
        subroutine Xo2v4ctl (NvGrp,NvSGrp,LunAux)
!
!!      drajver procesu na testovanie ktory W3/W4 file
!        treba na ktorom node
!        N.B. upraveny drajver o2v4 procesu
!
        use Para_Info, only: nProcs
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "parcc.fh"
!
        integer NvGrp,NvSGrp,LunAux
!
!        help variables
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
!##
        integer aGrp,bGrp,proc,i,j
        integer nJobs,addJobs,actJobs
!
!
!1      Inicializacia premennych (Predbezna)
!
        NaGrp=NvGrp
        NbeGrp=NvGrp
        NaSGrp=NvSGrp
        nbeSGrp=NvSGrp
!
!
!2      define all groups and ssungroup parameters
!
        call DefParo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,                   &
     &                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
!
!3        distribute work among nodes (def ABID)
!
        do proc=0,(nProcs-1)
        do aGrp=1,NaGrp
        do bGrp=1,NaGrp
          ABID(proc,aGrp,bGrp)=0
        end do
        end do
        end do
!
!
        if ((nProcs.eq.1).or.(NaGrp.eq.1)) then
!3.1        nP = <1>
!        all jobs to one node
!
          do aGrp=1,NaGrp
          do bGrp=1,NaGrp
            ABID(0,aGrp,bGrp)=1
          end do
          end do
!
!
        else
!3.2        nP = <2, inf)
!        N.B. cele zatial trochu odflaknute, lebo sa neberie do uvahy,
!        ze pre half joby je lepsie ked nadvazuju na plne joby s
!        rovnakymi indexami
!
!3.2.1        Full Jobs
!
          i=(NaGrp*(NaGrp-1))/2
          nJobs=int(i/nProcs)
          addJobs=mod(i,nProcs)
!
!          first nodes: 0-addJobs
!             will have int(N'(N'-1)/2 /nProc)+1 Full Jobs
!          rest nodes: addJobs-nProcs-1
!               will have int(N'(N'-1)/2 /nProc)
!
          proc=0
          actJobs=nJobs
          do aGrp=2,NaGrp
          do bGrp=1,aGrp-1
            ABID(proc,aGrp,bGrp)=1
            actJobs=actJobs-1
            if (actJobs.eq.-1) then
              proc=proc+1
              actJobs=nJobs
            else if (actJobs.eq.0) then
              if (addJobs.gt.0) then
                addJobs=addJobs-1
              else
                proc=proc+1
                actJobs=nJobs
              end if
            end if
          end do
          end do
!
!
!3.2.2        Half Jobs
!
!        distribution of Half Jobs
!          - first distribution - addjobs-nProcs-1
!          - second distribution - addjobs-nProcs-1
!          - other distributions - 0 - nProcs-1
!
          addJobs=mod(i,nProcs)
          proc=addJobs
          j=0
!
          do aGrp=1,NaGrp
            ABID(proc,aGrp,aGrp)=1
            proc=proc+1
            if (proc.eq.nProcs) then
              if (j.eq.0) then
!              first distribution - addjobs-nProcs-1
                proc=addJobs
                j=1
              else if (j.eq.1) then
!              second distribution - addjobs-nProcs-1
                proc=addJobs
                j=2
              else
!              other distributions - 0 - nProcs-1
                proc=0
              end if
            end if
          end do
!
        end if
!
!
!        Printing ABID
        do proc=0,nProcs-1
        if (printkey.ge.10) then
        write (6,*)   ' For myRank = ',proc
        end if
          do aGrp=1,NaGrp
          do bGrp=1,aGrp
          if (ABID(proc,aGrp,bGrp).eq.1) then
          if (printkey.ge.10) then
          write (6,*) '    aGrp,bGrp ',aGrp,bGrp
          end if
          end if
          end do
          end do
        end do
!
!
!
!4      A ideme na to
!
        call InsReqo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,                   &
     &             mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
!
!
        return
! Avoid unused argument warnings
        if (.false.) call Unused_integer(LunAux)
        end
