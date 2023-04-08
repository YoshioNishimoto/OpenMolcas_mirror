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
        subroutine o3v3ctl (wrk,wrksize,NvGrp,LunAux)
!
!       docasny drajver o3v3 procesov
!
        use Para_Info, only: nProcs
        implicit none
#include "chcc1.fh"
#include "parcc.fh"
#include "o3v3.fh"
#include "wrk.fh"
!
        integer NvGrp,LunAux,maxdim
        integer proc,beGrp,aGrp
        integer nJobs,addJobs,actJobs,actNode
!
!
!1      Def parameters for o3v3 processes
!
        call DefParo3v3 (NvGrp,maxdim)
!
!
!2        Distribute work among nodes (def BetaID, BeAID)
!
!2.1        vanish BetaID, BeAID
        do proc=0,nProcs-1
        do beGrp=1,NvGrp
          BetaID(proc,beGrp)=0
          do aGrp=1,NvGrp
            BeAID(proc,beGrp,aGrp)=0
          end do
        end do
        end do
!
!
        if (nProcs.eq.1) then
!2.2.1        single node
!
           do beGrp=1,NvGrp
             BetaID(0,beGrp)=1
            do aGrp=1,NvGrp
              BeAID(0,beGrp,aGrp)=1
            end do
           end do
!
        else
!2.2.2        multi node
!        (N.B. trochu odflaknute, dalo by sa este zohladnit ze tie nody,
!         ktore maju o jeden job navyse nebudu tie, kde je BetaID=1,
!        resp. nebudu tie, kde sa realizuje X0.1 prispevok a pod)
!
          nJobs=int((NvGrp*NvGrp)/nProcs)
          addJobs=mod((NvGrp*NvGrp),nProcs)
!
          actNode=0
          actJobs=nJobs
!
!
          do beGrp=1,NvGrp
          BetaID(actNode,beGrp)=1
          if (printkey.ge.10) then
          write (6,*) 'BetaID',actnode,beGrp
          end if
          do aGrp=1,NvGrp
            BeAID(actNode,beGrp,aGrp)=1
            if (printkey.ge.10) then
            write (6,*) 'BeAID',actnode,beGrp,aGrp
            end if
            actJobs=actJobs-1
            if (actJobs.eq.-1) then
              actNode=actNode+1
              actJobs=nJobs
            else if (actJobs.eq.0) then
              if (addJobs.gt.0) then
                addJobs=addJobs-1
              else
                actNode=actNode+1
                actJobs=nJobs
              end if
            end if
          end do
          end do
!
        end if
!@@
        if (printkey.ge.10) then
        do proc=0,nProcs-1
        do beGrp=1,NvGrp
          write (6,99) proc,beGrp,(BeAID(proc,beGrp,aGrp),aGrp=1,NvGrp)
        end do
        end do
        end if
99        format (1x,i3,1x,i2,5x,24(i1,1x))
!@@
!
!
!3      A ideme na to
!
        call o3v3jk (wrk(1),wrksize,NvGrp,maxdim,LunAux)
        if (printkey.gt.1) then
        write (6,*) ' o3v3jk done'
        end if
!
        call o3v3chol (wrk(1),wrksize,NvGrp,maxdim,LunAux)
        if (printkey.gt.1) then
        write (6,*) ' o3v3chol done'
        end if
!
        call o3v3t2 (wrk(1),wrksize,NvGrp,maxdim,LunAux)
        if (printkey.gt.1) then
        write (6,*) ' o3v3t2 done'
        end if
!
        return
        end
