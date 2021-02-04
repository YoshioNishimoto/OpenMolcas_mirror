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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      subroutine wgtinit(H)
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb
      implicit real(8) (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      real(8) H(Nstate,Nstate)


      if (IPRGLB.GE.DEBUG) then
        write(6,*)' Entered wgtinit.'
      end if

* Initialize array of weights with all zeros
      call dcopy_(Nstate**2,[0.0D0],0,WORK(LDWGT),1)

* Main loop over all states to compute the weights
      do I=1,Nstate

* If it is an XDW-CASPT2 calculation, the weights are computed
        if (IFDW.and.zeta.ge.0.0d0) then
          Ebeta = H(I,I)
* Compute normalization factor Wtot, i.e. the sum of all weights
          do J=1,Nstate
            Ealpha = H(J,J)
            Wtot = 0.0D0
            ! pref = 0.001D0
            do K=1,Nstate
              Egamma = H(K,K)
              Dag = abs(Ealpha - Egamma) + 1.0d-8
              if (abs(H(J,K)).le.1.0d-9) then
                Hag = 0.0d0
              else
                Hag = sqrt(abs(H(J,K)))
              end if
* Compute interaction strength parameter $\xi_{\alpha\gamma}$
*              if (K.ne.J) then
                xi_ag = abs(Dag/Hag)
*              else
* Set explicitly to zero when the numerator cancels out to avoid
* numerical issues
*                xi_ag = 0.0d0
*              end if
              Wtot = Wtot + exp(-zeta*xi_ag)
              write(6,*)'exp(-zeta*xi_ag) = ',exp(-zeta*xi_ag)
            end do
            IJ = (I-1) + Nstate*(J-1)
* Compute weight according to XDW prescription
            Dab = abs(Ealpha - Ebeta) + 1.0d-8
            if (abs(H(J,I)).le.1.0d-9) then
              Hab = 0.0d0
            else
              Hab = sqrt(abs(H(J,I)))
            end if
* Compute interaction strength parameter $\xi_{\alpha\beta}$
*            if (I.ne.J) then
              xi_ab = abs(Dab/Hab)
*            else
* Set explicitly to zero when the numerator cancels out to avoid
* numerical issues
*              xi_ab = 0.0d0
*            end if
            WORK(LDWGT+IJ) = exp(-zeta*xi_ab)/Wtot
            write(6,'(A,I1,A,I1,A)')'(',J,',',I,')'
            write(6,'(A,F18.12)')' Dab = ',Dab
            write(6,'(A,F18.12)')' Hab = ',Hab
            write(6,'(A,F18.12)')'xiab = ',xi_ab
          end do
          write(6,*)

* If it is an XMS-CASPT2 calculation, all the weights are equal,
* i.e. they all are 1/Nstate
        else if (IFXMS.and.(.not.IFDW)) then
          call dcopy_(Nstate**2,[1.0D0/Nstate],0,WORK(LDWGT),1)

* If it is a normal MS-CASPT2 or a (X)DW-CASPT2 with zeta->infinity
* the weight vectors are the standard unit vectors e_1, e_2, ...
        else
          WORK(LDWGT + (Nstate*(I-1)) + (I-1)) = 1.0d0
        end if

* End of loop over states
      end do

* In case it is a XDW calculation, print out the weights
      if (IFDW.and.(IPRGLB.ge.VERBOSE)) then
        if (IFEFOCK) then
          write(6,*)' Weights calculated with <I|H0|I>:'
        else
          write(6,*)' Weights calculated with <I|H|I>:'
        end if
        call prettyprint(WORK(LDWGT),Nstate,Nstate)
      end if


      return
      end
