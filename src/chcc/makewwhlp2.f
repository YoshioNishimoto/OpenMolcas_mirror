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
        subroutine MakeWwHlp2 (Ww,W1,                                   &
     &                       dima,dimb,dimab,dimbe,dimga,key)
!
!       this routine do:
!       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
!       for the case a"=b" , be".ne.ga"
!       N.B. algoritmus nieje prilis vymakany
!
!       parameter description
!       Ww     - array for Ww+(-) (O)
!       W1     - array for W1(a",be",b",ga") (I)
!       dimx   - dimension of a",b",ab",be",ga" (I)
!       key    - 1 - calc Ww+, 2 - calc Ww- (I)
!
        implicit none
        integer dima,dimb,dimab,dimbe,dimga,key
        real*8 Ww(1:dimab,1:dimbe,1:dimga)
        real*8 W1(1:dima,1:dimbe,1:dimb,1:dimga)
!
!       help variables
        integer a,b,ab,be,ga
!
!
        if (key.eq.1) then
          do ga=1,dimga
          do be=1,dimbe
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,be,ga)=W1(a,be,b,ga)+W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        else
          do ga=1,dimga
          do be=1,dimbe
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,be,ga)=W1(a,be,b,ga)-W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        end if
!
!        Cely clen ma Faktor 2, tu teda nevydelim 2
!        call mv0sv (dimab*dimbe*dimga,dimab*dimbe*dimga,
!    c              Ww(1,1,1),0.5d0)
!
        return
        end
