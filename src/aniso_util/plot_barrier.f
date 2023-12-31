************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine plot_barrier(nBlock,nMult,nDIM,E,M)

      IMPLICIT NONE

      Integer, parameter        :: wp=kind(0.d0)
      INTEGER, INTENT(in)   :: nBlock, nMult
      INTEGER, INTENT(in)   :: nDIM(nMult)
      COMPLEX (wp), INTENT(in) :: M(3,nMult,10,nMult,10)
      REAL (wp), INTENT(in) :: E(nBlock)
      ! local variables
      REAL (wp)         :: xstart, xend, dlt, xmin, xmax, emin, emax,
     &                     fact, F, ystart, yend, X, Y, Z, RAVE,
     &                     RAVEMIN, RAVEMAX
      REAL (wp)         :: gnuplot_version
      INTEGER           :: file_number, i, l, j, i1, j1, k
      INTEGER           :: LuPlt, LuData, file_size, StdOut
      LOGICAL           :: file_exist, is_file_open, execute_gnuplot_cmd
      CHARACTER(LEN=100):: line1, line2, fmtx, cdummy
      CHARACTER(LEN=100):: datafile_e, datafile_m,
     &                     plotfile, imagefile, epsfile
      LOGICAL           :: dbg
      Integer, external :: AixRm
      INTEGER           :: iErr, Length
      CHARACTER(LEN=1023) :: realname_plt,
     &                       realname_e_dat, realname_m_dat,
     &                       realname_png, realname_eps, gnuplot_CMD
      INTEGER, EXTERNAL :: IsFreeUnit

      dbg=.false.
      StdOut = 6
      iErr=0

      xmin=        MINVAL(DBLE( M(3,1:nMult,:,1:nMult,:) )) -
     &     0.10_wp*MAXVAL(DBLE( M(3,1:nMult,:,1:nMult,:) ))
      xmax=        MAXVAL(DBLE( M(3,1:nMult,:,1:nMult,:) )) +
     &     0.10_wp*MAXVAL(DBLE( M(3,1:nMult,:,1:nMult,:) ))
      dlt=0.01_wp*(xmax-xmin)

      F=MAXVAL(E(1:nBlock))
      fact=1.0_wp
      IF (F.lt.10.0_wp) THEN
        fact=1.0_wp
      !ELSE IF ( (F.ge.10.0_wp) .AND. (F.lt.100.0_wp) ) THEN
      !  fact=10.0_wp
      !ELSE IF ( (F.ge.100.0_wp) .AND. (F.lt.1000.0_wp) ) THEN
      !  fact=100.0_wp
      !ELSE IF ( (F.ge.1000.0_wp) .AND. (F.lt.10000.0_wp) ) THEN
      !  fact=1000.0_wp
      ELSE IF ( (F.ge.10000.0_wp) .AND. (F.lt.100000.0_wp) ) THEN
        fact=10000.0_wp
      ELSE IF ( (F.ge.100000.0_wp) .AND. (F.lt.1000000.0_wp) ) THEN
        fact=100000.0_wp
      END IF
      Emin=MINVAL(E(1:nBlock))/fact - 0.05_wp*MAXVAL(E(1:nBlock))/fact
      Emax=MAXVAL(E(1:nBlock))/fact + 0.05_wp*MAXVAL(E(1:nBlock))/fact


      RAVEMIN= 9999999999999999.0_wp
      RAVEMAX=-9999999999999999.0_wp

      l=0
      DO i=1,nMult
        DO i1=1,ndim(i)
          l=l+1
          k=0
          DO j=1,nMult
            DO j1=1,ndim(j)
              k=k+1
              if ((i==j).and.(i1==j1)) cycle
              X=ABS(M(1,i,i1,j,j1))
              Y=ABS(M(2,i,i1,j,j1))
              Z=ABS(M(3,i,i1,j,j1))
              RAVE=(X+Y+Z)/3.0_wp
              IF(RAVE < RAVEMIN) RAVEMIN=RAVE
              IF(RAVE > RAVEMAX) RAVEMAX=RAVE
            END DO
          END DO
        END DO
      END DO

      ! generate the GNUPLOT script in the $WorkDir
      line1=' '
      line2=' '
      gnuplot_CMD=' '
      file_exist=.false.
      is_file_open=.false.
      file_size=0
      execute_gnuplot_cmd =.false.
      gnuplot_version=0.0_wp

      ! check if the file lineOUT exists
      INQUIRE(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open,
     &        NUMBER=file_number)

      IF (file_exist) then
         IF (dbg) WRITE (StdOut,'(A)') 'file "lineOUT" exists in '//
     &                                 'WorkDir'
         IF (is_file_open) then
           IF (dbg) WRITE (StdOut,'(A)') 'file  "lineOUT" is opened'
           ! close the file:
           CLOSE (UNIT=file_number,STATUS='DELETE')
         END IF
         ! delete the file
         IF (dbg) WRITE (StdOut,'(A)') 'deleting the file...'
         iErr=AixRm("lineOUT")
         IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
      ELSE
         IF (dbg) WRITE (StdOut,'(A)') 'file "lineOUT" does not exist'//
     &                                 ' in WorkDir'
      END IF


      ! find the gnuplot
      IF (dbg) WRITE (StdOut,'(A)') 'inquire which GNUPLOT'

!#ifdef __INTEL_COMPILER
      CALL systemf ( "which gnuplot >> lineOUT", iErr )
      IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
!#else
!      CALL execute_command_line ( "which gnuplot >> lineOUT" )
!#endif

      INQUIRE(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open,
     &        NUMBER=file_number,SIZE=file_size)

      IF (dbg) WRITE (StdOut,*) 'File_number =',file_number
      IF (dbg) WRITE (StdOut,*) 'Is_file_open=',is_file_open
      IF (dbg) WRITE (StdOut,*) 'File_exist  =',file_exist
      IF (dbg) WRITE (StdOut,*) 'File_size   =',file_size

      IF (file_exist) then
        IF (file_size>0) then

          IF (dbg) WRITE (StdOut,'(A)') 'new file  "lineOUT"  exists'//
     &                                  ' in WorkDir'

          file_number=IsFreeUnit(45)
          Call molcas_open(file_number,"lineOUT")

          READ (file_number,'(A)') line1

          IF (dbg) WRITE (StdOut,*) 'line1=',line1
          IF (dbg) WRITE (StdOut,*) trim(line1)
          line2=trim(line1)
          IF (dbg) WRITE (StdOut,*) 'line2=',line2

          CLOSE(file_number)
          IF (dbg) WRITE (StdOut,*) 'Closing lineOUT file'
          FLUSH (StdOut)
          execute_gnuplot_cmd =.true.
        ELSE
          ! file_size =0
          WRITE (StdOut,'(A)') 'file  "lineOUT" has a size=0. '//
     &                         'gnuplot was not found on the system.'
          WRITE (StdOut,'(A)') 'plots will not be created.'
        END IF
      ELSE
         WRITE (StdOut,'(A)') 'file "lineOUT" does not exist in WorkDir'
      END IF
      ! remove file "lineOUT"
      iErr=AixRm("lineOUT")
      IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr


!!!!! check the version of the gnuplot:
      IF ( execute_gnuplot_cmd ) Then
        IF (dbg) WRITE (StdOut,'(A)') 'inquire which version of '//
     &                                'GNUPLOT is installed'
        ! attempt to execute the script
        WRITE (gnuplot_CMD,'(2A)') trim(line2),' --version > lineOUT'
        IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
!#ifdef __INTEL_COMPILER
        CALL systemf ( gnuplot_CMD, iErr )
        IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
!#else
!        CALL execute_command_line ( gnuplot_CMD )
!#endif
        file_number=IsFreeUnit(42)
        Call molcas_open(file_number,'lineOUT')
        READ (file_number,*) cdummy, gnuplot_version
        IF (dbg) WRITE (StdOut,'(A,F4.1)') 'gnuplot_version = ',
     &           gnuplot_version
        IF (abs(gnuplot_version)<0.1_wp) execute_gnuplot_cmd =.false.
        CLOSE (file_number)
        ! remove file "lineOUT"
        iErr=AixRm("lineOUT")
        IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
      END IF

!--------------------------------------------------------------------------------------------
! get the true real names of the files on disk:
      WRITE(datafile_e,'(3A)') 'BARRIER_ENE.dat'
      WRITE(datafile_m,'(3A)') 'BARRIER_TME.dat'
      WRITE(imagefile,'(3A)')  'BARRIER.png'
      WRITE(epsfile  ,'(3A)')  'BARRIER.eps'
      WRITE(plotfile ,'(3A)')  'BARRIER.plt'
      Call prgmtranslate(datafile_e ,realname_e_dat,Length)
      Call prgmtranslate(datafile_m ,realname_m_dat,Length)
      Call prgmtranslate(imagefile  ,realname_png,Length)
      Call prgmtranslate(epsfile    ,realname_eps,Length)
      Call prgmtranslate(plotfile   ,realname_plt,Length)
      IF (dbg) THEN
        WRITE(StdOut,'(3A)') 'realname_e_dat=',trim(realname_e_dat)
        WRITE(StdOut,'(3A)') 'realname_m_dat=',trim(realname_m_dat)
        WRITE(StdOut,'(3A)') 'realname_png=',trim(realname_png)
        WRITE(StdOut,'(3A)') 'realname_eps=',trim(realname_eps)
        WRITE(StdOut,'(3A)') 'realname_plt=',trim(realname_plt)
      END IF




!!!!! prepare the energy data file:
      INQUIRE(FILE=datafile_e,EXIST=file_exist,OPENED=is_file_open,
     &        NUMBER=file_number)
      IF(file_exist)  iErr=AixRm(trim(datafile_e))
      IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
      LuData=IsFreeUnit(75)
      Call molcas_open(LuData,datafile_e)
      IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(datafile_e)//'" file'
      IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(realname_e_dat)//
     &                          '" file'
      FLUSH(StdOut)
      ! write energies
      l=0
      DO i=1,nMult
        DO i1=1,ndim(i)
          l=l+1
          xstart =  dble(M(3,i,i1,i,i1))-dlt
          xend   =  dble(M(3,i,i1,i,i1))+dlt

          WRITE (LuData,'(3ES24.14)') xstart, E(l)/fact
          WRITE (LuData,'(3ES24.14)') xend  , E(l)/fact
          WRITE (LuData,'(A)') ' '
        END DO
      END DO
      IF (dbg) WRITE (StdOut,*) 'Writing into the "'//trim(datafile_e)
     &                           //'" file'
      IF (dbg) WRITE (StdOut,*) 'Writing into the "'//
     &                           trim(realname_e_dat)//'" file'
      CLOSE (LuData)
      IF (dbg) WRITE (StdOut,*) 'Closing the "'//trim(datafile_e)
     &                           //'" file'
      IF (dbg) WRITE (StdOut,*) 'Closing the "'//
     &                           trim(realname_e_dat)//'" file'
      FLUSH (StdOut)

!!!!!================================================================
!!!!! prepare the transition matrix elements file:
      INQUIRE(FILE=datafile_m,EXIST=file_exist,OPENED=is_file_open,
     &        NUMBER=file_number)
      IF(file_exist)  iErr=AixRm(trim(datafile_m))
      IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
      LuData=IsFreeUnit(76)
      Call molcas_open(LuData,datafile_m)
      IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(datafile_m)//'" file'
      IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(realname_m_dat)//
     &                          '" file'
      FLUSH(StdOut)
      ! write energies,
      l=0
      DO i=1,nMult
        DO i1=1,ndim(i)
          l=l+1
          k=0
          DO j=1,nMult
            DO j1=1,ndim(j)
              k=k+1

              if ((i==j).and.(i1==j1)) cycle

              X=ABS(M(1,i,i1,j,j1))
              Y=ABS(M(2,i,i1,j,j1))
              Z=ABS(M(3,i,i1,j,j1))
              RAVE=(X+Y+Z)/3.0_wp

              xstart = dble( M(3,i,i1,i,i1) )
              xend   = dble( M(3,j,j1,j,j1) )

              ystart = E(l)/fact
              yend   = E(k)/fact

              fmtx='(ES24.14,ES24.14,ES24.14)'
              IF(RAVE > RAVEMAX*0.005_wp) THEN
                 WRITE (LuData,fmtx) xstart,ystart,RAVE
                 WRITE (LuData,fmtx) xend  ,yend  ,RAVE
                 WRITE (LuData,'(A)') ' '
              END IF
            END DO
          END DO
        END DO
      END DO
      IF (dbg) WRITE (StdOut,*) 'Writing into the "'//trim(datafile_m)
     &                           //'" file'
      CLOSE (LuData)
      IF (dbg) WRITE (StdOut,*) 'Closing the "'//trim(datafile_m)
     &                           //'" file'
      FLUSH (StdOut)

!!!!!===============================================================
!!!!! generate the GNUPLOT script in the $WorkDir



      INQUIRE(FILE=plotfile,EXIST=file_exist,OPENED=is_file_open,
     &        NUMBER=file_number)
      IF(file_exist) iErr=AixRm(trim(plotfile))
      IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
      LuPlt=IsFreeUnit(85)
      Call molcas_open(LuPlt,plotfile)
      IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(plotfile)//'" file'
      IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(realname_plt)//
     &                          '" file'

      IF ( gnuplot_version < 5.0_wp ) Then
      !===  GNUPLOT VERSION 4 and below ==>>  generate EPS
         WRITE (LuPlt,'(A)') 'set terminal postscript eps enhanced '//
     &                       'color  size 3.0, 2.0 font "arial, 10"'
         WRITE (LuPlt,'(A)') 'set output "'//trim(realname_eps)//'" '
         WRITE (LuPlt,'(A)') 'set grid'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# Axes'
         WRITE (LuPlt,'(A)') '  set xlabel "Momentum / '//
     &                       '{/Symbol m}_{B}" font "Arial,12"'
         IF( nint(fact).eq.1) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy / '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.10) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x10) / '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.100) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x100) / '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.1000) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x1000)/ '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.10000) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x10000)/ '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.100000) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x100000) / '//
     &                         'cm^{-1}" font "Arial,12"'
         END IF
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A,F24.14,A,F24.14,A)')
     &                              '  set xrange[', xmin,':', xmax,']'
         WRITE (LuPlt,'(A,F24.14,A,F24.14,A)')
     &                              '  set yrange[', Emin,':', Emax,']'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# Tics for axes'
         WRITE (LuPlt,'(A)') '  set xtics nomirror'
         WRITE (LuPlt,'(A)') '  set ytics nomirror'
         WRITE (LuPlt,'(A)') '  set mxtics'
         WRITE (LuPlt,'(A)') '  set mytics'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# Margins'
         WRITE (LuPlt,'(A)') '  set bmargin at screen 0.20'
         WRITE (LuPlt,'(A)') '  set lmargin at screen 0.15'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# legend'
         WRITE (LuPlt,'(A)') '  set key off' !right bottom'
         WRITE (LuPlt,'(A)')

         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A,ES24.14,A,ES24.14,A)')
     &                 'set palette defined ( ',RAVEMIN," '#FFF0F0', ",
     &                                          RAVEMAX," 'red' )"

         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# actual plotting'
         WRITE (LuPlt,'(A)') 'plot  "'//trim(realname_m_dat)//
     &                       '" using 1:2:3'//
     &                       ' with lines lt 1  lw 1 lc palette , \'
         WRITE (LuPlt,'(A)') '      "'//trim(realname_e_dat)//
     &                       '" using 1:2  '//
     &                       ' with lines lt 1  lw 8 lc rgb "black" '
         WRITE (LuPlt,'(A)')


      ELSE IF ( (gnuplot_version >= 5.0_wp).AND.
     &          (gnuplot_version <  6.0_wp) ) Then
      !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
         WRITE (LuPlt,'(A)') 'set terminal pngcairo transparent '//
     &         'enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
         WRITE (LuPlt,'(A)') 'set output "'//trim(realname_png)//'" '
         WRITE (LuPlt,'(A)') 'set grid'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# Axes'
         WRITE (LuPlt,'(A)') '  set xlabel "Momentum / '//
     &                       '{/Symbol m}_{B}" font "Arial,12"'
         IF( nint(fact).eq.1) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy / '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.10) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x10) '//
     &                         '/ cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.100) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x100) '//
     &                         '/ cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.1000) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x1000)/ '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.10000) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x10000)/ '//
     &                         'cm^{-1}" font "Arial,12"'
         ELSE IF ( nint(fact).eq.100000) THEN
           WRITE (LuPlt,'(A)') '  set ylabel "Energy (x100000) / '//
     &                         'cm^{-1}" font "Arial,12"'
         END IF
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A,F24.14,A,F24.14,A)')
     &                              '  set xrange[', xmin,':', xmax,']'
         WRITE (LuPlt,'(A,F24.14,A,F24.14,A)')
     &                              '  set yrange[', Emin,':', Emax,']'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# Tics for axes'
         WRITE (LuPlt,'(A)') '  set xtics nomirror'
         WRITE (LuPlt,'(A)') '  set ytics nomirror'
         WRITE (LuPlt,'(A)') '  set mxtics'
         WRITE (LuPlt,'(A)') '  set mytics'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# Margins'
         WRITE (LuPlt,'(A)') '  set bmargin at screen 0.20'
         WRITE (LuPlt,'(A)') '  set lmargin at screen 0.15'
         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# legend'
         WRITE (LuPlt,'(A)') '  set key off' !right bottom'
         WRITE (LuPlt,'(A)')

         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A,ES24.14,A,ES24.14,A)')
     &                  'set palette defined ( ',RAVEMIN," '#FFF0F0', ",
     &                                           RAVEMAX," 'red' )"

         WRITE (LuPlt,'(A)')
         WRITE (LuPlt,'(A)') '# actual plotting'
         WRITE (LuPlt,'(A)') 'plot  "'//trim(realname_m_dat)//
     &                       '" using 1:2:3 '//
     &                       'with lines lt 1  lw 1 lc palette , \'
         WRITE (LuPlt,'(A)') '      "'//trim(realname_e_dat)//
     &                       '" using 1:2   '//
     &                       'with lines lt 1  lw 8 lc rgb "black" '
         WRITE (LuPlt,'(A)')
      ELSE
        WRITE (StdOut,*) 'GNUPLOT has version: ', gnuplot_version
        WRITE (StdOut,*) 'This version of GNUPLOT is not known and '//
     &                   'thus, unsupported.'
        execute_gnuplot_cmd=.false.
      END IF
      CLOSE (LuPlt)

      IF (execute_gnuplot_cmd) Then
        ! attempt to execute the script
        WRITE (gnuplot_CMD,'(5A)') trim(line2),' < ',trim(realname_plt)
        IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
!#ifdef __INTEL_COMPILER
        CALL systemf ( gnuplot_CMD, iErr )
        IF(dbg)  WRITE (StdOut,*) 'iErr = ',iErr
!#else
!        CALL execute_command_line ( gnuplot_CMD )
!#endif

        IF ( gnuplot_version < 5.0_wp ) Then
          INQUIRE(FILE=trim(realname_eps),EXIST=file_exist,
     &            OPENED=is_file_open,NUMBER=file_number)
          IF(file_exist) THEN
            WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//
     &                           '" was created in Working directory.'
          ELSE
            WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//
     &                       '" was NOT created in Working directory.'
          END IF
        ELSE
          INQUIRE(FILE=trim(realname_png),EXIST=file_exist,
     &            OPENED=is_file_open,NUMBER=file_number)
          IF(file_exist) THEN
            WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//
     &                           '" was created in Working directory.'
          ELSE
            WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//
     &                       '" was NOT created in Working directory.'
          END IF
        END IF
      END IF


      dbg=.false.

      RETURN
#ifdef _WARNING_WORKAROUND_
      IF (.FALSE.) CALL UNUSED_CHARACTER(cdummy)
#endif
      END SUBROUTINE plot_barrier
