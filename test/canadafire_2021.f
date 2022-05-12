! Main program to run the stuff in canadafire.f
      PROGRAM run_canadafire

        integer, dimension(:), allocatable :: monthin,dayin
        real, dimension(:), allocatable :: tin,hin,win,rin
        real, dimension(:), allocatable :: buiout,fmout,
     &       isiout,fwiout,dmcout,dcout,dsrout

!      integer, parameter :: ntime = 49
      integer, dimension(12) :: lmon
      real, dimension(12) :: el,fl

       ! Hardcode input filename
C Read input from STDIN
C      OPEN(UNIT=1,FILE='f32in.dat',status='old',action='read')
      OPEN(UNIT=2,FILE='f32out_2021.dat',
     &    status='replace',action='write')
      
      DATA LMON /31,28,31,30,31,30,31,31,30,31,30,31/
      DATA EL /6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,
     *6.0/
      DATA FL /-1.6,-1.6,-1.6,.9,3.8,5.8,6.4,5.0,2.4,.4,-1.6,-1.6/

      IUNIT=0    ! Hardcode to metric units
C
C     20 - METRIC FORMAT, 25 - ENGLISH FORMAT
C
   20 FORMAT(F4.1,2I4,F5.1)
   25 FORMAT(F4.0,2I4,F5.2)
C
C     READS IN STATION & YEAR.
C
      READ(5,30) TITLE
   30 FORMAT(20A4)

C
C     READS STARTING MONTH OF THE YEAR AND NUMBER OF DAYS IN STARTING
C     MONTH.
C
   55 READ(5,60) M, NDAYS, NTIME 
   60 FORMAT(I1,I2,I4)

      allocate(monthin(ntime))
      allocate(dayin(ntime))
      allocate(tin(ntime))
      allocate(hin(ntime))
      allocate(win(ntime))
      allocate(rin(ntime))
      allocate(buiout(ntime))
      allocate(fmout(ntime))
      allocate(isiout(ntime))
      allocate(fwiout(ntime))
      allocate(dmcout(ntime))
      allocate(dcout(ntime))
      allocate(dsrout(ntime))


      WRITE(2,65)
   65 FORMAT(1H1'PROGRAM NO.: F-32')
C      CALL DATE(DAT)
C      WRITE(2,70) DAT,TITLE
   70 FORMAT(1X,9A1///1X,20A4//)
!      WRITE(2,75) FO,PO,DOT
!   75 FORMAT(' INITIAL VALUES FOR FFMC:',F5.1,', DMC:',F5.1,', DC:',F5.1,
!     *//)
 
      jday=1
      DO J=M,12
        NN=LMON(J)
        IF(J.EQ.M) then
          IDAYS=LMON(J)-NDAYS+1
        ELSE
          IDAYS=1
        end if

        IAST=1
        DO I=IDAYS,NN
          READ(5,20,END=295) T,IH,IW,R
          W=IW
          TX=T
          H=IH
          RAIN=R

          monthin(jday) = J
          dayin(jday) = i
          tin(jday) = T
          hin(jday) = h
          win(jday) = w
          rin(jday) = rain
          jday = jday + 1
        end do
      end do

 290  print *,'done 1'
! Now we write our own stuff!
 295  print *,'Done reading2'

      print *,monthin
      call canadafire(monthin,dayin,
     &     tin,hin,win,rin,
     &     buiout,fmout,isiout,fwiout,dmcout,dcout,dsrout)

      WRITE(2,105)
  105 FORMAT(///,1X1X,' DATE  TEMP  RH  WIND   RAIN   FFMC   DMC',
     *   '     DC   ISI   BUI   FWI     DSR'/)

      print *,monthin
      do jday=1,size(monthin)
        j=monthin(jday)
        i=dayin(jday)
        TX=tin(jday)
        IH=hin(jday)
        IW=win(jday)
        RAIN=rin(jday)
        ffm=fmout(jday)
        DMC=dmcout(jday)
        DC=dcout(jday)
        SI=isiout(jday)
        BUI=buiout(jday)
        FWI=fwiout(jday)
        DSR=dsrout(jday)
        WRITE(2,285) J,I,TX,IH,IW,RAIN,FFM,DMC,DC,SI,BUI,FWI,DSR
  285   FORMAT(1X,2I3,F6.1,I4,I6,F7.1,F7.1,F6.1,F7.1,3F6.1,F8.2) 
      end do

      CONTAINS

      include 'canadafire.f'


      end program
