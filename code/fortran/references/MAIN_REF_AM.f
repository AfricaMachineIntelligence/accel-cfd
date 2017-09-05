C             *** ABHIJIT MUKHERJEE ***
      PROGRAM MAIN
!
                IMPLICIT NONE
                INCLUDE "VAR3D.DEC"
C             ***LOCAL VARIABLE***
                INTEGER I,J,K
!             
!*******LOOP VARIABLE***********************************
!
!              INTEGER IS
!*******************************************************
!             
                PI=ACOS(-1.E00)
                CALL INPUT3D
C             PRINT *, "0"
                CALL SETUP3D
C             PRINT *, "1"
                CALL SETLS3D
!
                CALL INI3D
                ITS=0
                TIME=0.
                DTSET=DT
                IF (START.EQ.1) THEN
                CALL INREAD3D
                ENDIF
!
!              CALL PBUB3D                                                                                                                                                                     !########
!
C             PRINT *, "2"
                CALL SETHF3D
C             PRINT *, "3"
                CALL SETPROP3D
C             PRINT *, "4"
                CALL SETK3D
C             PRINT *, "5"
C             PRINT *, "HELLO1", REL
               
C             *** FOR TECPLOT ***
                OPEN (60, FILE='hfluxS.dat')
                WRITE (60,*) 'VARIABLES="X(mm)","Z(mm)","NuwS","TW",
     6 "T2","LSW","LS1","LS2","NumS"'
                OPEN (62, FILE='hfluxN.dat')
                WRITE (62,*) 'VARIABLES="X(mm)","Z(mm)","NuwN","TW",
     6 "T2","LSW","LS1","LS2","NumN"'
                OPEN (64, FILE='hfluxT.dat')
                WRITE (64,*) 'VARIABLES="X(mm)","Y(mm)","NuwT","TW",
     6 "T2","LSW","LS1","LS2","NumT"'
                OPEN(70, FILE='ebal.dat')
                WRITE(70,*) 'VARIABLES="IT","Time(ms)","Qw","Qi",
     6"Qm","Qwt","Qit","Qv","Ql","Qd","Qie","Qic"'
                OPEN(75, FILE='nu.dat')
                WRITE(75,*) 'VARIABLES="IT","Time(ms)","NuWS","NMS",
     6 "NuWN","NMN","NuWT","NMT","NuT","NuC",'
                OPEN(80, FILE='vol.dat')
                WRITE(80,*) 'VARIABLES="IT","Time(ms)","ED(mm)","BD(mm)","L",
     6 "L1","L2","R0","P","K","S","T"'
                OPEN(82, FILE='volmult.dat')
                WRITE(82,*) 'VARIABLES="IT","Time(ms)","ED1","BD1","ED2","BD2"
     6 ,"ED3","BD3"'
                OPEN(85, FILE='error.dat')
                WRITE(85,*) 'VARIABLES="IT","Error","Erroravg","IS"'
!
                OPEN (87, FILE='force.dat')
                WRITE (87,*) 'VARIABLES="ITS","TIME","PXF","PXB",
     6 "PYF","PYB","PZF","PZB","PX","PY","PZ"'
!
!
                OPEN(91, FILE='utdma.dat')
                WRITE(91,*) 'VARIABLES="IT","IS","ITER",
     6         "RP0","RPK","RPR","GAMMAP"'
                OPEN(92, FILE='vtdma.dat')
                WRITE(92,*) 'VARIABLES="IT","IS","ITER",
     6         "RP0","RPK","RPR","GAMMAP"'
                OPEN(93, FILE='wtdma.dat')
                WRITE(93,*) 'VARIABLES="IT","IS","ITER",
     6         "RP0","RPK","RPR","GAMMAP"'
                OPEN(94, FILE='ptdma.dat')
                WRITE(94,*) 'VARIABLES="IT","IS","ITER",
     6         "RP0","RPK","RPR","GAMMAP"'
                OPEN(89, FILE='SRptdma.dat')
                WRITE(89,*) 'VARIABLES="IT","IS","ITER",
     6         "RP0","RPK","RPR","GAMMAP"'
                OPEN(95, FILE='ttdma.dat')
                WRITE(95,*) 'VARIABLES="IT","IS","ITER",
     6         "RP0","RPK","RPR","GAMMAP"'
                OPEN(96, FILE='reini.dat')
                WRITE(96,*) 'VARIABLES="IT","IS","ITREINI",
     6         "MITER","EMAX","EAVG"'
                OPEN(97, FILE='time.dat')
                WRITE(97,*) 'VARIABLES="IT","IS","VELMAX",
     6         "DT","DTMIN"'
                OPEN(98, FILE='mbal.dat')
                WRITE(98,*) 'VARIABLES="IT","Time(ms)","DV"
     6,"DBASE","ML1(mg)","ML2(mg)","DML(mg)","DMV(mg)"'
!              OPEN(99, FILE='newptdma.dat')
!              WRITE(99,*) 'VARIABLES="IT","IS","ITER",
!     6       "RP0","RPK","RPR","GAMMAP"'
!
!              OPEN(89, FILE='tdma.dat')
!              WRITE(89,*) 'VARIABLES="X","Y","Z","AE","AW","AN",
!     6       "AS","AT","AB","AP","BB","PHI"'
 
!
!
C             *************************************************
 
C             *** OUTPUT FOR FIRST TIMESTEP ***
                IT=ITS
                CALL EVAPORATION3D
                CALL OUTPUT3D
                CALL OUTINT13D                                                                                             
                CALL CHECKVOL3D
                CALL CHECKVOLMULT3D
                CALL MICRO3D
                CALL EBALANCE3D
                CALL CHECKHT3D
C             **************
                IF (START.EQ.1) THEN
                                CALL INTVEL3D
                                CALL LSET3D
                                CALL REINI3D
                ENDIF
                               
C             *************************************************
                ITS=ITS+1
                DO IT=ITS,ITMAX
                TIME=TIME + DT
!                                             
                CALL SETHF3D
                CALL SETPROP3D
                CALL SETK3D
                CALL MICRO3D
C
!              PRINT *, "HELLO2", REL
!              *******SIMPLER-ALGORITHM*******************************
                                ERROROLD=0.
                                DO IS=1,ITERMAX
!                              *** SIMPLER SUBROUTINE ***
                                                CALL SIMPLER3D              ! REMEMBER NO PRESSURE CORRECTION
!                              **************************
                !                              PRINT *, "HELLO3", REL
                                                CALL UVELOCITY3D
                !                              PRINT *, "HELLO3A", REL
                                                CALL VVELOCITY3D
                !                              PRINT *, "HELLO3B", REL
                                                CALL WVELOCITY3D
                !                              PRINT *, "HELLO3C", REL
                                                CALL UPDATEUVEL3D
                !                              PRINT *, "HELLO3D", REL
                                                CALL UPDATEVVEL3D
                !                              PRINT *, "HELLO3E", REL              
                                                CALL PRESSURE3D
                !                              PRINT *, "HELLO3F", REL
                !                              CALL PCORR3D
                !                              PRINT *, "HELLO3G", REL
                                                CALL VELOCITYCORR3D
                !                              PRINT *, "HELLO3H", REL
                                                CALL TEMP3D
                !                              PRINT *, "HELLO3I", REL
C                             PRINT *, "HELLO3"
!
                                WRITE(*,*) 'TS',IT,ERROR,IE,JE,KE,IS
!
                                ERRR=ABS(ERROR-ERROROLD)/ERROR
                                IF ((ERROR < EPS).OR.(ERRR<0.01)) EXIT
                                ERROROLD=ERROR
!
                                END DO ! Simple-loop
C
!                              ***ASSIGN VELSTAR TO VELSTGR ************
                                CALL FINALVEL3D
!                              ***UPDATE TEMPERATURE******************
                                CALL UPDATETEMP3D
!********************************************************
                                WRITE(85,8500) IT,ERROR,ERRORAVG,IS
8500                       FORMAT(I6,2(E20.8),I6)
C
                                WRITE(*,*) 'TS',IT,ERROR,IE,JE,KE,IS
!
!
                IF (MOD(IT,OT/5) .EQ. 0) THEN
C             *** VOLUME CHECK ***
                                CALL CHECKVOL3D
                                CALL CHECKVOLMULT3D
C             *** ENERGY BALANCE ***
                                CALL EBALANCE3D
                                CALL CHECKHT3D
C             *** FORCE BALANCE ***
                                CALL FORCE3D                                                  
!
                                OPEN(50, FILE='stop')
                                READ(50,*)
                                READ(50,*) STOP, ITERMAX, DTSET, PBUB
                                CLOSE(50)
                                                IF (STOP.NE.0) THEN
                                                CALL OUTPUT3D
                                                CALL OUTINT13D                                             
                                                EXIT
                                                END IF
                ENDIF
!
                IF (MOD(IT,OT) .EQ. 0) THEN
                                WRITE(*,*) 'TS',IT,ERROR,IE,JE,KE,IS
                                CALL OUTPUT3D
                                CALL OUTINT13D                                                                                             
                ENDIF
!
C             **********************************************************
C
!              **************************
C                             PRINT *, "HELLO4a"
                                CALL INTVEL3D                                 
                !              PRINT *, "HELLO4b"
                                CALL LSET3D                                                                                                      
                !              PRINT *, "HELLO4c"
                                CALL REINI3D                                                                                    
                !              PRINT *, "HELLO4d"
C                             ***************************
                                                END DO ! Time-loop
!
                END
