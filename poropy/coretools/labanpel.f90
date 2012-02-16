! laban-pel
! modified to read command line for process rank (gives unique unit
! numbers for a process), and i/o files.
!                                        
      PROGRAM LABAN 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      REAL SINGL
      double precision  INITTIME
      integer u_inp, u_out
!      ------------------                                               
      PARAMETER (MEMMA = 250000, MEMMH = 500, MEMMS = 100) 
!      ------------------                                               
      DIMENSION A (MEMMA), HELL (MEMMH), SINGL (MEMMS) 
      COMMON / ADRS / AD (50) 
      COMMON / CLAB / WW (100) 
      COMMON / SING / MEMH, MEMS 
      COMMON / PWR / EPF, CPOW, VCORE 
      COMMON / VERSI / IVERS 

      EQUIVALENCE (WW (13), NQ), (WW (34), ITER), (WW (42), XKEF),      &
      (WW (60), IFLUX), (WW (64), IFIL), (WW (70), NQTRUE), (WW (75),   &
      MAVAIL), (WW (78), IBUCK), (WW (79), JFIL)                        

      MAVAIL = MEMMA 
      MEMH = MEMMH 
      MEMS = MEMMS 
!                                                                       
! IVERS =0 -> IBM MAINFRAME VERSION                                     
!       =1 -> 386 PC VERSION (XENIX/LPI FORTRAN;  DOS/LAHEY FORTRAN)    
!PC                                                                     
!       IVERS=0                                                         
      IVERS = 1 
      INITTIME = WTIME()
!PC                                                                     
!                                                                       
!  OPEN FILES                                                           
!                                                                       
!      --------------------                                             
      IF (IVERS.EQ.1) CALL FILOP 

                                                                                                                                    
   10 CALL INPUT (SINGL, A, IEXIT, HELL) 

!      --------------------                                             
      IF (IEXIT.EQ.2) GOTO 200 
!         
                                                       
      ITER = 0 

      !>>>>>>>>>>>>>>>MAIN LOOP>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
                                                           !
   20 ITER = ITER + 1                                      ! 
      XKOLD = XKEF                                         !
                                                           !      
      CALL RESCAL (A, 1)                                   !
      CALL ITERA (A, IEXIT, HELL)                          !
                                                           !  
      IF (IEXIT.EQ.2) GOTO 20                              !
                                                           !
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
                     
      CALL OUTPUT (A, 0, HELL) 

!      --------------------                                             
      IF (IFLUX.EQ.0.AND.IBUCK.EQ.0) GOTO 10 
!                                                                       
      XKEF = XKOLD 
!      --------------------                                             
      CALL RESCAL (A, NQTRUE) 
      CALL OUTPUT (A, 1, HELL) 
!      --------------------                                             
      GOTO 10 
  200 CONTINUE 

      !IF (IVERS.EQ.1) THEN 
      !   CLOSE (5) 
      !   CLOSE (6) 
      !   CLOSE (IFIL) 
      !   CLOSE (JFIL) 
      !ENDIF 
      !write (0,*) "ELAPSED TIME = ", WTIME() - INITTIME, " SECONDS."
      STOP 
      END PROGRAM LABAN    

function wtime ( )
!*****************************************************************************
!
! WTIME returns a reading of the wall clock time.
!
!  Discussion:
!
!    To get the elapsed wall clock time, call WTIME before and after a given
!    operation, and subtract the first reading from the second.
!
!    This function is meant to suggest the similar routines:
!
!      "omp_get_wtime ( )" in OpenMP,
!      "MPI_Wtime ( )" in MPI,
!      and "tic" and "toc" in MATLAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) WTIME, the wall clock reading, in seconds.
!
  implicit none
  integer           clock_max
  integer           clock_rate
  integer           clock_reading
  double precision  wtime
  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = dble(clock_reading) / dble(clock_rate)
  return
end

                         
!=======================================================================
!FILOP                       4 APR  1990         FILOP               LAB
!                                                                       
!                            *********                                  
!                            * FILOP *                                  
!                            *********                                  
!                                                                       
! OPENS THE INPUT, OUTPUT AND SCRATCH FILES FOR THE PC VERSION          
! NOTE THAT THE INSTALL FILE 'LABANPEL.INS' MUST EXIST                  
!                                                                       
!      --------------------------                                       
      SUBROUTINE FILOP 
!      --------------------------                                       
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      CHARACTER FILNAM * 32, FNAME * 32, BLANK * 32 
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (64), IFIL) 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / FILECM / FILNAM (4) 

      EQUIVALENCE (LUNFIL (1), INST) 
      EQUIVALENCE (LUNFIL (2), NPRNT) 
      EQUIVALENCE (LUNFIL (3), NIN) 
      EQUIVALENCE (LUNFIL (4), JFIL) 
      character(len=30)  :: inputfile, outputfile
      character(30)      :: tmp
      integer node, io

      if ( COMMAND_ARGUMENT_COUNT() .lt. 3 ) then
          print *, "*** ERROR: user input file not specified ***"
          return
      else
          node = 0
          call get_command_argument(1,tmp);
          read(tmp, '(i10)') node
          INST  = 1+5*node  
          NPRNT = 2+5*node  
          NIN   = 3+5*node 
          JFIL  = 0+5*node  
          IFIL  = 4+5*node 
          call get_command_argument(2,inputfile);
          open (unit = LUNFIL(3), file = inputfile, action = "read", &
                status = "old", position = "rewind", iostat = io)
          if (io > 0) stop "*** ERROR: user input file not found ***"
          call get_command_argument(3,outputfile);
          open (unit = LUNFIL(2), file = outputfile, action = "write", &
                status = "unknown", position = "rewind", iostat = io)
          if (io > 0) stop "*** ERROR: user output file not found ***"
      end if


      END SUBROUTINE FILOP                          
!INPUT                       4  DEC 1983         INPUT               LAB
!                                                                       
!                            *********                                  
!                            * INPUT *                                  
!                            *********                                  
!                                                                       
!      --------------------------                                       
      SUBROUTINE INPUT (SINGL, C, IEXIT, HELP) 
!      --------------------------                                       
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION KEFF 
      REAL SINGL ( * ) 
      DIMENSION C ( * ), HELP ( * ) 
      DIMENSION ITERM (4), ITITLE (16) 
      CHARACTER ITERM * 5, ITITLE * 5 
!                                                                       
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / ADRS / AD (50) 
      COMMON / DINV / IDIAG 
      COMMON / SING / MEMH, MEMS 
      COMMON / PWR / EPF, CPOW, VCORE 

!                                                                       
      EQUIVALENCE (AD (01), LNDX), (AD (02), LNDY), (AD (03), LHX),     &
      (AD (04), LHY), (AD (05), LCHI), (AD (06), LXSEC), (AD (07),      &
      LALB), (AD (08), LONE), (AD (09), LMAP), (AD (10), LCHECK),       &
      (AD (15), LF), (AD (16), LA), (AD (17), LT), (AD (18), LB),       &
      (AD (19), LRES), (AD (21), LOLD), (AD (22), LCUR), (AD (23),      &
      LFLUX), (AD (24), LPOW), (AD (25), LOUT), (AD (30), LHX0),        &
      (AD (31), LHY0), (AD (32), LMAP0), (AD (33), LALB0), (AD (50),    &
      LAST)                                                             
!                                                                       
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (03), NX2), (WW (04)&
      , NY2), (WW (05), NNX), (WW (06), NNY), (WW (07), NXY2), (WW (08),&
      NG), (WW (09), NGNG), (WW (10), NG4), (WW (11), NL), (WW (12),    &
      NLR), (WW (13), NQ), (WW (15), NGL), (WW (17), NGL4), (WW (18),   &
      NGLB), (WW (19), NGLA), (WW (20), NGLF), (WW (21), NGLR), (WW (22)&
      , NGLR2), (WW (24), LH1), (WW (25), LH2), (WW (27), MIXCAL),      &
      (WW (28), MIXALB), (WW (29), MIXTOT), (WW (30), KURTOT), (WW (31),&
      LUXTOT), (WW (32), MAXIN), (WW (33), MAXOUT), (WW (36), EPSR),    &
      (WW (37), EPSK), (WW (38), EPSI), (WW (41), KACC), (WW (42),      &
      KEFF), (WW (44), OMEGA), (WW (55), TIMIT), (WW (56), TIMER),      &
      (WW (57), TIMIN), (WW (58), TIMUT), (WW (59), ICUR), (WW (60),    &
      IFLUX), (WW (62), MCASE), (WW (63), NP1), (WW (66), NNX2),        &
      (WW (67), NNY2), (WW (70), NQTRUE), (WW (71), NFSIZ), (WW (72),   &
      NHSIZ), (WW (73), MIXXSE), (WW (74), MIXEXT), (WW (75), MAVAIL),  &
      (WW (77), IADJNT), (WW (78), IBUCK), (WW (79), JFIL)              
!                                                                       
      DATA ITERM / 'LABAN', 'NEWL ', 'NEWX ', 'END  ' / 
!                                                                       
!         READ INDATA                                                   
!                                                                       
      READ (LUNFIL(3), 1) (ITITLE (I), I = 1, 16) 
      DO 10 M = 1, 4 
         MCASE = M 
         IF (ITITLE (1) .EQ.ITERM (M) ) GOTO 12 
   10 END DO 
   12 IF (MCASE.NE.4) WRITE (LUNFIL(2), 8) 
      WRITE (LUNFIL(2), 2) (ITITLE (I), I = 1, 16) 
      IF (MCASE.EQ.4) GOTO 60 
!                                                                       
      IF (MCASE.EQ.1) READ (LUNFIL(3), * ) NG, NNX, NNY, MIXXSE, IDIAG, MIXEXT, &
      MIXALB, KEFF                                                      
!                                                                       
      READ (LUNFIL(3), * ) L, NLR, NQ, ICUR, IFLUX, IBUCK, NGRID, IADJNT, JFIL 
      READ (LUNFIL(3), * ) NNORM, CPOW, EPF, KACC, MAXIN, MAXOUT, EPSR, EPSI,   &
      EPSK, OMEGA                                                       
!                                                                       
!         SET DEFAULT VALUES                                            
!                                                                       
      IF (JFIL.EQ.5.OR.JFIL.EQ.6.OR.JFIL.EQ.10.OR.JFIL.LE.0) JFIL = 0 
      IF (MAXIN.EQ.0) MAXIN = 150 
      IF (MAXOUT.EQ.0) MAXOUT = 5 
      IF (KACC.EQ.1) OMEGA = 1.0D0 
      IF (CPOW.LT.1.D-1) CPOW = 0.D0 
      IF (EPF.LE.0.D0) EPF = 3.2D-11 
      IF (EPSR.LT.1.D-10) EPSR = 1.D-5 
      IF (EPSI.LT.1.D-10) EPSI = 1.D-5 
      IF (EPSK.LT.1.D-10) EPSK = 1.D-5 
      IF (OMEGA.LT.1.D-10) OMEGA = 1.0D0 
      IF (KEFF.LT.1.D-10) KEFF = 1.0D0 
      IF (NLR.LT.L.OR.NG.EQ.1) NLR = L 
      IF (MIXEXT.GT.0) NQ = 0 
      IF (MIXEXT.GT.0) IFLUX = 0 
      IF (NQ.EQ.0) NGRID = 0 
      IF (NG.GT.2) IBUCK = 0 
      IF (IDIAG.LT.0.OR.IDIAG.GT.2) IDIAG = 0 
      KDIAG = IDIAG 
      IF (IDIAG.EQ.2) KDIAG = 1 
!                                                                       
      WRITE (LUNFIL(2), 7) NG, MIXXSE, KEFF, IFLUX, L, MIXEXT, EPSR, ICUR, NLR, &
      MIXALB, EPSI, NGRID, NQ, EPSK, NNORM, NNX, MAXOUT, KACC, NNY,     &
      MAXIN, OMEGA, IDIAG, JFIL, IBUCK, CPOW, EPF                       
      IF (IADJNT.EQ.1) WRITE (LUNFIL(2), 2050) 
!                                                                       
!         READ NODE PARTITIONS IN X AND Y (SUBROUTINE SUBDIV)           
!                                                                       
      IF (MCASE.EQ.2) GOTO 15 
      LNDX = 1 
      LNDY = LNDX + NNX 
      LAST = LNDY + NNY 
!      -----------------------------                                    
      CALL SUBDIV (C (LNDX), C (LNDY) ) 
!      -----------------------------                                    
!                                                                       
!  COMMON VARIABLES                                                     
!                                                                       
   15 NG4 = NG + 6 + NG * KDIAG 
      NGNG = NG * NG 
      NL = L + 1 
      NQ = NQ + 1 
      NQTRUE = NQ 
      NLR = NLR + 1 
      NGL = NG * NL 
      NGL4 = 4 * NGL 
      NGLB = NGL * NGL 
      NGLA = NGNG * ( (NL * NL + 1) / 2) 
      NGLR = NG * NLR 
      NGLR2 = 2 * NGLR * NGLR 
      NX2 = NX + 2 
      NY2 = NY + 2 
      NXY2 = NX2 * NY2 
      NNX2 = NNX + 2 
      NNY2 = NNY + 2 
      NP1 = NGRID+1 
      MIXCAL = MIXXSE-MIXEXT 
      TIMER = 0.D0 
      TIMIT = 0.D0 
      TIMIN = 0.D0 
      TIMUT = 0.D0 
!                                                                       
!         ALLOCATION OF ARRAYS                                          
!                                                                       
      LHX0 = LNDY + NNY 
      LHY0 = LHX0 + NNX 
      LCHI = LHY0 + NNY 
      LXSEC = LCHI 
      LMAP0 = LXSEC + NG * NG4 * MIXXSE 
      LALB0 = LMAP0 + NNX2 * NNY2 + 1 
      LHX = LALB0 + MIXALB * NGNG + 1 
      LHY = LHX + NX2 
      LMAP = LHY + NY2 
      LALB = LMAP + 7 * NXY2 
      LONE = LALB + MIXALB * NGLB 
      LCHECK = LONE+NGLR * NGLR 
      LAST1 = LCHECK + 6 * NG 
!                                                                       
!         CALL INDATL FOR FURTHER DATA INPUT.  CALL MAPALL FOR DETERMI- 
!         NATION OF HOW RESPONSE MATRICES AND CURRENTS ARE TO BE STORED.
!                                                                       
      IF (MCASE.EQ.2) GOTO 20 
!      -------------------------------------------------------------    
      CALL INDATA (SINGL, C (LMAP0), C (LNDX), C (LNDY), C (LHX),       &
      C (LHY), C (LHX0), C (LHY0), C (LXSEC), NG, NG4, NNX2, HELP)      
      CALL MAPCOR (C (LMAP), C (LMAP0), C (LNDX), C (LNDY), NX2, NNX2) 
   20 CALL ALBINL (C (LALB0), C (LALB), NG, NGL) 
      CALL MAPALL (C (LMAP), C (LHX), C (LHY), NXY2) 
!      -------------------------------------------------------------    
!                                                                       
!  ALLOCATION OF ARRAYS                                                 
!                                                                       
      NQQ = NQ * NQ 
      NGQ = NG * NQQ 
      NGLF = NGNG * NQ * ( (NQ * NL + 1) / 2) 
      LQX = MAX (NLR, NQ) 
      NFSIZ = 2 * NGNG * NLR * NQQ 
      LH1 = NGNG * MAX (NLR * NLR, 6 + 2 * LQX + MAX (4, 2 + LQX) ) 
      LH2 = NGNG * NLR * MAX (NLR, NQQ) 
      NHSIZ = 2 * NGNG + 2 * NGLR2 + NFSIZ + LH1 + LH2 
      NRES = 4 * NGLR2 + NFSIZ + NHSIZ 
      NOUT = NP1 * NQ + 2 * NGQ + 3 * NGNG + 2 * NG + NP1 * NP1 
!                                                                       
      LF = LAST1 
      LA = LF + MIXTOT * NGLF 
      LT = LA + MIXTOT * NGLA 
      LB = LT + MIXTOT * NGLA 
      LRES = LB + MIXTOT * NGLB 
      LOLD = LRES + NRES 
      LCUR = LOLD+KURTOT 
      LFLUX = LCUR + KURTOT 
      LPOW = LFLUX + NG * LUXTOT 
      LOUT = LPOW + NX * NY 
      LAST = LOUT + NOUT + 2 * NG + NG 
!                                                                       
      MEMHR = NG * NL * 4 
      MEMSR = NG * NG4 
      WRITE (LUNFIL(2), 2100) LAST, MAVAIL, MEMHR, MEMH, MEMSR, MEMS 
      IF (LAST.LE.MAVAIL.AND.MEMHR.LE.MEMH.AND.MEMSR.LE.MEMS) GOTO 30 
      STOP 'LABAN-PEL FAILED - TOO LITTLE MEMORY' 
!                                                                       
!  CALL RERESP TO READ EXTERNAL RESPONSE MATRICES                       
!                                                                       
   30 IF (MIXEXT.EQ.0) GOTO 40 
!      ----------------------------------------------------------       
      CALL RERESP (NXY2, C (LMAP), C (LONE), C (LA), C (LT), C (LB),    &
      C (LF), C (LRES) )                                                
!      ----------------------------------------------------------       
!                                                                       
   40 WRITE (LUNFIL(2), 2110) 
      IEXIT = 1 
      GOTO 90 
   60 IEXIT = 2 
   90 RETURN 
!                                                                       
    1 FORMAT  (16A5) 
    2 FORMAT  (/1X,16A5/118('.')) 
    3 FORMAT  (16I5) 
    5 FORMAT  (7I5,2F10.0) 
    6 FORMAT  (4I5,6F10.0) 
    7 FORMAT (/,' NR OF GROUPS   NG =', I3,4X,' CROSS-SECTION MIX =',I3,&
     &4X,' K-EFF START GUESS =', F9.6,3X,' FLUX PRINTING     =', I3/    &
     & ' CURR APPR ORDER L =', I3,4X,' EXTERNAL RESP MIX =', I3,4X,     &
     & ' ACCU RESP  CALC   =', F9.6,3X,' CURRENT PRINTING  =', I3/      &
     & ' RESP CALC ORDER   =', I3,4X,' BOUND. ALBEDO MIX =', I3,4X,     &
     & ' ACCU POWER ITER   =', F9.6,3X,' NO OF GRID POINTS =', I3/      &
     & ' FLUX APPR ORDER   =', I3,4X,27X,' ACCU KEFF  ITER   =', F9.6,  &
     &3X,' NORMALIZATION     =', I3/                                    &
     & ' HORIS NR OF NODES =', I3,4X,27X,' MAX NR K-EFF ITER =', I4,8X, &
     & ' ACCELERATION OPT. =', I3/                                      &
     & ' VERTI NR OF NODES =', I3,4X,27X,' MAX NR POWER ITER =', I4,8X, &
     & ' OMEGA (SOR)       =', F9.6/                                    &
     & ' D-MATRIX OPTION   =', I3,4X,' POWER-MAP EXPORTED=', I3,4X,     &
     & ' FLUX IN BUCK. REP.=', I3/                                      &
     & ' CORE POWER =',1PE12.6,' WATTS'/                                &
     & ' ENERGY RELEASED PER FISSION =',1PE12.6,' (JOULES/FISS)')       
    8 FORMAT  (///,25X,'*************'/                                 &
     &             25X,'* LABAN-PEL *    RELEASED OCT 1991'/            &
     &             25X,'*************'/)                                
 2050 FORMAT  (//1X,23('=')/' ADJOINT SOLUTION SOUGHT'/1X,23('=')) 
!                                                                       
 2100 FORMAT  (//1X,100('*')/                                           &
     &        ' LABAN-PEL MEMORY REQUIREMENTS' //                       &
     &        ' MEMORY NEEDED IN A                =',I6/                &
     &        ' MEMORY AVAILABLE FOR A   (MAVAIL) =',I6/                &
     &        ' MEMORY NEEDED IN HELL             =',I6/                &
     &        ' MEMORY AVAILABLE FOR HELL  (MEMH) =',I6/                &
     &        ' MEMORY NEEDED IN SINGL            =',I6/                &
     &        ' MEMORY AVAILABLE FOR SINGL (MEMS) =',I6/)               
!                                                                       
 2110 FORMAT  (///1X,30('*'),'  INPUT FINISHED  ',30(1H*)) 
!                                                                       
      END SUBROUTINE INPUT                          
!SUBDIV                      26 MAI 1982         INPUTA              LAB
!                                                                       
!                            **********                                 
!                            * SUBDIV *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE SUBDIV (NDX, NDY) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION NDX ( * ), NDY ( * ) 
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      EQUIVALENCE (WW (1), NX), (WW (2), NY), (WW (5), NNX), (WW (6),   &
      NNY)                                                              
!                                                                       
!  READ ARRAYS FOR SUBDIVISION OF NODES.  DEFAULT VALUES = 1            
!                                                                       
      READ (LUNFIL(3), * ) (NDX (I), I = 1, NNX) 
      READ (LUNFIL(3), * ) (NDY (I), I = 1, NNY) 

      DO I = 1, NNX 
        IF (NDX (I) .EQ. 0) THEN
            write(0, *) "ndx(i) = ", NDX(I)
            NDX (I) = 1 
        END IF
      END DO

      DO 6 I = 1, NNY 
    6 IF (NDY (I) .EQ.0) NDY (I) = 1 
      WRITE (LUNFIL(2), 102) (NDX (I), I = 1, NNX) 
      WRITE (LUNFIL(2), 104) (NDY (I), I = 1, NNY) 
!                                                                       
!  CALCULATE TOTAL NUMBER OF NODES IN X- AND Y-DIRECTIONS.              
!                                                                       
      NX = 0 
      DO 10 I = 1, NNX 
   10 NX = NX + NDX (I) 
      NY = 0 
      DO 20 I = 1, NNY 
   20 NY = NY + NDY (I) 
      RETURN 
!                                                                       
  102 FORMAT  (//' NDX =',35I3) 
  104 FORMAT  ( /' NDY =',35I3) 
      END SUBROUTINE SUBDIV                         
!INDATA                      4  DEC 1983         INPUTA              LAB
!                                                                       
!                            **********                                 
!                            * INDATA *                                 
!                            **********                                 
!                                                                       
!      ------------------------------------------------------           
      SUBROUTINE INDATA (S, MAP, NDX, NDY, HX, HY, HX0, HY0, XSEC, MG,  &
      NG4, MNX2, HELP)                                                  
!      ------------------------------------------------------           
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      REAL S (1) 
      DIMENSION MAP (MNX2, * ), NDX ( * ), NDY ( * ), HX ( * ), HY ( * )&
      , XSEC (MG, NG4, * ), HX0 ( * ), HY0 ( * ), HELP ( * )            
      CHARACTER LTITLE * 32 
      COMMON / XSFIL / LTITLE 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / DINV / IDIAG 
!                                                                       
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (5), NNX), (WW (06),&
      NNY), (WW (62), MCASE), (WW (73), MIXXSE), (WW (79), JFIL)        
!                                                                       
!  ----------------                                                     
!  READ NODE SIZES.                                                     
!  ----------------                                                     
      NG = MG 
      NGG4 = NG * NG4 
      IF (MCASE.NE.1) GOTO 1 
      READ (LUNFIL(3), * ) (HX0 (I), I = 1, NNX) 
      READ (LUNFIL(3), * ) (HY0 (I), I = 1, NNY) 
      WRITE (LUNFIL(2), 102) (HX0 (I), I = 1, NNX) 
      WRITE (LUNFIL(2), 104) (HY0 (I), I = 1, NNY) 
!                                                                       
!  -----------------------------------------                            
!  CALCULATE NODE SIZES FOR SUBDIVIDED NODES.                           
!  -----------------------------------------                            
    1 JX = NX + 1 
      DO 4 J = 1, NNX 
         J0 = NNX + 1 - J 
         FNDX = NDX (J0) 
         DO 2 JJ = 1, NDX (J0) 
    2    HX (JX - JJ) = HX0 (J0) / FNDX 
    4 JX = JX - NDX (J0) 
!                                                                       
      IY = NY + 1 
      DO 6 I = 1, NNY 
         I0 = NNY + 1 - I 
         FNDY = NDY (I0) 
         DO 5 II = 1, NDY (I0) 
    5    HY (IY - II) = HY0 (I0) / FNDY 
    6 IY = IY - NDY (I0) 
      INDX = 0 
      HX (INDX) = 1.D0 
      HY (INDX) = 1.D0 
      HX (NX + 1) = 1.D0 
      HY (NY + 1) = 1.D0 
!                                                                       
!  ---------------------------                                          
!  READ MIXTURE ALLOCATION MAP.                                         
!  ---------------------------                                          
      IF (MCASE.NE.1) RETURN 
      WRITE (LUNFIL(2), 118) (J, J = 1, NNX + 1) 
      WRITE (LUNFIL(2), 119) 
      NNY2 = NNY + 2 
      DO 8 I = 1, NNY2
         IM = I - 1 
         READ (LUNFIL(3), 120) (MAP (J, I), J = 1, NNX + 2) 
    8 WRITE (LUNFIL(2), 122) IM, (MAP (J, I), J = 1, NNX + 2) 
!                                                                       
!  -------------------                                                  
!  READ CROSS SECTIONS.                                                 
!  -------------------                                                  
      WRITE (LUNFIL(2), 107) 
!                                                                       
      DO 20 K = 1, MIXXSE 
         READ (LUNFIL(3), 142) MIXTP, IFLXS, IDENT, LTITLE 
         IFILE = IABS (IFLXS) 
         IF (MIXTP.NE.K) THEN 
            WRITE (LUNFIL(2), 2050) K, MIXTP 
            STOP 'MIXTURES NOT IN CHRONOLOGICAL ORDER' 
         ENDIF 
         IF (IFILE.EQ.6.OR.IFILE.EQ.10.OR.IFILE.EQ.JFIL) THEN 
            WRITE (LUNFIL(2), 2040) JFIL 
            STOP 'FILE ERROR' 
         ENDIF 
         IF (IFLXS.EQ.5) THEN 
!                                                                       
!  READ CROSS SECTIONS FROM FILE 5 (SAME FILE AS OTHER INDATA)          
!                                                                       
            CALL XSREAD (XSEC (1, 1, K), IFLXS, NG, NG4) 
!                                                                       
         ELSE 
!                                                                       
!  READ CROSS SECTIONS FROM STANDARD/LABAN FILE  (BINARY)               
!                                                                       
            CALL INMIX (S, XSEC (1, 1, K), IDENT, IFLXS, NG, NGG4) 
!                                                                       
         ENDIF 
         IF (IFLXS.EQ.9999) GOTO 99 
!                                                                       
!  PRINT XSECTS                                                         
!                                                                       
         IF (IDIAG.EQ.0) WRITE (LUNFIL(2), 112) K, IDENT, LTITLE 
         IF (IDIAG.EQ.1) WRITE (LUNFIL(2), 113) K, IDENT, LTITLE 
         IF (IDIAG.EQ.2) WRITE (LUNFIL(2), 114) K, IDENT, LTITLE 
         DO 15 I = 1, NG 
            WRITE (LUNFIL(2), 110) (XSEC (I, J, K), J = 1, NG4) 
   15    END DO 
!                                                                       
!  CHECK FISSION SPECTRUM SUMS TO UNITY ELSE ADJUST GROUP 1 VALUE       
!                                                                       
         CALL CHIC (XSEC (1, 4, K), NG, K) 
!                                                                       
!  INVERT D MATRIX IF IDIAG=2                                           
!                                                                       
         IF (IDIAG.EQ.2) THEN 
            IADD = NG4 - NG + 1 
            CALL MINV (XSEC (1, IADD, K), NG, DET, HELP, HELP (1 + NG) ) 
         ENDIF 
!                                                                       
!  REPLACE INPUT DIFFUSION COEFFICIENTS USING COLUMN SUMS OF DINV       
!                                                                       
         IF (IDIAG.NE.0) THEN 
            IADD = NG4 - NG 
            DO 17 KGR = 1, NG 
               IAS = IADD+KGR 
               SUM = 0D0 
               DO 16 IGR = 1, NG 
                  SUM = SUM + XSEC (IGR, IAS, K) 
   16          END DO 
               XSEC (KGR, 1, K) = 1D0 / SUM 
   17       END DO 
         ENDIF 
   20 END DO 
      RETURN 
!                                                                       
!  ERROR RETURN                                                         
!                                                                       
   99 WRITE (LUNFIL(2), 2020) MIXTP, IDENT 
      STOP 'ERROR WHEN READING CROSS SECTIONS' 
!                                                                       
  102 FORMAT  ( /' HX  =',10F10.4,(/6X,10F10.4)) 
  104 FORMAT  ( /' HY  =',10F10.4,(/6X,10F10.4)) 
  107 FORMAT  (//' CROSS SECTIONS (PRINTED COLUNMWISE FOR D,ABS,NUFIS,CH&
     &I,SCAT,FIS,DISC,MMAT)'/1X,14('*')/)                               
  108 FORMAT  (8F10.0) 
  110 FORMAT  (6X,12F10.6) 
  112 FORMAT  (/' MIX',I4,4X,'IDENT =',I5,1X,'<<',A32,'>>',2X,' DIAGONAL&
     & D MATRIX'/ 1X,7('-'))                                            
  113 FORMAT  (/' MIX',I4,4X,'IDENT =',I5,1X,'<<',A32,'>>',2X,' FULL INV&
     &ERSE D MATRIX'/ 1X,7('-'))                                        
  114 FORMAT  (/' MIX',I4,4X,'IDENT =',I5,1X,'<<',A32,'>>',2X,' FULL D M&
     &ATRIX'/ 1X,7('-'))                                                
  118 FORMAT  (//' MIXTURE ALLOCATION' /1X,18('*') /'   I/J  0' ,39I3) 
  119 FORMAT  (1X) 
  120 FORMAT  (30I4) 
  122 FORMAT  (1X,I3,2X,40I3) 
  142 FORMAT  (3I5,2X,A32) 
 2020 FORMAT  (//1X,100('*')/' MIXTURE=',I3,5X,'IDENT=',I4,3X,          &
     &        'COULD NOT BE FOUND IN CROSS SECTION LIBRARY')            
 2040 FORMAT  (//1X,100('*')/' IFLXS SHOULD NOT EQUAL 6, 10 OR ',I2) 
 2050 FORMAT  (//1X,100('*')/' MIXTURE NUMBER',I3,' HAS MIXTURE ID',I3, &
     &        '.   THESE NUMBERS SHOULD BE THE SAME')                   
      END SUBROUTINE INDATA                         
!*                                                                      
!CHIC             ************                                          
!----EZM          *  CHIC    *         MAY 1991                         
!                 ************                                          
!                                                                       
!     CHECKS THAT FISSION SPECTRUM IS NORMALISED TO UNITY               
!     IF NOT, GROUP 1 ENTRY IS ADJUSTED                                 
!                                                                       
      SUBROUTINE CHIC (CHI, NG, MIX) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION CHI (NG) 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      SUM = 0D0 
      DO 10 I = 1, NG 
         SUM = SUM + CHI (I) 
   10 END DO 
      DIF = 1D0 - SUM 
      IF (DIF.NE.0) THEN 
         CHI (1) = CHI (1) + DIF 
         WRITE (LUNFIL(2), 1000) MIX, DIF 
      ENDIF 
      RETURN 
 1000 FORMAT (/1X,'WARNING!  FISSION SPECTRUM FOR MIXTURE',I4,' DOES NOT&
     & SUM T0 1., GROUP 1 VALUE ADJUSTED BY',F10.6)                     
      END SUBROUTINE CHIC                           
!*                                                                      
!XSREAD           ************                                          
!----EZM          *  XSREAD  *         AUG 1989                         
!                 ************                                          
!                                                                       
!     XSREAD READS THE XSECTS FROM THE SAME FILE AS OTHER INPUT DATA    
!                                                                       
      SUBROUTINE XSREAD (XSEC, NF, NG, NG4) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION XSEC (NG, NG4) 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      READ (LUNFIL(3), *, END = 99, ERR = 99) ( (XSEC (I, J), J = 1, NG4),      &
      I = 1, NG)                                                        
      RETURN 
   99 NF = 9999 
      RETURN 
      END SUBROUTINE XSREAD                         
!*                                                                      
!INMIX            ************                                          
!----EZM          *  INMIX   *         JUNE 1989                        
!                 ************                                          
!                                                                       
!     INMIX READS THE XSECTS FROM INPUT BINARY FILE (STANDARD OR LABAN) 
!     NOTE: THE INPUT BINARY FILE MUST BE SINGLE PRECISION !!!!         
!     ----                                                              
!     ON A STANDARD FILE:                                               
!     RECORD 1: NGB                                                     
!     RECORD 2: CHI                                                     
!     RECORD 3: DUMMY (GROUP STRUCTURE)                                 
!     RECORDS FOR EACH MIXTURE:                                         
!     ID, DESCRIPTION                                                   
!     XSECTS ARE GIVEN GROUPWISE FOR D,SIGA,NUSIGF,FIS,DISCON,AVFLUX,   
!     SCATT(ROWWISE)                                                    
!                                                                       
!     ON A LABAN FILE:                                                  
!     XSECTS ARE GIVEN IN A STRING OF D,SIGA,NUSIGF,CHI,SCATT(COLUMN    
!     AFTER COLUMN), FIS, DISCON, DMAT (COLUMN AFTER COLUMN)            
!                                                                       
      SUBROUTINE INMIX (S, XSEC, IDENT, NF, NG, NGG4) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      REAL S (NGG4) 
      DIMENSION XSEC (NGG4) 
      COMMON / DINV / IDIAG 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / VERSI / IVERS 
      CHARACTER LTITLE * 32 
      COMMON / XSFIL / LTITLE 
!                                                                       
      NFF = IABS (NF) 
!                                                                       
      IF (IVERS.EQ.1) THEN 
         OPEN (NFF, FILE = LTITLE, STATUS = 'OLD', FORM = 'UNFORMATTED',&
         ERR = 999)                                                     
!LPI       OPEN(NFF,FILE=LTITLE,STATUS='OLD',FORM='UNFORMATTED',        
!LPI &          ORGANIZATION='DYNAMIC',ERR=999)                         
      ELSE 
         REWIND NFF 
      ENDIF 
!                                                                       
      IF (NF.GT.0) THEN 
!                                                                       
!  READ CROSS SECTIONS FROM STANDARD FILE                               
!                                                                       
         READ (NFF, END = 99, ERR = 99) LG 
         IF (IDIAG.NE.0) THEN 
            WRITE (LUNFIL(2), 1000) IDIAG 
            STOP 'ONLY DIAGONAL D WITH STANDARD FILE' 
         ENDIF 
         IADD = 1 
         IADABS = IADD+LG 
         IADNUF = IADABS + LG 
         IADCHI = IADNUF + LG 
         IADSCM = IADCHI + LG 
         IADFIS = IADSCM + LG * LG 
         IADDIS = IADFIS + LG 
         READ (NFF, END = 99, ERR = 99) (S (IADCHI + I - 1), I = 1, LG) 
         READ (NFF) 
    1    READ (NFF, END = 99, ERR = 99) LDENT 
         IF (LDENT.EQ.IDENT.AND.LG.EQ.NG) GOTO 2 
!                                                                       
!        SKIP RECORDS (DESIRED CROSS SECTION WAS NOT FOUND - TRY AGAIN) 
!                                                                       
         DO 111 I = 1, LG 
            READ (NFF, END = 99, ERR = 99) 
            READ (NFF, END = 99, ERR = 99) 
  111    END DO 
         GOTO 1 
    2    DO 211 I = 1, LG 
            READ (NFF, END = 99, ERR = 99) S (IADD+I - 1), S (IADABS +  &
            I - 1), S (IADNUF + I - 1), S (IADFIS + I - 1), S (IADDIS + &
            I - 1)                                                      
            READ (NFF, END = 99, ERR = 99) (S (IADSCM + (I - 1) * LG +  &
            J - 1), J = 1, LG)                                          
  211    END DO 
!                                                                       
      ELSE 
!                                                                       
!  READ CROSS SECTIONS FROM LABAN FILE                                  
!                                                                       
    3    READ (NFF, END = 99, ERR = 99) LDENT, LG, JDIAG 
         IF (LDENT.EQ.IDENT.AND.LG.EQ.NG) GOTO 4 
!                                                                       
!        SKIP RECORDS (DESIRED CROSS SECTION WAS NOT FOUND - TRY AGAIN) 
!                                                                       
         READ (NFF, END = 99, ERR = 99) 
         GOTO 3 
    4    IF (JDIAG.NE.IDIAG) THEN 
            WRITE (LUNFIL(2), 2000) IDIAG, JDIAG 
            STOP 'ERROR WITH D MATRIX SPECIFICATION' 
         ENDIF 
         READ (NFF, END = 99, ERR = 99) S 
      ENDIF 
!                                                                       
!  CONVERT TO DOUBLE PRECISION                                          
!                                                                       
      DO 5 I = 1, NGG4 
         XSEC (I) = DBLE (S (I) ) 
    5 END DO 
      IF (IVERS.EQ.1) CLOSE (NFF) 
      RETURN 
!                                                                       
   99 NF = 9999 
      IF (IVERS.EQ.1) CLOSE (NFF) 
      RETURN 
  999 WRITE (LUNFIL(2), 3000) NFF, LTITLE 
      STOP 'XSEC FILE NOT FOUND' 
 1000 FORMAT (//1X,100('*')/' ONLY DIAGONAL D MATRIX AVAILABLE ON STANDA RD FILE, IDIAG ON CARD 2 SHOULD = 0')                             
 2000 FORMAT (//1X, 100 ('*')  / ' IDIAG =', I2, ' ON INPUT CARD 2 WHILE IDIAG =', I2, ' ON LABAN FILE.  THEY SHOULD BE THE SAME!')               
 3000 FORMAT(' FILE NO ',I2,'   FILE NAME ',A32,' NOT FOUND') 
      END SUBROUTINE INMIX                          
!MAPCOR                      26 MAI 1982         INPUTA              LAB
!                                                                       
!                            **********                                 
!                            * MAPCOR *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE MAPCOR (MAP, MAPIN, NDX, NDY, MX2, MNX2) 
!                                                                       
!  FIND MIXTURE ALLOCATION MAP FOR SUBDIVIDED NODES.                    
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP (MX2, * ), MAPIN (MNX2, * ), NDX ( * ), NDY ( * ) 
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (03), NX2), (WW (04)&
      , NY2), (WW (05), NNX), (WW (06), NNY)                            
!                                                                       
      DO 10 I = 1, NY2 
         DO 10 J = 1, NX2 
   10 MAP (J, I) = 0.D0 
!                                                                       
      JX = 1 
      DO 20 J = 1, NNX 
         DO 18 JJ = 1, NDX (J) 
            MAP (JX + JJ, 1) = MAPIN (J + 1, 1) 
   18    MAP (JX + JJ, NY2) = MAPIN (J + 1, NNY + 2) 
   20 JX = JX + NDX (J) 
      IY = 1 
      DO 30 I = 1, NNY 
         DO 28 II = 1, NDY (I) 
            MAP (1, IY + II) = MAPIN (1, I + 1) 
   28    MAP (NX2, IY + II) = MAPIN (NNX + 2, I + 1) 
   30 IY = IY + NDY (I) 
!                                                                       
      IY = 1 
      DO 40 I = 1, NNY 
         JX = 1 
         DO 38 J = 1, NNX 
            DO 36 II = 1, NDY (I) 
               DO 36 JJ = 1, NDX (J) 
   36       MAP (JX + JJ, IY + II) = MAPIN (J + 1, I + 1) 
   38    JX = JX + NDX (J) 
   40 IY = IY + NDY (I) 
!                                                                       
      DO 50 I = 2, NY + 1 
         DO 50 J = 2, NX + 1 
            IF (MAP (J, I) ) 48, 50, 50 
   48       IF (MAP (J + 1, I) .LE.0.AND.MAP (J - 1, I) .LE.0.AND.MAP ( &
            J, I + 1) .LE.0.AND.MAP (J, I - 1) .LE.0) MAP (J, I)        &
            = 0                                                         
   50 CONTINUE 
      RETURN 
      END SUBROUTINE MAPCOR                         
!ALBINL                      4  DEC 1983         INPUT               LAB
!                                                                       
!                            **********                                 
!                            * ALBINL *                                 
!                            **********                                 
!                                                                       
!      -----------------------------------                              
      SUBROUTINE ALBINL (ALB0, ALB, MG, MGL) 
!      -----------------------------------                              
!                                                                       
!         READING OF ALBEDO MATRICES.                                   
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      CHARACTER IBC (3) * 5, IIBC * 5 
      DIMENSION ALB0 (MG, MG, * ), ALB (MGL, MGL, * ) 
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      EQUIVALENCE (WW (08), NG), (WW (11), NL), (WW (15), NGL), (WW (28)&
      , MIXALB), (WW (62), MCASE)                                       
!                                                                       
      DATA IBC / 'BLACK', 'WHITE', 'ALBED' / 
!                                                                       
      IF (MCASE.NE.1) GOTO 70 
      WRITE (LUNFIL(2), 100) 
      DO 50 K = 1, MIXALB 
         READ (LUNFIL(3), 102) IIBC 
         DO 12 M = 1, 3 
            L = M 
            IF (IIBC.EQ.IBC (M) ) GOTO 14 
   12    END DO 
         WRITE (LUNFIL(2), 104) IIBC 
         STOP 'ERROR WHEN READING ALBEDO MATRICES' 
!                                                                       
   14    GOTO (20, 20, 30), L 
!                                                                       
!         BLACK AND WHITE BOUNDARIES.                                   
!                                                                       
   20    XALB = 0.D0 
         IF (L.EQ.2) XALB = 1.D0 
         DO 24 I = 1, NG 
            DO 22 J = 1, NG 
   22       ALB0 (I, J, K) = 0.D0 
   24    ALB0 (I, I, K) = XALB 
         GOTO 40 
!                                                                       
!         GROUP TO GROUP ALBEDO MATRICES (ONLY P-0 MOMENTS).            
!         READ ROW-WISE.                                                
!                                                                       
   30    READ (LUNFIL(3), * ) ( (ALB0 (I, J, K), J = 1, NG), I = 1, NG) 
!                                                                       
!         PRINT P-0 MOMENTS OF ALBEDO MATRIX                            
!                                                                       
   40    WRITE (LUNFIL(2), 114) 
         DO 42 I = 1, NG 
   42    WRITE (LUNFIL(2), 114) K, (ALB0 (I, J, K), J = 1, NG) 
   50 END DO 
!                                                                       
   70 DO 90 K = 1, MIXALB 
         DO 72 I = 1, NGL 
            DO 72 J = 1, NGL 
   72    ALB (I, J, K) = 0.D0 
         DO 76 I = 1, NGL, NL 
            II = (I - 1) / NL + 1 
            DO 76 J = 1, NGL, NL 
               JJ = (J - 1) / NL + 1 
               XALB = ALB0 (II, JJ, K) 
               DO 74 M = 1, NL 
                  ALB (J + M - 1, I + M - 1, K) = XALB 
   74          XALB = - XALB 
   76    CONTINUE 
   90 END DO 
      RETURN 
!                                                                       
  100 FORMAT  (//' ALBEDO MATRICES'/                                    &
     &        1X,15('*')/'   MIX')                                      
  102 FORMAT  (A5) 
  104 FORMAT  (//' ERROR IN READING ALBEDO MATRICES',3X,A5) 
  114 FORMAT  (1X,I4,3X,12F10.6,/(8X,12F10.6)) 
!                                                                       
      END SUBROUTINE ALBINL                         
!MAPALL                      26 MAI 1982         INPUT               LAB
!                                                                       
!                            **********                                 
!                            * MAPALL *                                 
!                            **********                                 
!                                                                       
!      ----------------------------------                               
      SUBROUTINE MAPALL (MAP, HX, HY, NXY2) 
!      ----------------------------------                               
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP (NXY2, 7), HX ( * ), HY ( * ) 
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (03), NX2), (WW (04)&
      , NY2), (WW (08), NG), (WW (11), NL), (WW (15), NGL), (WW (27),   &
      MIXCAL), (WW (29), MIXTOT), (WW (30), KURTOT), (WW (31), LUXTOT), &
      (WW (62), MCASE)                                                  
!                                                                       
!         KURTOT = HALF THE NR OF CURRENT EXPANSION COEFF (NR OF UNKNOWN
!         LUXTOT = NR OF NODES.                                         
!         MIXTOT = NR OF DIFFERENT SETS OF RESPONSE MATRICES NEEDED.    
!                                                                       
!         THE ARRAY MAP CONTAINS INFORMATION AS FOLLOWS:                
!         MAP(I,1) : CROSS SECTION MIXTURE FOR CELL I.                  
!         MAP(I,2) : ADRS OF Y-DIRECTION RESP MATRICES FOR CELL I.      
!         MAP(I,3) : ADRS OF X-DIRECTION RESP MATRICES FOR CELL I.      
!         MAP(I,4) : ADRS OF INGOING CURRENT  AT CELL I'S NORTH SIDE.   
!         MAP(I,5) : ADRS OF INGOING CURRENT  AT CELL I'S EAST  SIDE.   
!         MAP(I,6) : ADRS OF INGOING CURRENT  AT CELL I'S SOUTH SIDE.   
!         MAP(I,7) : ADRS OF INGOING CURRENT  AT CELL I'S WEST  SIDE.   
!                                                                       
      EPSH = 1.D-6 
      NGLB = NGL * NGL 
      IF (MCASE.NE.2) MIXTOT = 0 
      LUXTOT = 0 
      NOD = 0 
      DO 10 I = 1, NY2 
         DO 10 J = 1, NX2 
            NOD = NOD+1 
            MIX = MAP (NOD, 1) 
            IF (MIX) 2, 10, 4 
    2       MAP (NOD, 2) = ( - MAP (NOD, 1) - 1) * NGLB 
            GOTO 10 
    4       LUXTOT = LUXTOT + 1 
            IF (MCASE.EQ.2) GOTO 10 
!                                                                       
            L = NX2 + 1 
            HY1 = HY (I - 1) 
            HX1 = HX (J - 1) 
    6       L = L + 1 
            IF (L.EQ.NOD) GOTO 8 
            IF (MAP (L, 1) .NE.MAP (NOD, 1) ) GOTO 6 
            IL = L / NX2 
            JL = L - IL * NX2 - 1 
            IF (DABS (HX (JL) - HX1) .GT.EPSH.OR.DABS (HY (IL) - HY1)   &
            .GT.EPSH) GOTO 6                                            
            MAP (NOD, 2) = MAP (L, 2) 
            MAP (NOD, 3) = MAP (L, 3) 
            GOTO 10 
    8       MAP (NOD, 2) = MIXTOT 
            IF (DABS (HY1 - HX1) .GT.EPSH) MIXTOT = MIXTOT + 1 
            MAP (NOD, 3) = MIXTOT 
            MIXTOT = MIXTOT + 1 
   10 CONTINUE 
!                                                                       
      DO 11 LSIDE = 4, 7 
         DO 11 NOD = 1, NXY2 
   11 MAP (NOD, LSIDE) = - 1 
!                                                                       
      KURTOT = 0 
      NOD = 0 
      ISW = 1 
      DO 20 I1 = 1, NY2 
         JSW = ISW 
         DO 18 J1 = 1, NX2 
            NOD = NOD+1 
            IF (JSW) 18, 18, 12 
   12       IF (I1.EQ.1) GOTO 13 
            N1 = NOD-NX2 
            IF (MAP (NOD, 1) .LE.0.AND.MAP (N1, 1) .LE.0) GOTO 13 
            MAP (NOD, 4) = KURTOT 
            MAP (N1, 6) = KURTOT 
            KURTOT = KURTOT + NGL 
   13       IF (J1.EQ.NX2) GOTO 14 
            N2 = NOD+1 
            IF (MAP (NOD, 1) .LE.0.AND.MAP (N2, 1) .LE.0) GOTO 14 
            MAP (NOD, 5) = KURTOT 
            MAP (N2, 7) = KURTOT 
            KURTOT = KURTOT + NGL 
   14       IF (I1.EQ.NY2) GOTO 15 
            N3 = NOD+NX2 
            IF (MAP (NOD, 1) .LE.0.AND.MAP (N3, 1) .LE.0) GOTO 15 
            MAP (NOD, 6) = KURTOT 
            MAP (N3, 4) = KURTOT 
            KURTOT = KURTOT + NGL 
   15       IF (J1.EQ.1) GOTO 18 
            N4 = NOD-1 
            IF (MAP (NOD, 1) .LE.0.AND.MAP (N4, 1) .LE.0) GOTO 18 
            MAP (NOD, 7) = KURTOT 
            MAP (N4, 5) = KURTOT 
            KURTOT = KURTOT + NGL 
!                                                                       
   18    JSW = - JSW 
   20 ISW = - ISW 
!                                                                       
      RETURN 
      END SUBROUTINE MAPALL                         
!PERESP                      4  DEC 1983         INPUT               LAB
!                                                                       
!                            **********                                 
!                            * RERESP *                                 
!                            **********                                 
!                                                                       
!      ---------------------------------------------                    
      SUBROUTINE RERESP (NXY2, MAP, ONE, A, T, B, F, HELP) 
!      ---------------------------------------------                    
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP (NXY2, 7), ONE ( * ), A ( * ), T ( * ), B ( * ),    &
      F ( * ), HELP ( * )                                               
      CHARACTER ITITLE * 60 
!                                                                       
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      EQUIVALENCE (WW (08), NG), (WW (11), NL), (WW (18), NGLB),        &
      (WW (19), NGLA), (WW (20), NGLF), (WW (27), MIXCAL), (WW (73),    &
      MIXXSE), (WW (74), MIXEXT)                                        
!                                                                       
!      -----------------------                                          
      CALL ONECAL (NG, NL, ONE) 
!      -----------------------                                          
      WRITE (LUNFIL(2), 2000) 
!                                                                       
!  SWEEP OVER ALL SETS OF RESPONSE MATRIX MIXTURES.                     
!                                                                       
      DO 100 KK = 1, MIXEXT 
!                                                                       
         READ (LUNFIL(3), 1000, ERR = 5) MIX, L0, NANT, ITITLE 
         WRITE (LUNFIL(2), 2010) MIX, L0, NANT, ITITLE 
         L1 = L0 + 1 
         NLX = MAX (NL, L1) 
!                                                                       
         IF (MIX.LE.MIXCAL.OR.MIX.GT.MIXXSE) GOTO 5 
         IF (NANT.LT.1.OR.NANT.GT.2) GOTO 5 
         GOTO 8 
!                                                                       
!  ERROR STOP                                                           
!                                                                       
    5    STOP 'ERRONEOUS RESP MATRIX MIXTURE' 
!                                                                       
!  FIND NODE WHERE THIS MIXTURE TYPE IS USED                            
!                                                                       
    8    DO 10 N = 1, NXY2 
            NOD = N 
   10    IF (MIX.EQ.MAP (N, 1) ) GOTO 20 
         NGANG = NANT 
         GOTO 40 
!                                                                       
   20    MY = MAP (NOD, 2) 
         MX = MAP (NOD, 3) 
         MYA = MY * NGLA + 1 
         MXA = MX * NGLA + 1 
         MYB = MY * NGLB + 1 
         MXB = MX * NGLB + 1 
         MYF = MY * NGLF + 1 
         MXF = MX * NGLF + 1 
!                                                                       
!  READ RESPONSE MATRICES A-Y, T-Y, B-Y, F-Y.                           
!                                                                       
!      --------------------------------------------------------         
         CALL INRESP (1, 'REFLECT.', NG, NL, L1, NLX, ONE, HELP, A (MYA)&
         )                                                              
         CALL INRESP (2, 'TRANSM. ', NG, NL, L1, NLX, ONE, HELP, T (MYA)&
         )                                                              
         CALL INRESP (3, 'SIDE-TR.', NG, NL, L1, NLX, ONE, HELP, B (MYB)&
         )                                                              
         CALL INRESP (4, 'TOT FLUX', NG, NL, L1, NLX, ONE, HELP, F (MYF)&
         )                                                              
!      --------------------------------------------------------         
!                                                                       
         N0 = 1 
         IF (MX.NE.MY) N0 = 2 
         NGANG = 1 
!                                                                       
!  THERE ARE NANT (=1 OR 2) SETS OF RESP MATRICES STORED FOR THIS MIXTUR
!  THERE SHOULD BE N0 SETS.  IF NANT.NE.N0  THEN TAKE SPECIAL ACTIONS.  
!                                                                       
         IF (N0.EQ.1.AND.NANT.EQ.1) GOTO 100 
         IF (N0.EQ.2.AND.NANT.EQ.1) GOTO 30 
         IF (N0.EQ.1.AND.NANT.EQ.2) GOTO 40 
!                                                                       
!  READ RESPONSE MATRICES A-X, T-X, B-X, F-X.                           
!                                                                       
!      --------------------------------------------------------         
         CALL INRESP (1, 'REFLECT.', NG, NL, L1, NLX, ONE, HELP, A (MXA)&
         )                                                              
         CALL INRESP (2, 'TRANSM. ', NG, NL, L1, NLX, ONE, HELP, T (MXA)&
         )                                                              
         CALL INRESP (3, 'SIDE-TR.', NG, NL, L1, NLX, ONE, HELP, B (MXB)&
         )                                                              
         CALL INRESP (4, 'TOT FLUX', NG, NL, L1, NLX, ONE, HELP, F (MXF)&
         )                                                              
!      --------------------------------------------------------         
         GOTO 100 
!                                                                       
!  SET BY FORCE  A-X EQUAL TO A-Y, ETC, SINCE A-X IS MISSING IN INDATA. 
!                                                                       
!      -------------------------------                                  
   30    CALL FMOVT (A (MXA), A (MYA), NGLA) 
         CALL FMOVT (T (MXA), T (MYA), NGLA) 
         CALL FMOVT (B (MXB), B (MYB), NGLB) 
         CALL FMOVT (F (MXF), F (MYF), NGLF) 
!      -------------------------------                                  
         GOTO 100 
!                                                                       
!  READ RESPONSE MATRIX SET BUT DO NOT STORE RESULTS.                   
!                                                                       
   40    DO 50 M = 1, NGANG 
            DO 50 N = 1, 4 
               NN = NG * L1 
               IF (N.EQ.4) NN = NG 
!      READ(5,1010)                                                     
               DO 50 I = 1, NN 
   50    READ (LUNFIL(3), * ) 
!                                                                       
  100 END DO 
!                                                                       
 1000 FORMAT  (3I5, 5X, A60) 
 1010 FORMAT  (8F10.0) 
 2000 FORMAT  (1H1,'EXTERNAL RESPONSE MATRICES'/1X,26('*')) 
 2010 FORMAT  (//' MIX =',I4,4X,'L =',I2,4X,'SETS =',I2,4X, A60/        &
     &        1X,30(1H-))                                               
!                                                                       
      RETURN 
      END SUBROUTINE RERESP                         
!INRESP                      1  NOV 1982         INPUT               LAB
!                                                                       
!                            **********                                 
!                            * INRESP *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE INRESP (IVAL, TEXT, NG, NL, L1, NLX, ONE, R, A) 
!                                                                       
!* FULL MATRIX IS READ INTO R                             3 MAJ 1985    
!* CONDENSED MATRIX (WITHOUT ZEROS) IS STORED INTO A      3 MAJ 1985    
!* **************************************************                   
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      CHARACTER TEXT * 8 
      DIMENSION ONE (NL, NG, NL, * ), R (NLX, NG, NLX, * ), A ( * ) 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      WRITE (LUNFIL(2), 2000) TEXT 
      N = 0 
      NLI = NL 
      IF (IVAL.EQ.4) NLI = 1 
!      READ(5,1000)                                                     
!                                                                       
!  READ MATRIX                                                          
!                                                                       
      DO 70 IG = 1, NG 
         EM = 1.D0 
         DO 60 IL = 1, NLI 
            J0 = 1 
            IF (IL.GT.L1) GOTO 10 
            READ (LUNFIL(3), * ) ( (R (IL, IG, JL, JG), JL = 1, L1), JG = 1, NG) 
!      READ(5,1000)  ((R(IL,IG,JL,JG), JL=1,L1), JG=1,NG)               
!                                                                       
!  NL = ORDER OF APPROXIMATION TO BE USED LATER IN ITERATIONS.          
!  L1 = ACTUAL ORDER OF APPROXIMATION FOR MATRICES JUST FOUND.          
!  IF  L1 > NL  THEN JUST IGNORE THE SUPERFLUOUS MATRIX ELEMENTS.       
!  IF  L1 < NL  THEN FILL MISSING POSITIONS WITH ZEROS. FOR MATRICES A A
!               SET DIAGONAL EQUAL TO ZEROTH ORDER RESPONSE.            
!                                                                       
            J0 = L1 + 1 
   10       DO 50 JG = 1, NG 
               IF (J0.GT.NL) GOTO 30 
               DO 20 JL = J0, NL 
                  R (IL, IG, JL, JG) = 0.D0 
                  IF (IL.NE.JL) GOTO 20 
                  IF (IVAL.EQ.1) R (IL, IG, IL, JG) = R (1, IG, 1, JG)  &
                  * EM                                                  
                  IF (IVAL.EQ.2) R (IL, IG, IL, JG) = R (1, IG, 1, JG) 
   20          END DO 
!                                                                       
!  RESTORE THE MATRIX FOUND IN ARRAY A().                               
!                                                                       
   30          DO 40 JL = 1, NL 
                  ON = ONE (IL, IG, JL, JG) 
                  IF (ON.LT.0.D0.AND. (IVAL.LE.2.OR.IVAL.EQ.4) ) GOTO   &
                  40                                                    
                  N = N + 1 
                  A (N) = R (IL, IG, JL, JG) 
   40          END DO 
   50       END DO 
            WRITE (LUNFIL(2), 2010) ( (R (IL, IG, JL, JG), JL = 1, NL), JG = 1, &
            NG)                                                         
   60    EM = - EM 
!                                                                       
!      SHOULD LAST ROWS OF RESPONSE MATRIX BE SKIPPED ?                 
         IF (IVAL.EQ.4.OR.NLI.GE.L1) GOTO 70 
         DO 65 IL = NLI + 1, L1 
   65    READ (LUNFIL(3), * ) 
   70 END DO 
!                                                                       
 1000 FORMAT  (8F10.0) 
 2000 FORMAT  (/1X,A8,' MATRIX') 
 2010 FORMAT  (1X,10F10.5) 
!                                                                       
      RETURN 
      END SUBROUTINE INRESP                         
!ITERA                       1  NOV 1982         ITERA               LAB
!                                                                       
!                            ***************                            
!                            * ITERA/LABAN *                            
!                            ***************                            
!                                                                       
!      ITERATION LINK FOR CALCULATION OF PARTIAL CURRENTS.              
!      --------------------------------------------------               
!                                                                       
      SUBROUTINE ITERA (C, IEXIT, HELL) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION LAMBDA 
      DIMENSION C ( * ), HELL ( * ) 
!                                                                       
      COMMON / CLAB / WW (100) 
      COMMON / ADRS / AD (50) 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      EQUIVALENCE (AD (05), LCHI), (AD (06), LXSEC), (AD (07), LALB),   &
      (AD (08), LONE), (AD (09), LMAP), (AD (10), LCHECK), (AD (15),    &
      LF), (AD (16), LA), (AD (17), LT), (AD (18), LB), (AD (19),       &
      LAH), (AD (21), LOLD), (AD (22), LCUR), (AD (23), LFLUX), (AD (24)&
      , LPOW), (AD (50), LAST)                                          
!                                                                       
      EQUIVALENCE (WW (07), NXY2), (WW (08), NG), (WW (10), NG4),       &
      (WW (11), NL), (WW (15), NGL), (WW (30), KURTOT), (WW (31),       &
      LUXTOT), (WW (33), MAXOUT), (WW (34), ITER), (WW (35), ITIN),     &
      (WW (37), EPSK), (WW (42), XK), (WW (44), OMEGA), (WW (46),       &
      XK1), (WW (47), XK2), (WW (48), XK3), (WW (49), EV1), (WW (50),   &
      EV2), (WW (51), EV3), (WW (52), LAMBDA), (WW (53), EIG), (WW (54),&
      ITTOT), (WW (55), TIMIT), (WW (56), TIMER), (WW (62), MCASE),     &
      (WW (65), XKRMM), (WW (76), DELTA)                                
!                                                                       
!                                                                       
!         INITIAL GUESS OF CURRENTS.                                    
      IF (ITER.GT.1) GOTO 10 
!      -----------------------------                                    
      CALL INGELA (C (LCUR), C (LOLD) ) 
!      -----------------------------                                    
      IF (MCASE.EQ.1) LAMBDA = 1.0D0 
      ITTOT = 0 
      WRITE (LUNFIL(2), 31) 
!                                                                       
!         GO TO CONACC TO PREPARE FOR ITERATION.                        
!                                                                       
!      -----------------------------                                    
   10 CALL CONACC (C (LCUR), C (LOLD) ) 
      CALL ONECAL (NG, NL, C (LONE) ) 
!      -----------------------------                                    
!                                                                       
!         CALL POWIT FOR FINDING CURRENT SOLUTION BY POWER ITERATION.   
!                                                                       
!      -------------------------------------------------                
      CALL POWIT (NXY2, 0, C (LMAP), C (LONE), C (LA), C (LT), C (LB),  &
      C (LALB), C (LCUR), C (LOLD), HELL, NGL)                          
!      -------------------------------------------------                
      ITTOT = ITTOT + ITIN 
!                                                                       
!         FIND ESTIMATE OF K-EFF (XK) BY REGULA FALSI OR PARABOLA FIT.  
!                                                                       
      WRITE (LUNFIL(2), 32) ITER, ITIN, XK, LAMBDA, DELTA, OMEGA 
      EPS7 = 1.D-7 
      IF (ITER.EQ.1) GOTO 20 
      XK1 = 1.D0 / XK 
      EV1 = 1.D0 / LAMBDA 
      XK0 = XK1 
      E12 = EV1 - EV2 
      IF (DABS (E12) .LT.EPS7) GOTO 17 
      IF (ITER.GT.2) GOTO 14 
   12 XK0 = XK1 + (1.D0 - EV1) / E12 * (XK1 - XK2) 
      GOTO 15 
!                                                                       
   14 E23 = EV2 - EV3 
      E31 = EV3 - EV1 
      IF (DABS (E23) .LT.EPS7.OR.DABS (E31) .LT.EPS7) GOTO 12 
      XK0 = (1.D0 - EV2) / E23 * (XK2 - XK3) - (1.D0 + EV1 - EV2 - EV3) &
      / E12 * (XK1 - XK2)                                               
      XK0 = XK1 + (1.D0 - EV1) / E31 * XK0 
   15 IF (XK0.LE.0.D0) XK0 = 0.5D0 * XK1 
   17 XK = 1.D0 / XK0 
      XK3 = XK2 
      XK2 = XK1 
      EV3 = EV2 
      EV2 = EV1 
!                                                                       
!         CHECK CONVERGENCE OF OUTER ITERATION.                         
!                                                                       
      DELTA1 = DABS (1.D0 - LAMBDA) 
      DELTA2 = DABS (1.D0 - XK0 / XK1) 
      DELTA3 = MAX (DELTA1, DELTA2) 
      IF (DELTA3.LT.EPSK.AND.ITIN.LT.10) GOTO 20 
      IF (ITER.LT.MAXOUT) GOTO 60 
      WRITE (LUNFIL(2), 18) 
!                                                                       
!         CALL POWIT FOR CALCULATION OF CURRENTS ENTERING BLACK CELLS   
!         (CURRENTS IN WHITE CELLS ARE ALREADY KNOWN),                  
!         LEAKAG FOR CALCULATION OF THE LEAKAGE,                        
!         FLUXL FOR CALCULATION OF THE FLUX,                            
!         AND CHECKL FOR CHECK OF THE NEUTRON BALANCE.                  
!         IN CHECKL ALSO CALCULATE XKRMM, THE SECOND GUESS OF K-EFF.    
!                                                                       
!      -----------------------------------------------------------------
   20 CALL POWIT (NXY2, 1, C (LMAP), C (LONE), C (LA), C (LT), C (LB),  &
      C (LALB), C (LOLD), C (LOLD), HELL, NGL)                          
      CALL LEAKAG (NXY2, C (LMAP), C (LCHECK), C (LCUR), C (LOLD) ) 
      CALL FLUXL (NXY2, C (LMAP), C (LONE), C (LF), C (LCUR), C (LOLD), &
      C (LFLUX), HELL, NGL)                                             
      CALL CHECKL (NG, NG4, C (LMAP), C (LXSEC), C (LCHECK), C (LFLUX) ) 
!      -----------------------------------------------------------------
!      MEMT = LAST                                                      
      IF (ITER.GT.1) GOTO 23 
      XK2 = 1.D0 / XK 
      EV2 = 1.D0 / LAMBDA 
      XK = XKRMM 
      IF (MAXOUT.GT.1) GOTO 60 
      WRITE (LUNFIL(2), 29) LAMBDA, XKRMM 
      GOTO 24 
   23 WRITE (LUNFIL(2), 30) XK, XKRMM 
   24 WRITE (LUNFIL(2), 34) ITTOT, KURTOT 
!  24  WRITE(6,34) ITTOT,TIMER,TIMIT,MEMT,KURTOT                        
!                                                                       
      IEXIT = 1 
      GOTO 90 
   60 IEXIT = 2 
   90 RETURN 
!                                                                       
   18 FORMAT  (/1X,125(1H*)/' WARNING : NO K-EFF ITERATION CONVERGENCE' &
     &         /1X,125(1H*))                                            
   29 FORMAT  (////1X,22(1H*)/' * LAMBDA =',F10.6,' *       K-RMM =',   &
     &         F10.6/1X,22(1H*))                                        
   30 FORMAT  (////1X,21(1H*)/' * K-EFF =',F10.6,' *',/1X,21(1H*)//     &
     &              1X,'  K-EFF FROM BALANCE =',F10.6)                  
   31 FORMAT  (1X,'ITERATION PROCESS'/1X,17(1H*)//                      &
     &         1X,'ITER NO  IN-ITER    K-EFF     EIGEN     DELTA',      &
     &         '     OMEGA')                                            
   32 FORMAT  (1X,I6,I9,2X,5F10.6) 
   34 FORMAT  (///' TOTAL NR OF INNER ITERATIONS =',I5/                 &
     &            ' NUMBER OF UNKNOWNS (KURTOT)  =',I5)                 
!  34  FORMAT (///' TOTAL NR OF INNER ITER =',I5/' TIME RESPONSE MAT',  
!    &        6X,1H= ,F8.2,' SEC'/ ' TIME ITERATION',9X,                
!    &        1H= ,F8.2/' MEMORY',                                      
!    &        17X,1H= ,I6/' KURTOT',17X,1H= ,I5)                        
!                                                                       
      END SUBROUTINE ITERA                          
!INGELA                      26 MAI 1982         ITERA               LAB
!                                                                       
!                            ***************                            
!                            * INGELA/ITER *                            
!                            ***************                            
!                                                                       
      SUBROUTINE INGELA (CUR, X) 
!                                                                       
!         INITIAL GUESS OF THE CURRENT DISTRIBUTION.                    
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION CUR ( * ), X ( * ) 
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (07), NXY2),        &
      (WW (08), NG), (WW (11), NL), (WW (26), NLOLD), (WW (30), KURTOT),&
      (WW (31), LUXTOT), (WW (61), PNORM), (WW (62), MCASE), (WW (64),  &
      IFIL)                                                             
!                                                                       
      DO 10 I = 1, KURTOT 
   10 CUR (I) = 0.D0 
      IF (MCASE.EQ.2) GOTO 20 
      CC = 1D0 / DSQRT (DBLE (KURTOT) ) 
      DO 12 I = 1, KURTOT, NL 
   12 CUR (I) = CC 
      GOTO 30 
!                                                                       
!  READ START GUESSES OF CURRENTS FROM FILE IFIL                        
!                                                                       
   20 NLM = MIN (NL, NLOLD) 
      NSID = KURTOT / NL 
      IC = 0 
!                                                                       
      REWIND IFIL 
      READ (IFIL) 
      READ (IFIL) 
      READ (IFIL) 
      READ (IFIL) 
!                                                                       
      DO 24 N = 1, NSID 
         READ (IFIL) (X (I), I = 1, NLOLD) 
         DO 22 I = 1, NLM 
   22    CUR (I + IC) = X (I) / PNORM 
   24 IC = IC + NL 
   30 RETURN 
      END SUBROUTINE INGELA                         
!POWIT                       26 MAI 1982         ITERA               LAB
!                                                                       
!                            **************                             
!                            * POWIT/ITER *                             
!                            **************                             
!                                                                       
      SUBROUTINE POWIT (NXY2, NEXTRA, MAP, ONE, A, T, B, ALB, CUR, OLD, &
      HELL, NHH)                                                        
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION LAMBDA 
      DIMENSION MAP (NXY2, 7), ONE ( * ), A ( * ), T ( * ), B ( * ),    &
      ALB ( * ), CUR ( * ), OLD ( * ), HELL (NHH, 4)                    
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (03), NX2), (WW (04)&
      , NY2), (WW (08), NG), (WW (11), NL), (WW (15), NGL), (WW (18),   &
      NGLB), (WW (19), NGLA), (WW (27), MIXCAL), (WW (35), ITIN),       &
      (WW (39), ICONV), (WW (52), LAMBDA)                               
!                                                                       
!         POWIT FINDS THE SOLUTION BY POWER ITERATION OF 2-CYCLIC MATRIX
!                                                                       
      ICONV = 0 
!         FIND LOCATION OF CURRENTS IN ARRAY MAP.                       
      ITIN = 0 
    1 ISW0 = 1 
      ITIN = ITIN + 1 
    2 NOD = 0 
      ISW = ISW0 
      DO 50 I1 = 1, NY2 
         JSW = ISW 
         DO 40 J1 = 1, NX2 
            NOD = NOD+1 
            IF (JSW) 40, 40, 4 
    4       IF (MAP (NOD, 1) ) 30, 40, 8 
!                                                                       
!                                                                       
    8       KUR1 = MAP (NOD, 4) 
            KUR2 = MAP (NOD, 5) 
            KUR3 = MAP (NOD, 6) 
            KUR4 = MAP (NOD, 7) 
            DO 10 I = 1, NGL 
               HELL (I, 1) = CUR (I + KUR1) 
               HELL (I, 2) = CUR (I + KUR2) 
               HELL (I, 3) = CUR (I + KUR3) 
   10       HELL (I, 4) = CUR (I + KUR4) 
!                                                                       
!         THE VERY INNER LOOP OF LABANM.  CURRENTS ARE CALCULATED.      
!         HELL(NGL,1) ARE THE INGOING CURRENTS AT NORTH SIDE, ETC.      
!         CUR IS THE RESULTING OUTGOING CURRENT.                        
!         NOTICE THAT INGOING AND OUTGOING CURRENTS ARE STORED IN THE SA
!                                                                       
            MY = MAP (NOD, 2) 
            MX = MAP (NOD, 3) 
            MYA = MY * NGLA 
            MXA = MX * NGLA 
            MYB = MY * NGLB 
            MXB = MX * NGLB 
            IONE = 0 
            DO 20 I = 1, NGL 
               S1 = 0.D0 
               S2 = 0.D0 
               S3 = 0.D0 
               S4 = 0.D0 
               DO 18 J = 1, NGL 
                  MYB = MYB + 1 
                  MXB = MXB + 1 
                  IONE = IONE+1 
                  ON = ONE (IONE) 
                  BYMINS = ON * B (MYB) 
                  BXMINS = ON * B (MXB) 
                  S1 = S1 + B (MXB) * HELL (J, 2) + BXMINS * HELL (J, 4) 
                  S2 = S2 + B (MYB) * HELL (J, 3) + BYMINS * HELL (J, 1) 
                  S3 = S3 + B (MXB) * HELL (J, 4) + BXMINS * HELL (J, 2) 
                  S4 = S4 + B (MYB) * HELL (J, 1) + BYMINS * HELL (J, 3) 
                  IF (ON) 18, 18, 17 
   17             MYA = MYA + 1 
                  MXA = MXA + 1 
                  S1 = S1 + A (MYA) * HELL (J, 1) + T (MYA) * HELL (J,  &
                  3)                                                    
                  S2 = S2 + A (MXA) * HELL (J, 2) + T (MXA) * HELL (J,  &
                  4)                                                    
                  S3 = S3 + A (MYA) * HELL (J, 3) + T (MYA) * HELL (J,  &
                  1)                                                    
                  S4 = S4 + A (MXA) * HELL (J, 4) + T (MXA) * HELL (J,  &
                  2)                                                    
   18          END DO 
               CUR (I + KUR1) = S1 
               CUR (I + KUR2) = S2 
               CUR (I + KUR3) = S3 
   20       CUR (I + KUR4) = S4 
            GOTO 40 
!                                                                       
!         CALCULATION FOR FRINGE CELLS (BOUNDARY REFLECTION).           
!                                                                       
   30       DO 38 LSIDE = 4, 7 
               KUR = MAP (NOD, LSIDE) 
               IS = KUR + 1 
               IF (IS) 38, 38, 31 
   31          MA = MAP (NOD, 2) 
               DO 32 I = 1, NGL 
   32          HELL (I, 1) = CUR (I + KUR) 
               IE = KUR + NGL 
               DO 36 I = IS, IE 
                  S1 = 0.D0 
                  DO 34 J = 1, NGL 
                     MA = MA + 1 
   34             S1 = S1 + ALB (MA) * HELL (J, 1) 
   36          CUR (I) = S1 * LAMBDA 
   38       END DO 
!                                                                       
   40    JSW = - JSW 
   50 ISW = - ISW 
      IF (NEXTRA.EQ.1) RETURN 
      ISW0 = - ISW0 
      IF (ISW0.LT.0) GOTO 2 
!                                                                       
!         CALL CONAC1 FOR CONVERGENCE ACCELERATION AND TO CHECK         
!         IF THE ITERATION HAS CONVERGED.                               
!                                                                       
!      -----------------                                                
      CALL CONAC1 (CUR, OLD) 
!      -----------------                                                
      IF (ICONV.EQ.1) RETURN 
      GOTO 1 
      END SUBROUTINE POWIT                          
!CONACC                      1  NOV 1982         ITERA               LAB
!                                                                       
!                            ************************                   
!                            * CONACC/CONAC1 - ITER *                   
!                            ************************                   
!                                                                       
      SUBROUTINE CONACC (CUR, OLD) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION LAMBDA 
      DIMENSION CUR ( * ), OLD ( * ) 
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (11), NL), (WW (30), KURTOT), (WW (32), MAXIN),   &
      (WW (35), ITIN), (WW (38), EPSIT), (WW (39), ICONV), (WW (41),    &
      KACC), (WW (44), OMEGA), (WW (45), ALFA), (WW (52), LAMBDA),      &
      (WW (53), EIG), (WW (76), DELTA)                                  
      DATA EPSALP / 0.02D0 / 
!      DATA EPSALP /0.05D0/                                             
!                                                                       
      ITIN = 0 
      EIG2 = 0.D0 
      ALFA = 0.D0 
      OMEG = 1.D0 
      DO 2 I = 1, KURTOT 
    2 OLD (I) = CUR (I) 
      RETURN 
!                                                                       
!      ------------------                                               
      ENTRY CONAC1 (CUR, OLD) 
!      ------------------                                               
!                                                                       
!         EIGENVALUE IS CALCULATED.                                     
!         EMAX=UPPER BOUND ON EIGENVALUE                                
!         EMIN=LOWER BOUND ON EIGENVALUE                                
!         EIG = EIGENVALUE ESTIMATE IN L-2 NORM.                        
!                                                                       
      EIG = 0.D0 
      EMAX = 0.D0 
      EMIN = 1.0D+20 
      DO 6 I = 1, KURTOT, NL 
         IF (DABS (OLD (I) ) .LT.1.D-20) GOTO 6 
         XKX = DABS (CUR (I) / OLD (I) ) 
!      XKX=CUR(I)/OLD(I)                                                
         EMAX = MAX (EMAX, XKX) 
         EMIN = MIN (EMIN, XKX) 
    6 EIG = EIG + CUR (I) * CUR (I) 
      EIG = DSQRT (EIG) 
!                                                                       
!         CONVERGENCE ACCELERATION BY SOR METHOD.  RENORMALIZATION.     
!                                                                       
      XX = 1.D0 / EIG 
      DO 8 I = 1, KURTOT 
    8 CUR (I) = OLD (I) + OMEG * (XX * CUR (I) - OLD (I) ) 
!                                                                       
!         ESTIMATE THE EIGENVALUE BY GEOMETRICAL EXTRAPOLATION.         
!                                                                       
      EPS7 = 1.D-7 
      IF (ITIN.LT.4) GOTO 18 
      XX = 2.D0 * EIG2 - EIG - EIG3 
      IF (DABS (XX) .LT.EPS7) GOTO 18 
      ALFA = (EIG2 - EIG3) / XX 
      IF (DABS (ALFA) .LT.EPS7) GOTO 18 
      EIGAS = EIG2 + ALFA * (EIG - EIG2) 
      IF (DABS (1.D0 - ALFOLD / ALFA) .GE.EPSALP) GOTO 18 
      IF (KACC) 18, 18, 10 
!                                                                       
!         CONVERGENCE ACCELERATION BY GEOMETRICAL EXTRAPOLATION         
!         (ONLY IF SOLUTION IN EXPONENTIAL MODE).                       
!                                                                       
   10 DO 12 I = 1, KURTOT 
   12 CUR (I) = OLD (I) + ALFA * (CUR (I) - OLD (I) ) 
      EIG = EIGAS 
   18 EIG3 = EIG2 
      EIG2 = EIG 
      ALFOLD = ALFA 
      OMEG = OMEGA 
      LAMBDA = DSQRT (EIG) 
      DO 20 I = 1, KURTOT 
   20 OLD (I) = CUR (I) 
!                                                                       
!         CHECK IF ITERATION HAS CONVERGED.                             
!                                                                       
      DELTA = DABS (1.D0 - EMIN / EMAX) 
      IF (DELTA.LT.EPSIT) ICONV = 1 
      IF (ITIN.EQ.1) ICONV = 0 
      IF (ITIN.GE.MAXIN) ICONV = 1 
!                                                                       
      RETURN 
      END SUBROUTINE CONACC                         
!FLUXL                       26 MAI 1982         ITERA               LAB
!                                                                       
!                            *********                                  
!                            * FLUXL *                                  
!                            *********                                  
!                                                                       
      SUBROUTINE FLUXL (NXY2, MAP, ONE, F, CURA, CURB, FLUX, X, NHH) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP (NXY2, 7), ONE ( * ), F ( * ), CURA ( * ), CURB ( * &
      ), FLUX ( * ), X (NHH, 4)                                         
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (03), NX2), (WW (04), NY2), (WW (08), NG),        &
      (WW (11), NL), (WW (15), NGL), (WW (20), NGLF)                    
!                                                                       
!         FLUXL CALCULATES THE TOTAL FLUX IN EACH GROUP AND CELL.       
!                                                                       
!         NODE SYSTEM IS DIVEDED INTO A CHECKERBOARD PATTERN :          
!           + - + - + - + - + .....                                     
!           - + - + - + - + - .....                                     
!           + - + - + - + - + .....                                     
!           ......                                                      
!         CURA()  IS  IN-CURRENT IN  "+"  NODES.                        
!         CURB()  IS  IN-CURRENT IN  "-"  NODES.                        
!                                                                       
      KFLUX = 0 
      NOD = 0 
      ISW = 1 
      DO 24 I1 = 1, NY2 
         JSW = ISW 
         DO 22 J1 = 1, NX2 
            NOD = NOD+1 
            IF (MAP (NOD, 1) ) 22, 22, 2 
    2       DO 4 LSIDE = 1, 4 
               KUR = MAP (NOD, LSIDE+3) 
               DO 4 I = 1, NGL 
                  XX = CURA (I + KUR) 
                  IF (JSW.LT.0) XX = CURB (I + KUR) 
    4       X (I, LSIDE) = XX 
!                                                                       
            MFY = MAP (NOD, 2) * NGLF 
            MFX = MAP (NOD, 3) * NGLF 
            DO 18 I = 1, NG 
               S1 = 0.D0 
               DO 16 J = 1, NGL 
                  IF (ONE (J) ) 16, 16, 15 
   15             MFY = MFY + 1 
                  MFX = MFX + 1 
                  S1 = S1 + F (MFY) * (X (J, 1) + X (J, 3) ) + F (MFX)  &
                  * (X (J, 2) + X (J, 4) )                              
   16          END DO 
               KFLUX = KFLUX + 1 
   18       FLUX (KFLUX) = S1 
!                                                                       
   22    JSW = - JSW 
   24 ISW = - ISW 
!                                                                       
      RETURN 
      END SUBROUTINE FLUXL                          
!LEAKAG                      26 MAI 1982         ITERA               LAB
!                                                                       
!                            **********                                 
!                            * LEAKAG *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE LEAKAG (NXY2, MAP, XLEAK, CURA, CURB) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION LAMBDA 
      DIMENSION MAP (NXY2, 7), XLEAK ( * ), CURA ( * ), CURB ( * ) 
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (03), NX2), (WW (04), NY2), (WW (08), NG),        &
      (WW (11), NL), (WW (15), NGL), (WW (30), KURTOT), (WW (52),       &
      LAMBDA)                                                           
!                                                                       
!         MULTIPLY CURRENTS IN BLACK NODES WITH THE EIGENVALUE.         
!                                                                       
      XX = 1.D0 / LAMBDA 
      DO 10 I = 1, KURTOT 
   10 CURB (I) = XX * CURB (I) 
!                                                                       
!         CALCULATION OF LEAKAGE IN EACH GROUP (XLEAK) AND TOTAL LEAKAGE
!                                                                       
      NN = 6 * NG 
      DO 2 I = 1, NN 
    2 XLEAK (I) = 0.D0 
      ESC = 0.D0 
      NOD = 0 
      ISW = 1 
      DO 30 I1 = 1, NY2 
         JSW = ISW 
         DO 28 J1 = 1, NX2 
            NOD = NOD+1 
            IF (MAP (NOD, 1) ) 15, 28, 15 
   15       DO 27 LSIDE = 4, 7 
               JS = MAP (NOD, LSIDE) + 1 
               IF (MAP (NOD, 1) ) 21, 16, 16 
   16          IF (JSW) 17, 17, 19 
   17          DO 18 I = 1, NG 
                  XLEAK (I) = XLEAK (I) + CURB (JS) - LAMBDA * CURA (JS) 
   18          JS = JS + NL 
               GOTO 27 
   19          DO 20 I = 1, NG 
                  XLEAK (I) = XLEAK (I) + CURA (JS) - LAMBDA * CURB (JS) 
   20          JS = JS + NL 
               GOTO 27 
   21          IF (JS) 27, 27, 22 
   22          DO 23 I = 1, NG 
                  ESC = ESC + DBLE (JSW) * (CURA (JS) - CURB (JS) ) 
   23          JS = JS + NL 
   27       END DO 
   28    JSW = - JSW 
   30 ISW = - ISW 
      NG1 = NG + 1 
      XLEAK (NG1) = ESC 
!                                                                       
      RETURN 
      END SUBROUTINE LEAKAG                         
!CHECKL                      31 OCT 1983         ITERA               LAB
!                                                                       
!                            **********                                 
!                            * CHECKL *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE CHECKL (LG, NG4, MAP, XSEC, CHECK, FLUX) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP ( * ), XSEC (LG, NG4, * ), CHECK (LG, 6), FLUX (LG, &
      * )                                                               
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (08), NG), (WW (11),&
      NL), (WW (15), NGL), (WW (27), MIXCAL), (WW (30), KURTOT),        &
      (WW (34), ITER), (WW (42), XK), (WW (65), XKRMM), (WW (77),       &
      IADJNT)                                                           
!                                                                       
!         CHECKL CALCULATES FOR EACH GROUP TOTAL ABSORPTION, REMOVAL,   
!         INSCATTERING, FISSION, AND THE NEUTRON BALANCE RATIO.         
!         CHECKL ALSO CALCULATES XKRMM (THE SECOND GUESS OF K-EFF)      
!         FROM NEUTRON BALANCE CONSIDERATIONS.                          
!                                                                       
      XESC = CHECK (1, 2) 
      CHECK (1, 2) = 0.D0 
      N = 0 
      NOD = NX + 1 
      DO 20 I1 = 1, NY 
         NOD = NOD+2 
         DO 20 J1 = 1, NX 
            NOD = NOD+1 
            MIX = MAP (NOD) 
            IF (MIX) 20, 20, 10 
   10       N = N + 1 
!                                                                       
            DO 16 I = 1, NG 
               REM = 0.D0 
               SCAT = 0.D0 
               FIS = 0.D0 
               DO 14 J = 1, NG 
!      ORDINARY SOLUTION:                                               
                  IA = I 
                  JA = J 
                  IF (IADJNT.NE.1) GOTO 12 
!      ADJOINT SOLUTION:                                                
                  IA = J 
                  JA = I 
   12             REM = REM + XSEC (JA, IA + 4, MIX) 
                  SCAT = SCAT + XSEC (IA, JA + 4, MIX) * FLUX (J, N) 
   14          FIS = FIS + XSEC (IA, 4, MIX) / XK * XSEC (JA, 3, MIX)   &
               * FLUX (J, N)                                            
               CHECK (I, 2) = CHECK (I, 2) - XSEC (I, 2, MIX) * FLUX (I,&
               N)                                                       
               CHECK (I, 3) = CHECK (I, 3) - REM * FLUX (I, N) 
               CHECK (I, 4) = CHECK (I, 4) + SCAT 
   16       CHECK (I, 5) = CHECK (I, 5) + FIS 
   20 CONTINUE 
!                                                                       
      XPRO = 0.D0 
      XABS = 0.D0 
      DO 18 I = 1, NG 
         XPRO = XPRO + CHECK (I, 5) 
         XABS = XABS - CHECK (I, 2) 
   18 CHECK (I, 6) = - (CHECK (I, 4) + CHECK (I, 5) ) / (CHECK (I, 1)   &
      + CHECK (I, 2) + CHECK (I, 3) ) - 1.D0                            
      XKRMM = XK * XPRO / (XABS + XESC) 
!                                                                       
      RETURN 
      END SUBROUTINE CHECKL                         
!OUTPUT                      31 OCT 1983         OUTPUT              LAB
!                                                                       
!                            ****************                           
!                            * OUTPUT/LABAN *                           
!                            ****************                           
!                                                                       
!         OUTPUT LINK OF LABAN                                          
!         --------------------                                          
!                                                                       
      SUBROUTINE OUTPUT (C, JPR, HELL) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION C ( * ), HELL ( * ) 
      COMMON / CLAB / WW (100) 
      COMMON / ADRS / AD (50) 
      COMMON / PWR / EPF, CPOW, VCORE 
!                                                                       
      EQUIVALENCE (AD (01), LNDX), (AD (02), LNDY), (AD (03), LHX),     &
      (AD (04), LHY), (AD (05), LCHI), (AD (06), LXSEC), (AD (09),      &
      LMAP), (AD (10), LCHECK), (AD (15), LF), (AD (16), LA), (AD (21), &
      LOLD), (AD (22), LCUR), (AD (23), LFLUX), (AD (24), LPOW),        &
      (AD (25), LOUT)                                                   
!                                                                       
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (07), NXY2),        &
      (WW (08), NG), (WW (09), NGNG), (WW (10), NG4), (WW (13), NQ),    &
      (WW (15), NGL), (WW (16), NGQ), (WW (17), NGL4), (WW (30),        &
      KURTOT), (WW (31), LUXTOT), (WW (61), PNORM), (WW (63), NP1),     &
      (WW (68), ANORM), (WW (69), NNORM)                                
!                                                                       
      IF (JPR) 10, 10, 50 
!                                                                       
!  CALCULATION OF NODAL AVERAGE POWER,                                  
!  PRINTING OF NEUTRON BALANCE, POWER, AVERAGE FLUXES AND CURRENTS.     
!                                                                       
!      ------------------------------------------------------           
   10 CALL POWERL (NG, NG4, NY, C (LMAP), C (LFLUX), C (LPOW), C (LHX), &
      C (LHY), C (LXSEC) )                                              
!      ------------------------------------------------------           
      IF (CPOW.NE.0.D0) THEN 
         TNORM = CPOW 
      ELSE 
         TNORM = PNORM 
      ENDIF 
      IF (NNORM.EQ.1) TNORM = ANORM 
!                                                                       
      CALL NORM (C (LPOW), NX * NY, PNORM) 
      CALL NORM (C (LFLUX), NG * LUXTOT, TNORM) 
      CALL NORM (C (LCUR), KURTOT, TNORM) 
      CALL NORM (C (LOLD), KURTOT, TNORM) 
      CALL NORM (C (LCHECK), 6 * NG, TNORM) 
!                                                                       
!      --------------------------------------------------------------   
      CALL EXPORT (C (LMAP), C (LCUR), C (LOLD), C (LFLUX), C (LPOW),   &
      C (LHX), C (LHY), NG, NY)                                         
!      --------------------------------------------------------------   
      CALL UT (C (LMAP), C (LNDX), C (LNDY), C (LHX), C (LHY), C (LCUR),&
      C (LOLD), C (LFLUX), C (LPOW), C (LCHECK), NG, NG4, NXY2, C (     &
      LXSEC) )                                                          
!      --------------------------------------------------------------   
      RETURN 
!                                                                       
!  CALCULATION OF DETAILED FLUX DISTRIBUTION.                           
!                                                                       
   50 LPLEG = LOUT 
      LFLU = LPLEG + NP1 * NQ 
      LGRID = LFLU + NGQ 
      LBMAT = LGRID+NP1 * NP1 
      LEVEC = LBMAT + NGNG 
      LPSI = LEVEC + NGNG * 2 
      LREST = LPSI + NGQ 
      IF (NP1.EQ.1) GOTO 60 
!      -------------------------                                        
      CALL LEGPOL (C (LPLEG), NQ) 
!      ---------------------------------------------------------------  
   60 CALL EDIT (C (LMAP), C (LF), C (LHX), C (LHY), C (LCUR), C (LOLD),&
      C (LFLU), C (LPLEG), C (LGRID), C (LXSEC), C (LBMAT), C (LEVEC),  &
      C (LPSI), C (LREST), NGL, NXY2, NG, NG4, HELL)                    
!      ---------------------------------------------------------------  
      RETURN 
!                                                                       
      END SUBROUTINE OUTPUT                         
!POWERL                      26 MAI 1982         OUTPUT              LAB
!                                                                       
!                            **********                                 
!                            * POWERL *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE POWERL (LG, NG4, LY, MAP, FLUX, POWER, HX, HY, XSEC) 
!                                                                       
!         POWERL CALCULATES THE AVERAGE POWER FOR EACH NODE.            
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP ( * ), FLUX (LG, * ), POWER (LY, * ), HX ( * ),     &
      HY ( * ), XSEC (LG, NG4, * )                                      
      COMMON / CLAB / WW (100) 
      COMMON / PWR / EPF, CPOW, VCORE 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (08), NG), (WW (11),&
      NL), (WW (30), KURTOT), (WW (61), PNORM), (WW (68), ANORM)        
!                                                                       
      POWTOT = 0.D0 
      ANORM = 0.D0 
      VOL = 0.D0 
      NOD = 0 
      K = NX + 1 
      DO 6 I = 1, NY 
         K = K + 2 
         DO 6 J = 1, NX 
            K = K + 1 
            MIX = MAP (K) 
            POWER (I, J) = 0.D0 
            IF (MIX) 6, 6, 2 
    2       NOD = NOD+1 
            V = HX (J) * HY (I) 
            POW = 0.D0 
            ABSO = 0.D0 
            DO 4 IG = 1, NG 
               ABSO = ABSO + XSEC (IG, 2, MIX) * FLUX (IG, NOD) 
!   4  POW =POW  + XSEC(IG,3,MIX)*FLUX(IG,NOD)                          
    4       POW = POW + XSEC (IG, 5 + NG, MIX) * FLUX (IG, NOD) 
            ANORM = ANORM + ABSO 
            IF (DABS (POW) .LT.1.D-20) GOTO 6 
            POW = EPF * POW 
            POWER (I, J) = POW / V 
            POWTOT = POWTOT + POW 
            VOL = VOL + V 
    6 CONTINUE 
!                                                                       
      ANORM = 1.D0 / ANORM 
      PNORM = DSQRT (DBLE (KURTOT / NL) ) 
      IF (DABS (POWTOT) .GT.1.D-12) PNORM = VOL / POWTOT 
      IF (CPOW.NE.0.D0) CPOW = CPOW / POWTOT 
      VCORE = VOL 
!                                                                       
      RETURN 
      END SUBROUTINE POWERL                         
!EXPORT                      26 MAI 1982         OUTPUT              LAB
!                                                                       
!                            **********                                 
!                            * EXPORT *                                 
!                            **********                                 
!                                                                       
! SAVE DATA ON FILE IFIL                                                
! SAVE RECORD: KEFF,NODES,(POW(I),I=1,NODES),(VOL(I),I=1,NODES)         
!              ON FILE JFIL                                             
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE EXPORT (MAP, CURA, CURB, FLUX, POWER, HX, HY, LG, LY) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MAXPV = 200) 
      DIMENSION POW (MAXPV), VOL (MAXPV), HX ( * ), HY ( * ) 
      DIMENSION MAP ( * ), CURA ( * ), CURB ( * ), FLUX (LG, * ),       &
      POWER (LY, * )                                                    
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (07), NXY2),        &
      (WW (08), NG), (WW (11), NL), (WW (26), NLOLD), (WW (30), KURTOT),&
      (WW (31), LUXTOT), (WW (64), IFIL), (WW (42), XKEF), (WW (62),    &
      MCASE), (WW (79), JFIL)                                           
      CHARACTER FILNAM * 32, FNAME * 32 
      COMMON / FILECM / FILNAM (4) 
      COMMON / VERSI / IVERS 
!                                                                       
      NLOLD = NL 
      IF (IFIL.LE.0) GOTO 100 
      REWIND IFIL 
!                                                                       
      WRITE (IFIL) (WW (I), I = 1, 100) 
      WRITE (IFIL) ( (POWER (I, J), J = 1, NX), I = 1, NY) 
      WRITE (IFIL) ( (FLUX (I, J), J = 1, LUXTOT), I = 1, NG) 
      WRITE (IFIL) (MAP (I), I = 1, NXY2 * 7) 
!                                                                       
      NP = NL - 1 
      DO 60 N = 1, KURTOT, NL 
   60 WRITE (IFIL) (CURA (I), I = N, N + NP) 
      DO 70 N = 1, KURTOT, NL 
   70 WRITE (IFIL) (CURB (I), I = N, N + NP) 
!                                                                       
  100 IF (JFIL.EQ.0) RETURN 
      IF (IVERS.EQ.1.AND.MCASE.EQ.1) THEN 
         FNAME = FILNAM (4) 
         OPEN (JFIL, FILE = FNAME, FORM = 'UNFORMATTED') 
!LPI     OPEN(JFIL,FILE=FNAME,FORM='UNFORMATTED',ORGANIZATION='DYNAMIC')
      ENDIF 
      IF (LUXTOT.LE.MAXPV) GOTO 101 
      WRITE (LUNFIL(2), 1030) MAXPV, LUXTOT 
      RETURN 
  101 CONTINUE 
      NCOUNT = 0 
      K = NX + 1 
      DO 103 I = 1, NY 
         K = K + 2 
         DO 103 J = 1, NX 
            K = K + 1 
            IF (MAP (K) ) 103, 103, 104 
  104       NCOUNT = NCOUNT + 1 
            POW (NCOUNT) = POWER (I, J) 
            VOL (NCOUNT) = HX (J) * HY (I) 
  103 CONTINUE 
      WRITE (JFIL) XKEF, NCOUNT, (POW (I), I = 1, NCOUNT), (VOL (I),    &
      I = 1, NCOUNT)                                                    
      WRITE (LUNFIL(2), 1010) JFIL 
!                                                                       
      RETURN 
 1010 FORMAT(/1X,'**** POWER MAP WAS WRITTEN TO UNIT ',I3,' ****') 
 1030 FORMAT(//1X,120(1H*)/1H0,I5,'>200 (TOO MANY NODES-CANNOT EXPORT'  &
     &,'POWER-MAP ON FILE JFIL'/1X,'INCREASE PARAMETER <MAXPV> IN SUBR',&
     &'<EXPORT> FROM',I7,' TO',I7,' AND RECOMPILE'/1H0,120(1H*))        
      END SUBROUTINE EXPORT                         
!UT                          19 SEP 1983         OUTPUT              LAB
!                                                                       
!                            ******                                     
!                            * UT *                                     
!                            ******                                     
!                                                                       
      SUBROUTINE UT (MAP, NDX, NDY, HX, HY, CURA, CURB, FLUX, POWFLU,   &
      CHECK, LG, NG4, NXY2, XSEC)                                       
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP (NXY2, 7), NDX ( * ), NDY ( * ), CURA ( * ),        &
      CURB ( * ), FLUX (LG, * ), POWFLU ( * ), CHECK (LG, 6), HX ( * ), &
      HY ( * ), XSEC (LG, NG4, * )                                      
      DIMENSION CURIN (4) 
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / PWR / EPF, CPOW, VCORE 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (03), NX2), (WW (04)&
      , NY2), (WW (08), NG), (WW (11), NL), (WW (59), ICUR), (WW (61),  &
      PNORM)                                                            
!                                                                       
!  PRINT NEUTRON BALANCE                                                
!                                                                       
      EFACT = 1.D-11 
      WRITE (LUNFIL(2), 200) 
      WRITE (LUNFIL(2), 202) 
      DO 1 I = 1, NG 
         CHECK (I, 6) = CHECK (I, 6) / PNORM 
    1 WRITE (LUNFIL(2), 204) I, (CHECK (I, J), J = 1, 6) 
!                                                                       
!  OUTPUT OF NODE AVERAGE POWERS RELATIVE TO CORE AVERAGE POWER         
!                                                                       
      WRITE (LUNFIL(2), 120) 
!      ------------------------------                                   
      CALL KARTA (POWFLU, NDX, NDY, NY) 
!      ------------------------------                                   
!                                                                       
!  OUTPUT OF NODE INTEGRATED POWERS RELATIVE TO CORE TOTAL POWER        
!                                                                       
      WRITE (LUNFIL(2), 121) 
      K = NX + 1 
      DO 11 I = 1, NY 
         K = K + 2 
         DO 11 J = 1, NX 
            K = K + 1 
            IJ = (J - 1) * NY + I 
            IF (MAP (K, 1) ) 11, 11, 10 
   10       IF (POWFLU (IJ) .EQ.0.D0) GOTO 11 
            POWFLU (IJ) = POWFLU (IJ) * HX (J) * HY (I) / VCORE 
   11 CONTINUE 
!      ------------------------------                                   
      CALL KARTA (POWFLU, NDX, NDY, NY) 
!      ------------------------------                                   
!                                                                       
!  OUTPUT OF AVERAGE FLUXES.                                            
!                                                                       
      DO 20 IG = 1, NG 
         WRITE (LUNFIL(2), 100) IG 
         NOD = 0 
         K = NX + 1 
         DO 14 I = 1, NY 
            K = K + 2 
            DO 14 J = 1, NX 
               K = K + 1 
               IJ = (J - 1) * NY + I 
               POWFLU (IJ) = 0.D0 
               IF (MAP (K, 1) ) 14, 14, 12 
   12          NOD = NOD+1 
               POWFLU (IJ) = FLUX (IG, NOD) / HX (J) / HY (I) 
!----EZM------OCT. 1986---USING SET HETEROGENEITY FACTORS               
               MIX = MAP (K, 1) 
               DISCON = XSEC (IG, 6 + NG, MIX) 
               IF (DISCON.LE.0.D0) DISCON = 1.D0 
               POWFLU (IJ) = EFACT * POWFLU (IJ) / DISCON 
!----EZM                                                                
   14    CONTINUE 
!      ------------------------------                                   
   20 CALL KARTA (POWFLU, NDX, NDY, NY) 
!      ------------------------------                                   
!                                                                       
!  OUTPUT OF CURRENTS.                                                  
!                                                                       
      IF (ICUR.EQ.0) GOTO 55 
      NNNLL = 1 
      IF (ICUR.GT.1) NNNLL = NL 
      DO 50 IG = 1, NG 
         WRITE (LUNFIL(2), 104) IG 
         DO 50 IL = 1, NNNLL 
            JSP = (IG - 1) * NL + IL 
            LLL = IL - 1 
            WRITE (LUNFIL(2), 106) LLL 
            NOD = 0 
            ISW = 1 
            DO 48 I1 = 1, NY2 
               WRITE (LUNFIL(2), 204) 
               JSW = ISW 
               DO 46 J1 = 1, NX2 
                  NOD = NOD+1 
                  IF (MAP (NOD, 1) ) 32, 46, 32 
   32             DO 40 LSIDE = 1, 4 
                     XX = 0.D0 
                     HH = HX (J1 - 1) 
                     IF (LSIDE.EQ.2.OR.LSIDE.EQ.4) HH = HY (I1 - 1) 
                     JS = MAP (NOD, LSIDE+3) 
                     IF (JS + 1) 40, 40, 34 
   34                JS = JS + JSP 
                     IF (JSW) 36, 36, 38 
   36                XX = CURB (JS) 
                     GOTO 40 
   38                XX = CURA (JS) 
!  40  CURIN(LSIDE)=XX/HH                                               
   40             CURIN (LSIDE) = XX / HH * EFACT 
                  IM = I1 - 1 
                  JM = J1 - 1 
                  WRITE (LUNFIL(2), 108) IM, JM, (CURIN (I), I = 1, 4) 
   46          JSW = - JSW 
   48       ISW = - ISW 
   50 CONTINUE 
!                                                                       
   55 RETURN 
!                                                                       
  100 FORMAT  (//1X,'AVERAGE FLUX (UNITS OF 10+11) GROUP',I3            &
     &         /1X,'**************************************')            
  101 FORMAT  (//' AVERAGE FLUX GROUP',I3/1X,21(1H*)) 
  104 FORMAT  (//1X,'IN-CURRENTS (UNITS OF 10+11) GROUP',I3,            &
     &         '.  ORTHONORMALIZED MOMENTS'/ 1X,40(1H-)//               &
     &         '  I  J',3X,'NORTH',7X,'EAST',8X,'SOUTH',7X,'WEST')      
  105 FORMAT  (//1X,'IN-CURRENTS GROUP',I3,                             &
     &         '.  ORTHONORMALIZED MOMENTS'/ 1X,20(1H-)//               &
     &         '  I  J',3X,'NORTH',7X,'EAST',8X,'SOUTH',7X,'WEST')      
  106 FORMAT  (/70X,I2,'-TH MOMENT') 
  108 FORMAT  (1X,I2,I3,4F12.7) 
  120 FORMAT  (//1X,'NODE AVERAGE POWERS RELATIVE TO CORE AVERAGE POWER'&
     &         /1X,'**************************************************')
  121 FORMAT  (//1X,'NODE INTEGRAL POWERS RELATIVE TO CORE TOTAL POWER'/&
     &          1X,'*************************************************') 
  200 FORMAT  (//' NEUTRON BALANCE GROUP-WISE.'/1X,26(1H*)) 
  202 FORMAT  (/' GROUP  LEAKAGE      ABSORPTION   REMOVAL      INSC',  &
     &         'ATTER    FISSION      BALANCE')                         
  204 FORMAT  (1X,I3,1X,0P5E13.4,0PF13.8) 
!                                                                       
      END SUBROUTINE UT                             
!KARTA                       19 SEP 1983         OUTPUT              LAB
!                                                                       
!                            *********                                  
!                            * KARTA *                                  
!                            *********                                  
!                                                                       
      SUBROUTINE KARTA (ARRAY, NDX, NDY, MY) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION ARRAY (MY, * ), NDX ( * ), NDY ( * ) 
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (05), NNX), (WW (06)&
      , NNY)                                                            
!                                                                       
!  PRINT A CORE MAP OF THE "ARRAY" DISTRIBUTION.                        
!  NNX=NUMBER OF NODES IN X-DIRECTION                                   
!  NNY=NUMBER OF NODES IN Y-DIRECTION                                   
!  NX=NUMBER OF NODE PARTITIONS (SUBNODES) IN X-DIRECTION               
!  NY=NUMBER OF NODE PARTITIONS (SUBNODES) IN Y-DIRECTION               
!  MY=NY                                                                
!                                                                       
      WRITE (LUNFIL(2), 100) (J, J = 1, NX) 
      WRITE (LUNFIL(2), 102) 
      DO 10 I = 1, NY 
         WRITE (LUNFIL(2), 104) I, (ARRAY (I, J), J = 1, NX) 
   10 END DO 
!                                                                       
!  IF THE CORE NODES ARE DIVIDED INTO SUBNODES, PRINT A CORE MAP OF THE 
!  ORIGINAL CORE.                                                       
!                                                                       
      IF (NNX.EQ.NX.AND.NNY.EQ.NY) RETURN 
      WRITE (LUNFIL(2), 100) (J, J = 1, NNX) 
      WRITE (LUNFIL(2), 102) 
      IY = 0 
      DO 40 I = 1, NNY 
         JX = 0 
         DO 30 J = 1, NNX 
            SUM = 0.D0 
            DO 20 II = 1, NDY (I) 
               IYY = IY + II 
               DO 20 JJ = 1, NDX (J) 
                  JXX = JX + JJ 
   20       SUM = SUM + ARRAY (IYY, JXX) 
            ARRAY (I, J) = SUM / DBLE (NDX (J) * NDY (I) ) 
   30    JX = JX + NDX (J) 
         WRITE (LUNFIL(2), 104) I, (ARRAY (I, J), J = 1, NNX) 
   40 IY = IY + NDY (I) 
!                                                                       
!                                                                       
  100 FORMAT  (/'  I/J',I6,10I10,/(5X,I6,10I10)) 
  102 FORMAT  (1X) 
  104 FORMAT  (1X,I2,2X,11F10.5,/(5X,11F10.5)) 
!                                                                       
      RETURN 
      END SUBROUTINE KARTA                          
!EDIT                        31 OCT 1983         OUTPUT              LAB
!                                                                       
!                            ********                                   
!                            * EDIT *                                   
!                            ********                                   
!                                                                       
      SUBROUTINE EDIT (MAP, F, HX, HY, CURA, CURB, FLUX, PLEG, GRID,    &
      XSEC, BMAT, EVEC, PSI, REST, LGL, NXY2, LG, NG4, X)               
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION F ( * ), CURA ( * ), CURB ( * ), FLUX ( * ), MAP (NXY2, &
      7), PLEG ( * ), HX ( * ), HY ( * ), GRID ( * ), XSEC (LG, NG4,    &
      * ), BMAT ( * ), EVEC ( * ), PSI ( * ), REST ( * ), X (LGL, 4)    
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (03), NX2), (WW (04)&
      , NY2), (WW (08), NG), (WW (11), NL), (WW (13), NQ), (WW (15),    &
      NGL), (WW (16), NGQ), (WW (20), NGLF), (WW (23), NQQ), (WW (30),  &
      KURTOT), (WW (31), LUXTOT), (WW (60), IFLUX), (WW (63), NP1),     &
      (WW (64), IFIL), (WW (78), IBUCK)                                 
!                                                                       
!  READ CURA AND CURB FROM IFIL.                                        
!                                                                       
      REWIND IFIL 
      READ (IFIL) 
      READ (IFIL) 
      READ (IFIL) 
      READ (IFIL) 
!                                                                       
      NP = NL - 1 
      DO 10 N = 1, KURTOT, NL 
   10 READ (IFIL) (CURA (I), I = N, N + NP) 
      DO 20 N = 1, KURTOT, NL 
   20 READ (IFIL) (CURB (I), I = N, N + NP) 
!                                                                       
!  LOOP=1   FLUX IN GROUP    REPRESENTATION                             
!  LOOP=2   FLUX IN BUCKLING REPRESENTATION                             
!                                                                       
      DO 200 LOOP = 1, 2 
         GOTO (22, 24), LOOP 
   22    IF (IFLUX.EQ.0) GOTO 200 
         WRITE (LUNFIL(2), 1010) NQ, NQ 
         GOTO 26 
   24    IF (IBUCK.EQ.0) GOTO 200 
         WRITE (LUNFIL(2), 1020) NQ, NQ 
   26    IF (NP1.GT.1) WRITE (LUNFIL(2), 1030) NP1, NP1 
         WRITE (LUNFIL(2), 1040) 
!                                                                       
!  PREPARE CALCULATION OF DETAILED FLUX DISTRIBUTION IN "FLUX2".        
!                                                                       
         ISW = 1 
         NOD = 0 
         DO 90 I1 = 1, NY2 
            JSW = ISW 
            DO 80 J1 = 1, NX2 
               NOD = NOD+1 
               MIX = MAP (NOD, 1) 
               IF (MIX) 80, 80, 30 
   30          DO 34 LSIDE = 1, 4 
                  KUR = MAP (NOD, LSIDE+3) 
                  DO 34 I = 1, NGL 
                     XX = CURA (I + KUR) 
                     IF (JSW.LT.0) XX = CURB (I + KUR) 
   34          X (I, LSIDE) = XX 
               MFY = MAP (NOD, 2) * NGLF + 1 
               MFX = MAP (NOD, 3) * NGLF + 1 
!      --------------------------------------------------------------   
               CALL FLUX2 (F (MFY), F (MFX), X (1, 1), X (1, 2),        &
               X (1, 3), X (1, 4), FLUX, NQ)                            
!      --------------------------------------------------------------   
!                                                                       
               HH = HX (J1 - 1) * HY (I1 - 1) 
               DO 36 I = 1, NGQ 
   36          FLUX (I) = FLUX (I) / HH 
!                                                                       
!  GET FLUX IN BUCKLING REPRESENTATION                                  
!                                                                       
               IF (LOOP.EQ.1) GOTO 38 
!      -------------------------------------------------------          
               CALL BUCKRE (XSEC (1, 1, MIX), FLUX, BMAT, EVEC, PSI,    &
               REST, NG, NG4, NQQ)                                      
!      -------------------------------------------------------          
!                                                                       
!  PRINT FLUX MATRIX AND GRID POINT VALUES                              
!                                                                       
   38          DO 40 IG = 1, NG 
                  II = IG 
                  NN = (IG - 1) * NQQ + 1 
!      --------------------------------------------------               
   40          CALL OPED (I1 - 1, J1 - 1, II, FLUX (NN), PLEG, GRID, NQ,&
               NP1)                                                     
!      --------------------------------------------------               
   80       JSW = - JSW 
   90    ISW = - ISW 
  200 END DO 
!                                                                       
 1010 FORMAT  (1X,'DETAILED FLUX DISTRIBUTION  -  ORTHONORMALIZED  ',   &
     &        'FLUX MATRIX  (',I2,'X' ,I2,')  IN GROUP ',               &
     &        'REPR. PRINTED COLUMNWISE '/1X,26('*'))                   
 1020 FORMAT  (1X,'DETAILED FLUX DISTRIBUTION  -  ORTHONORMALIZED  ',   &
     &        'FLUX MATRIX  (',I2,'X' ,I2,')  IN BUCKLING ',            &
     &        'REPR. PRINTED COLUMNWISE '/1X,26('*'))                   
 1030 FORMAT  (32X,'POINT VALUES GIVEN FOR A  (',I2,'X', I2,            &
     &         ')  GRID.')                                              
 1040 FORMAT  (/'  I  J  G') 
!                                                                       
      RETURN 
      END SUBROUTINE EDIT                           
!FLUX2                       26 MAI 1982         OUTPUT              LAB
!                                                                       
!                            *********                                  
!                            * FLUX2 *                                  
!                            *********                                  
!                                                                       
      SUBROUTINE FLUX2 (FY, FX, CUR1, CUR2, CUR3, CUR4, FLUX, LQ) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION FLUX (LQ, LQ, * ), CUR1 ( * ), CUR2 ( * ), CUR3 ( * ),  &
      CUR4 ( * ), FX ( * ), FY ( * )                                    
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (08), NG), (WW (11), NL), (WW (13), NQ), (WW (16),&
      NGQ)                                                              
!                                                                       
!  CALCULATION OF THE (NQ*NQ) FLUX EXPANSION COEFFICINT MATRIX FOR      
!  ACTUAL NODE.                                                         
!                                                                       
      CALL FAT (FLUX, 0.D0, NGQ) 
!                                                                       
      IY = 0 
      DO 90 IG = 1, NG 
         EJ = 1.D0 
         DO 80 J = 1, NQ 
            EI = 1.D0 
            DO 70 I = 1, NQ 
               EIJ = EI * EJ 
               SUMX = 0.D0 
               SUMY = 0.D0 
               JJ = 0 
               DO 60 JG = 1, NG 
                  EM = 1.D0 
                  DO 50 M = 1, NL 
                     JJ = JJ + 1 
                     IF (EM * EI) 50, 50, 40 
   40                IY = IY + 1 
                     SUMY = SUMY + FY (IY) * (EIJ * CUR1 (JJ) + CUR3 (  &
                     JJ) )                                              
                     SUMX = SUMX + FX (IY) * (EJ * CUR2 (JJ) + EI *     &
                     CUR4 (JJ) )                                        
   50             EM = - EM 
   60          END DO 
               FLUX (I, J, IG) = FLUX (I, J, IG) + SUMY 
               FLUX (J, I, IG) = FLUX (J, I, IG) + SUMX 
   70       EI = - EI 
   80    EJ = - EJ 
   90 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE FLUX2                          
!OPED                        31 OCT 1983         OUTPUT              LAB
!                                                                       
!                            ********                                   
!                            * OPED *                                   
!                            ********                                   
!                                                                       
      SUBROUTINE OPED (IM, JM, IG, FLUX, PLEG, GRID, LQ, NP1) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION FLUX ( * ), PLEG (LQ, * ), GRID (NP1, * ) 
!                                                                       
      COMMON / CLAB / WW (100) 
      COMMON / LUNSCM / LUNFIL (4) 
      EQUIVALENCE (WW (13), NQ), (WW (23), NQQ) 
!                                                                       
!  PRINT FLUX MATRIX.                                                   
!                                                                       
      WRITE (LUNFIL(2), 1000) IM, JM, IG, (FLUX (L), L = 1, NQQ) 
      IF (NP1.EQ.1) RETURN 
!                                                                       
!  CALCULATE FLUX POINTWISE AND FIND MAXIMUM FLUX IN NODE.              
!                                                                       
      FMAX = - 1.D+10 
      DO 20 IX = 1, NP1 
         DO 20 JY = 1, NP1 
            FSUM = 0.D0 
            DO 18 K = 1, NQ 
               DO 18 L = 1, NQ 
                  KL = (L - 1) * NQ + K 
   18       FSUM = FSUM + FLUX (KL) * PLEG (K, IX) * PLEG (L, JY) 
            GRID (NP1 + 1 - JY, IX) = FSUM 
            IF (FSUM.LT.FMAX) GOTO 20 
            FMAX = FSUM 
            IMAX = IX 
            JMAX = JY 
   20 CONTINUE 
!                                                                       
!  WRITE POINTWISE FLUX.                                                
!                                                                       
      WRITE (LUNFIL(2), 1010) FMAX, IMAX, JMAX 
      DO 40 JY = 1, NP1 
   40 WRITE (LUNFIL(2), 1020) (GRID (JY, IX), IX = 1, NP1) 
!                                                                       
 1000 FORMAT  (/1X,I2,2I3,11F9.3,/(9X,11F9.3)) 
 1010 FORMAT  (12X,10H---  MAX =,F9.4,3X,10H FOR X,Y = ,2I3,2X,3H---) 
 1020 FORMAT  ( 9X,11F9.3) 
!                                                                       
      RETURN 
      END SUBROUTINE OPED                           
!BUCKRE                      31 OCT 1983         OUTPUT              LAB
!                                                                       
!                            **********                                 
!                            * BUCKRE *                                 
!                            **********                                 
!                                                                       
!  PURPOSE : TO CONVERT FLUX MATRIX FROM GROUP REPRESENTATION TO        
!            BUCKLING REPRESENTATION.                                   
!  INCLUDED: FULL 1/D MATRIX FEB. 1988                                  
!            FOR THE BUCKLING MATRIX.                                   
!  REMARK  : IT IS HERE PRESUMED THAT NO IMAGINARY EIGENVALUES APPEAR   
!            FOR THE BUCKLING MATRIX: THUS NOT VALID FOR MORE THAS 2 GRP
!                                                                       
!      ---------------------------------------------------------------  
      SUBROUTINE BUCKRE (XSEC, FLUX, BMAT, EVEC, PSI, REST, NG, NG4,    &
      NQQ)                                                              
!      ---------------------------------------------------------------  
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION XSEC (NG, NG4), FLUX (NQQ, * ), BMAT (NG, * ), EVEC (NG,&
      * ), PSI (NQQ, * ), REST ( * )                                    
!                                                                       
      COMMON / CLAB / WW (100) 
      COMMON / DINV / IDIAG 
      EQUIVALENCE (WW (42), XK), (WW (77), IADJNT) 
!                                                                       
!  COMPUTE BUCKLING MATRIX                                              
!                                                                       
      DO 20 IG = 1, NG 
         REM = 0.D0 
         DIF = XSEC (IG, 1) 
         DO 15 JG = 1, NG 
!                                                                       
!      ORDINARY SOLUTION                                                
            IA = IG 
            JA = JG 
            IF (IADJNT.NE.1) GOTO 10 
!                                                                       
!      ADJOINT SOLUTION                                                 
            IA = JG 
            JA = IG 
!                                                                       
   10       REM = REM + XSEC (JA, IA + 4) 
            BMAT (IG, JG) = XSEC (IA, 4) / XK * XSEC (JA, 3) + XSEC (IA,&
            JA + 4)                                                     
            IF (IDIAG.EQ.0) BMAT (IG, JG) = BMAT (IG, JG) / DIF 
   15    END DO 
         IF (IDIAG.EQ.0) THEN 
            BMAT (IG, IG) = BMAT (IG, IG) - (XSEC (IG, 2) + REM)        &
            / DIF                                                       
         ELSE 
            BMAT (IG, IG) = BMAT (IG, IG) - XSEC (IG, 2) + REM 
         ENDIF 
   20 END DO 
!                                                                       
      IF (IDIAG.NE.0) THEN 
         CALL MULTAB (XSEC (1, 7 + NG), BMAT, REST, NG) 
         CALL EQUATE (BMAT, REST, NGNG) 
      ENDIF 
!                                                                       
!  COMPUTE EIGENVALUES. NORMALIZE SUCH THAT DIAGONAL EQUALS UNITY.      
!                                                                       
      CALL EIGENV (NG, BMAT, REST, REST (NG + 1), EVEC, EVEC (1, NG + 1)&
      , IND)                                                            
!                                                                       
      DO 30 JG = 1, NG 
         DO 30 IG = 1, NG 
   30 EVEC (IG, JG) = EVEC (IG, JG) / EVEC (JG, JG) 
!                                                                       
!  INVERT EIGENVECTOR MATRIX                                            
!                                                                       
      CALL MINV (EVEC, NG, DET, REST, REST (NG + 1) ) 
!                                                                       
!  COMPUTE FLUX IN BUCKLING REPRESENTATION                              
!                                                                       
      DO 50 N = 1, NQQ 
         DO 50 IG = 1, NG 
            SUM = 0.D0 
            DO 40 JG = 1, NG 
   40       SUM = SUM + EVEC (IG, JG) * FLUX (N, JG) 
   50 PSI (N, IG) = SUM 
!                                                                       
      DO 60 N = 1, NQQ 
         DO 60 IG = 1, NG 
   60 FLUX (N, IG) = PSI (N, IG) 
!                                                                       
      RETURN 
      END SUBROUTINE BUCKRE                         
!LEGPOL                      26 MAI 1982         OUTPUT              LAB
!                                                                       
!                            **********                                 
!                            * LEGPOL *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE LEGPOL (PLEG, LQ) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION PLEG (LQ, * ) 
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (13), NQ), (WW (63), NP1) 
!                                                                       
!  FIND THE LEGENDRE POLYNOMIAL VALUES AT THE NP1  POINTS  -1, ..., +1  
!  FOR THE FUNCTIONS UP TO ORDER (NQ-1).                                
!                                                                       
      DO 30 IX = 1, NP1 
         PLEG (1, IX) = 1.D0 
         IF (NQ.EQ.1) GOTO 30 
         X = - 1.D0 + DBLE (IX - 1) / DBLE (NP1 - 1) * 2.D0 
         PLEG (2, IX) = X 
         IF (NQ.EQ.2) GOTO 30 
         DO 10 L = 3, NQ 
            H = DBLE (2 * L - 3) * X * PLEG (L - 1, IX) - DBLE (L - 2)  &
            * PLEG (L - 2, IX)                                          
   10    PLEG (L, IX) = H / DBLE (L - 1) 
   30 END DO 
!                                                                       
      DO 40 L = 1, NQ 
         H = DSQRT (DBLE (2 * L - 1) ) 
         DO 40 IX = 1, NP1 
   40 PLEG (L, IX) = PLEG (L, IX) * H 
!                                                                       
      RETURN 
      END SUBROUTINE LEGPOL                         
!NORM                        26 MAI 1982         OUTPUT              LAB
!                                                                       
!                            ********                                   
!                            * NORM *                                   
!                            ********                                   
!                                                                       
      SUBROUTINE NORM (A, N, SCALAR) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ) 
      DO 10 I = 1, N 
   10 A (I) = SCALAR * A (I) 
      RETURN 
      END SUBROUTINE NORM                           
!FAT                         26 MAI 1982         OUTPUT              LAB
!                                                                       
!                            *******                                    
!                            * FAT *                                    
!                            *******                                    
!                                                                       
      SUBROUTINE FAT (A, SCALAR, N) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ) 
      DO 10 I = 1, N 
   10 A (I) = SCALAR 
      RETURN 
      END SUBROUTINE FAT                            
!RESCAL                      1  NOV 1982         RESCAL              LAB
!                                                                       
!                            ****************                           
!                            * RESCAL/LABAN *                           
!                            ****************                           
!                                                                       
!      CALCULATION OF MULTIGROUP RESPONSE MATRICES IN THIS LINK.        
!      --------------------------------------------------------         
!                                                                       
      SUBROUTINE RESCAL (C, MQ) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION C ( * ) 
      COMMON / WARN / IWARN 
      COMMON / CLAB / WW (100) 
      COMMON / ADRS / AD (50) 
      EQUIVALENCE (AD (03), LHX), (AD (04), LHY), (AD (06), LXSEC),     &
      (AD (08), LONE), (AD (09), LMAP), (AD (15), LF), (AD (16),        &
      LA), (AD (17), LT), (AD (18), LB), (AD (19), LRES)                
      EQUIVALENCE (WW (07), NXY2), (WW (08), NG), (WW (09), NGNG),      &
      (WW (10), NG4), (WW (11), NL), (WW (12), NLR), (WW (13), NQ),     &
      (WW (14), LQX), (WW (16), NGQ), (WW (20), NGLF), (WW (21),        &
      NGLR), (WW (22), NGLR2), (WW (23), NQQ)                           
!                                                                       
      NQ = MQ 
      NQQ = NQ * NQ 
      NGQ = NG * NQQ 
      NGLF = NGNG * NQ * ( (NL * NQ + 1) / 2) 
      LQX = MAX (NLR, NQ) 
      NFSIZ = 2 * NGNG * NLR * NQQ 
!                                                                       
      IWARN = 0 
!                                                                       
      LAH = LRES 
      LTH = LAH + NGLR2 
      LBH = LTH + NGLR2 
      LBA = LBH + NGLR2 
      LFH = LBA + NGLR2 
      LHELP = LFH + NFSIZ 
!      LAST1=LHELP+ 2*NGNG + 2*NGLR2 + NFSIZ + LH1 + LH2 + 2*NG +NG     
!                                                                       
!      -------------------------------------------------------          
      CALL RESADM (NG, NG4, NGLR, NGQ, NXY2, C (LMAP), C (LHX), C (LHY),&
      C (LXSEC), C (LA), C (LT), C (LB), C (LF), C (LONE), C (LAH),     &
      C (LTH), C (LBH), C (LBA), C (LFH), C (LHELP) )                   
!      -------------------------------------------------------          
!TEST                                                                   
!              WRITE(6,5555)                                            
!5555  FORMAT(1X,'AFTER RESADM')                                        
!      NGGL = NG*NG*NLR*NLR                                             
!      DO 9999 I=1,NGGL                                                 
!              IAH = LAH + I - 1                                        
!              ITH = LTH + I - 1                                        
!              IBH = LBH + I - 1                                        
!              IBA = LBA + I - 1                                        
!9999          WRITE(6,7777) C(IAH),C(ITH),C(IBH),C(IBA)                
!7777  FORMAT(1X,6(F10.6,1X))                                           
!8888  FORMAT(1X,8(F10.6,1X))                                           
!      DO 6666 I=1,NGGL                                                 
!              IFH = LFH + I - 1                                        
!6666          WRITE(6,8888) C(IFH)                                     
!              STOP                                                     
!TEST                                                                   
!                                                                       
      RETURN 
      END SUBROUTINE RESCAL                         
!RESADM                      1  NOV 1982         RESCAL              LAB
!                                                                       
!                            ****************                           
!                            * RESADM/LABAN *                           
!                            ****************                           
!                                                                       
      SUBROUTINE RESADM (LG, NGX, LGLR, NGQ, NXY2, MAP, HX, HY, XSEC, A,&
      T, B, F, ONE, AH, TH, BH, BA, FH, HELP)                           
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION MAP (NXY2, 7), HX ( * ), HY ( * ), XSEC (LG, NGX,       &
      * ), A ( * ), T ( * ), B ( * ), F ( * ), ONE (LGLR, * ), AH (LGLR,&
      * ), TH (LGLR, * ), BH (LGLR, * ), BA (LGLR, * ), FH (NGQ,        &
      * ), HELP ( * )                                                   
      COMMON / CLAB / WW (100) 
      EQUIVALENCE (WW (01), NX), (WW (02), NY), (WW (08), NG), (WW (09),&
      NGNG), (WW (11), NL), (WW (12), NLR), (WW (13), NQ), (WW (18),    &
      NGLB), (WW (19), NGLA), (WW (20), NGLF), (WW (21), NGLR), (WW (22)&
      , NGLR2), (WW (23), NQQ), (WW (24), LH1), (WW (25), LH2), (WW (27)&
      , MIXCAL), (WW (34), ITER), (WW (36), EPSR), (WW (42), XK),       &
      (WW (77), IADJNT)                                                 
      LOGICAL LSW 
!                                                                       
!         RESADM ADMINISTRATES THE CALCULATION OF RESPONSE MATRICES.    
!         ---------------------------------------------------------     
      LE = 1 + NGNG * 2 + NG * 2 + NG 
      LW = LE+2 * NGNG * NLR * NQQ 
      LWA = LW + NGLR2 
      LHE1 = LWA + NGLR2 
      LHE2 = LHE1 + LH1 
!      LAST = LHE2 + LH2                                                
!      LH1  = NGNG*MAX(NLR**2, 6+2*LQX + MAX(4,2+LQX))                  
!      LH2  = NGNG*NLR*MAX(NLR,NQQ)                                     
      LSW = ITER.EQ.1.OR.NQ.GT.1 
!                                                                       
      KA = 0 
      KB = 0 
      KF = 0 
      M0 = - 1 
!                                                                       
!         FIND FOR WHICH CELLS TO CALCULATE RESP MATRICES.              
!                                                                       
      NOD = NX + 1 
      DO 200 I = 1, NY 
         NOD = NOD+2 
         DO 200 J = 1, NX 
            NOD = NOD+1 
            MIX = MAP (NOD, 1) 
!                                                                       
            IF (MIX.LE.0) GOTO 200 
            M1 = MAP (NOD, 2) 
            IF (M1.LE.M0) GOTO 200 
            M0 = M1 
            LEQU = 0 
            IF (MAP (NOD, 2) .NE.MAP (NOD, 3) ) LEQU = 1 
!        IF MIXTURE IS VALID FOR EXTERNAL RESP MATRICES, THEN SKIP CALCU
            IF (MIX.GT.MIXCAL) GOTO 35 
!                                                                       
!       FIND DIFF COEFF DD0 TO BE USED IN PROVISIONAL BOUNDARY CONDITION
!                                                                       
            FIS = 0.D0 
            DD0 = 0.D0 
            DO 30 IG = 1, NG 
               FIS = FIS + XSEC (IG, 3, MIX) 
   30       DD0 = MAX (DD0, XSEC (IG, 1, MIX) ) 
!                                                                       
!        SHOULD CALCULATION FROM PREVIOUS OUTER ITER BE RETAINED FOR THI
            IF (LSW.OR.FIS.GT.1.D-10) GOTO 40 
!        YES.                                                           
   35       N2 = LEQU + 1 
            KA = KA + N2 * NGLA 
            KB = KB + N2 * NGLB 
            KF = KF + N2 * NGLF 
            GOTO 200 
   40       CONTINUE 
!                                                                       
!         CALL RESPMG FOR CALCULATION OF RESP MAT.                      
!                                                                       
!      ---------------------------------------------------------------- 
            CALL RESPMG (NG, NLR, NQ, HX (J), HY (I), DD0, XSEC (1, 1,  &
            MIX), XK, EPSR, LEQU, IADJNT, ONE, HELP (1), AH, TH, BH, BA,&
            FH, HELP (LE), HELP (LW), HELP (LWA), HELP (LHE1), HELP (   &
            LHE2) )                                                     
!      ---------------------------------------------------------------- 
!TEST                                                                   
!         CHECK RESPONSE MATRICES                                       
!         -----------------------                                       
!      IAD1=1+NGLR                                                      
!      CALL OUTFL(AH,NGLR,'A(Y)  ')                                     
!      IF(LEQU.NE.0) CALL OUTFL(AH(1,IAD1),NGLR,'A(X)  ')               
!      CALL OUTFL(TH,NGLR,'T(Y)  ')                                     
!      IF(LEQU.NE.0) CALL OUTFL(TH(1,IAD1),NGLR,'T(X)  ')               
!      CALL OUTFL(BH,NGLR,'B+(Y) ')                                     
!      IF(LEQU.NE.0) CALL OUTFL(BH(1,IAD1),NGLR,'B+(X) ')               
!      CALL OUTFL(BA,NGLR,'B-(Y) ')                                     
!      IF(LEQU.NE.0) CALL OUTFL(BA(1,IAD1),NGLR,'B-(X) ')               
!      IAD2=1+NGLR                                                      
!      LENGF=NGQ*NGLR                                                   
!      CALL OUTFLV(FH,LENGF,'F(Y)  ')                                   
!      IF(LEQU.NE.0) CALL OUTFLV(FH(1,IAD2),LENGF,'F(X)  ')             
!TEST                                                                   
!                                                                       
!         STORE MATRICES ROWWISE. NO STORING OF ZEROS.                  
!                                                                       
            J1S = 1 
            J1E = NGLR 
   50       DO 110 I1 = 1, NGLR 
               IF (I1 - (I1 - 1) / NLR * NLR - NL) 60, 60, 110 
   60          DO 100 J1 = J1S, J1E 
                  IF (J1 - (J1 - 1) / NLR * NLR - NL) 70, 70, 100 
   70             KB = KB + 1 
                  J2 = J1 - J1S + 1 
                  IF (ONE (I1, J2) ) 80, 80, 90 
   80             B (KB) = BA (I1, J1) 
                  GOTO 100 
   90             B (KB) = BH (I1, J1) 
                  KA = KA + 1 
                  A (KA) = AH (I1, J1) 
                  T (KA) = TH (I1, J1) 
  100          END DO 
  110       END DO 
!                                                                       
            DO 140 I1 = 1, NGQ 
               IQQ = MOD (I1 - 1, NQ) + 1 
               IQQ = MOD (IQQ, 2) 
               IF (IQQ.EQ.0) IQQ = - 1 
               FQQ = DBLE (IQQ) 
               DO 140 J1 = J1S, J1E 
                  IF (J1 - (J1 - 1) / NLR * NLR - NL) 120, 120, 140 
  120             J2 = J1 - J1S + 1 
                  IF (FQQ * ONE (1, J2) ) 140, 140, 130 
  130             KF = KF + 1 
                  F (KF) = FH (I1, J1) 
  140       CONTINUE 
            IF (LEQU.EQ.0) GOTO 200 
            J1S = NGLR + 1 
            J1E = 2 * NGLR 
            LEQU = 0 
            GOTO 50 
  200 CONTINUE 
!                                                                       
      RETURN 
      END SUBROUTINE RESADM                         
!RESPMG                      31 OCT 1983                             RES
!                                                                       
!                            **********                                 
!                            * RESPMG *                                 
!                            **********                                 
!                                                                       
!                                                                       
!   RESPMG CALCULATES MULTIGROUP RECTANGULAR DIFFUSION RESPONSE MATRICES
!                                                                       
!   INPUT DATA :                                                        
!      NG    = NR OF GROUPS.                                            
!      NL    = L+1 ; L=ORDER OF APPROXIMATION FOR RESPONSE MATRICES     
!      NQ    = FLUX EXPANSION ORDER  + 1.                               
!      HX    = WIDTH OF RECTANGLE IN X-DIRECTION.                       
!      HY    = WIDTH OF RECTANGLE IN Y-DIRECTION.                       
!      DD0   = D-VALUE USED IN PROVISIONAL BOUNDARY CONDITION.          
!      XSEC  = CONTAINS CROSS SECTIONS AS FOLLOWS.                      
!              FIRST COLUMN : GROUP DIFFUSION CONSTANTS.                
!              SECOND COLUMN : GROUP ABSORPTION CROSS SECTIONS.         
!              THIRD COLUMN : GROUP NU * FISSION CROSS SECTION.         
!              FOURTH COLUMN : FISSION SPECTRUM                         
!              NG COLUMNS BEYOND 4TH COLUMN : SCATTERING CROSS SECTION  
!              COLUMN NG + 5 : FISSION CROSS SECTION                    
!              COLUMN NG + 6 : DISCONTINUITY FACTORS                    
!              NG COLUMNS BEYOND COLUMN NG+6 : 1/D MATRIX IF IDIAG=1    
!      XK    = K-EFF.                                                   
!      EPS   = REQUESTED ACCURACY OF THE RESPONSE MATRICES.             
!      LEQU  = 0 IF HX=HY.  LEQU = 1  IF  HX.NE.HY .                    
!      IADJNT= 0 / 1   ORDINARY / ADJOINT SOLUTION.                     
!      VECT, W, WA, H1, H2, FH, ONE  ARE WORK MATRICES.                 
!                                                                       
!   OUTPUT DATA: IF LEQU=0 THEN X-DIRECTION MATRICES ARE NOT CALCULATED.
!      A = THE REFLECTION MATRICES A(Y) AND A(X) STORED COLUMNWISE AND  
!          SEQUENTIALLY.                                                
!      T = THE TRANSMISSION MATRICES T(Y) AND T(X).                     
!      B = THE RIGHT-SIDE-TRANSMISSION MATRICES B+(Y) AND B+(X).        
!          B HAS THE SAME STRUCTURE AS A.                               
!      BA = THE LEFT-SIDE-TRANSMISSION MATRICES B-(Y) AND B-(X).        
!           THE ENTRIES OF BA ARE NON-ZERO WHERE THE ENTRIES OF B ARE   
!           ZERO AND VICE VERSA.                                        
!      F = THE CURRENT-TO-FLUX MATRICES F(Y) AND F(X).                  
!                                                                       
!   DIMENSIONING OF ARRAYS :                                            
!      K=(NG*NL)**2                                                     
!      J=NL*(NG*NQ)**2                                                  
!      M1=MAX(K, NG*NG*(6+2*LQX+MAX(4,LQX+2)))  ,  LQX=MAX(NQ,NL)       
!      M2=NG*NG*NL*MAX(NL,NQQ)                                          
!      VECT(2*NG*NG+2*NG+NG), XSEC(NG*(NG4)),   CHI(NG)                 
!      A(2K),  T(2K),  B(2K),  BA(2K),  W(2K),  WA(2K),  ONE(K)         
!      F(2J),  FH(2J)                                                   
!      H1(M1), H2(M2)                                                   
!                                                                       
!      ---------------------------------------------------------------- 
      SUBROUTINE RESPMG (LG, LL, LQ, HX, HY, DD0, XSEC, XK, EPS, LEQU,  &
      IADJNX, ONE, VECT, A, T, B, BA, F, FH, W, WA, H1, H2)             
!      ---------------------------------------------------------------- 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION XSEC (LG, * ), VECT ( * ), A ( * ), T ( * ), B ( * ),   &
      BA ( * ), F ( * ), FH ( * ), W ( * ), WA ( * ), H1 ( * ), H2 ( * )&
      , ONE ( * )                                                       
      COMMON / CRMG / NG, NL, NQ, NQQ, NGL, NGQ, LQX, IADJNT 
      COMMON / DINV / IDIAG 
!                                                                       
      FSQRTF (I) = DSQRT (DBLE (2 * I - 1) ) 
!                                                                       
!  INITIATE                                                             
!                                                                       
      NG = LG 
      NGNG = NG * NG 
      NL = LL 
      NQ = LQ 
      NQQ = NQ * NQ 
      NGL = NG * NL 
      NGQ = NG * NQQ 
      LQX = MAX (NL, NQ) 
      IADJNT = IADJNX 
      IA2 = 2 
      IA3 = 3 
      IA4 = 4 
!                                                                       
!      -----------------------                                          
      CALL ONECAL (NG, NL, ONE) 
!      -----------------------                                          
!----EZM                                                                
!ALCULATE DD = (DD0*DINV - E)/(DD0*DINV + E)            -> H2           
      IF (IDIAG.NE.0) THEN 
         IM1 = 1 + NGNG 
         IM2 = IM1 + NGNG 
         IM3 = IM2 + NG 
         CALL EQUATE (W, XSEC (1, NG + 7), NGNG) 
         CALL MULTA (W, DD0, NGNG) 
         CALL SETID (WA, NG) 
         CALL DIFAB (W, WA, H1, NGNG) 
         CALL SUMAB (W, WA, H1 (IM1), NGNG) 
         CALL MINV (H1 (IM1), NG, DET, H1 (IM2), H1 (IM3) ) 
         CALL MULTAB (H1, H1 (IM1), H2, NG) 
      ENDIF 
!----EZM                                                                
!                                                                       
!  CALCULATE THE BUCKLING MATRIX AND STORE IN A()                       
!                                                                       
      IT = 1 
      DO 5 I = 1, NG 
         XBD = 0.25D0 / XSEC (I, 1) 
         IS = IT 
         REM = 0.D0 
         DO 1 K = 1, NG 
            IF (IADJNT.NE.1) REM = REM + XSEC (K, I + 4) 
    1    IF (IADJNT.EQ.1) REM = REM + XSEC (I, K + 4) 
         DO 4 J = 1, NG 
!                                                                       
!      ORDINARY SOLUTION :                                              
            IA = I 
            JA = J 
            IF (IADJNT.NE.1) GOTO 11 
!      ADJOINT SOLUTION :                                               
            IA = J 
            JA = I 
!                                                                       
   11       A (IS) = - XSEC (IA, IA4) / XK * XSEC (JA, IA3) - XSEC (IA, &
            JA + 4)                                                     
            IF (I - J) 3, 2, 3 
    2       A (IS) = A (IS) + XSEC (I, IA2) + REM 
    3       IF (IDIAG.EQ.0) A (IS) = A (IS) * XBD 
    4    IS = IS + NG 
    5 IT = IT + 1 
      IF (IDIAG.NE.0) THEN 
         CALL MULTAB (XSEC (1, 7 + NG), A, H1, NG) 
         CALL MULTA (H1, 0.25D0, NGNG) 
         CALL EQUATE (A, H1, NGNG) 
      ENDIF 
!                                                                       
!  CALCULATE EIGENVALUE SPECTRUM OF BUCKLING MATRIX                     
!                                                                       
      IADV1 = 1 + NG 
      IADV2 = IADV1 + NG 
      IADV3 = IADV2 + NGNG 
      IADIND = IADV3 + NGNG 
!      -------------------------------------------------------          
      CALL EIGENV (NG, A, VECT, VECT (IADV1), VECT (IADV2), VECT (IADV3)&
      , VECT (IADIND) )                                                 
!      -------------------------------------------------------          
      HHX = HX 
      HHY = HY 
      I1 = 0 
      I2 = 0 
      LA = 1 
      LF = 1 
!                                                                       
!   CALCULATE NECESSARY NR OF ROOTS FOR THE ACCURACY EPS.               
!                                                                       
    6 A0 = 0.25D0 * HHX / DD0 
      AA = A0 * A0 
      XX = 0.D0 
      DO 7 I = 1, NG 
         AG = 0.25D0 * HHX / XSEC (I, 1) 
         XJ = AA * DABS (5.D0 * A0 + AA - AG * AG) 
    7 IF (XJ.GT.XX) XX = XJ 
      JJ = (0.00083D0 * XX / EPS) **0.2D0 + 1.0D0 
      JJ = MIN (JJ, 32, IFIX (SNGL ( (1400.D0 * A0 * A0) **0.333D0 +    &
      1.D0) ) )                                                         
!                                                                       
!   CALL PRORES FOR CALCULATION OF PROVISIONAL RESPONSE MATRICES.       
!                                                                       
!      --------------------------------------------------               
      CALL PRORES (HHX, HHY, DD0, XSEC, VECT, VECT (IADV2), A (LA),     &
      T (LA), B (LA), FH (LF), H1, H2 (1 + NGNG), JJ)                   
!      --------------------------------------------------               
!                                                                       
!   CALL HETSMG FOR IMPROVEMENT OF THE PROVISIONAL RESP MATRICES.       
!                                                                       
!      -----------------------------------------------------------      
      CALL HETSMG (NGL, NGQ, HHX, HHY, DD0, XSEC, JJ, A (LA), B (LA),   &
      FH (LF), H2 (1 + NGNG) )                                          
!      -----------------------------------------------------------      
!                                                                       
!   CALL BALPRO TO FORCE NEUTRON CONSERVATION TO HOLD.                  
!                                                                       
!      ---------------------------------------------------------------- 
      CALL BALPRO (NG, NGL, NGQ, DD0, XSEC, XK, A (LA), T (LA), B (LA), &
      FH (LF), H1, H2 (1 + NGNG) )                                      
!      ---------------------------------------------------------------- 
!                                                                       
!TEST                                                                   
!      NGGL = NG*NG*NL*NL                                               
!      NGGQ = NG*NG*NQ*NQ*NL*2                                          
!              WRITE(6,1555)                                            
!1555  FORMAT(1X,'AFTER BALPRO IN RESPMG')                              
!      DO 1666 I=1,NGGL                                                 
!      IIII = LA + I - 1                                                
!1666          WRITE(6,1777) A(IIII),T(IIII),B(IIII)                    
!              WRITE(6,1888) (FH(I),I=1,NGGQ)                           
!TEST                                                                   
!   NORMALIZE THE RESPONSE MATRICES AND CALCULATE THE ERROR CURRENTS AND
!   STORE THEM IN W AND WA.                                             
!                                                                       
      IF (IDIAG.NE.0) GOTO 24 
      DO 20 M = 1, NG 
         EJ = 1.D0 
         DO 18 J = 1, NL 
            XJ = FSQRTF (J) 
            DO 16 N = 1, NG 
               DD = (DD0 - XSEC (N, 1) ) / (DD0 + XSEC (N, 1) ) 
               EI = 1.D0 
               DO 14 I = 1, NL 
                  XX = XJ * FSQRTF (I) 
                  EX = EI * XX 
                  DX = DD * XX 
                  EE = EJ * EI 
                  I1 = I1 + 1 
                  T (I1) = T (I1) * XX 
                  A (I1) = A (I1) * EI * XX 
                  IF (M.EQ.N.AND.I.EQ.J) A (I1) = A (I1) - EI 
                  IF (EE) 10, 10, 8 
    8             W (I1) = DX * B (I1) 
                  B (I1) = EX * B (I1) 
                  BA (I1) = 0.D0 
                  WA (I1) = 0.D0 
                  GOTO 14 
   10             WA (I1) = DX * B (I1) 
                  BA (I1) = EX * B (I1) 
                  B (I1) = 0.D0 
                  W (I1) = 0.D0 
   14          EI = - EI 
!                                                                       
               DO 16 I = 1, NQ 
                  XX = XJ * FSQRTF (I) 
                  DO 16 II = 1, NQ 
                     I2 = I2 + 1 
   16       FH (I2) = FH (I2) * XX * FSQRTF (II) 
   18    EJ = - EJ 
   20 END DO 
!----EZM                                                                
      GOTO 29 
!ALCULATE W(L,L') = DD * B(L,L')                                        
   24 CALL EQUATE (H1, B (LA), NGL * NGL) 
      DO 49 M = 1, NG 
         EJ = 1.D0 
         DO 48 J = 1, NL 
            XJ = FSQRTF (J) 
            DO 46 N = 1, NG 
               EI = 1.D0 
               I11 = (J - 1) * NGL + (M - 1) * NGL * NL 
               DO 44 I = 1, NL 
                  XX = XJ * FSQRTF (I) 
                  EX = EI * XX 
                  EE = EJ * EI 
                  I1 = I1 + 1 
                  I11 = I11 + 1 
                  T (I1) = T (I1) * XX 
                  A (I1) = A (I1) * EI * XX 
                  IF (M.EQ.N.AND.I.EQ.J) A (I1) = A (I1) - EI 
                  SOM = 0.D0 
                  DO 59 IG = 1, NG 
                     IDD = (IG - 1) * NG + N 
                     DD = H2 (IDD) 
                     DX = DD * XX 
                     I3 = I11 + NL * (IG - 1) 
                     SOM = SOM + DX * H1 (I3) 
   59             END DO 
                  IF (EE) 60, 60, 58 
   58             W (I1) = SOM 
                  B (I1) = EX * B (I1) 
                  BA (I1) = 0.D0 
                  WA (I1) = 0.D0 
                  GOTO 44 
   60             WA (I1) = SOM 
                  BA (I1) = EX * B (I1) 
                  B (I1) = 0.D0 
                  W (I1) = 0.D0 
   44          EI = - EI 
!                                                                       
               DO 46 I = 1, NQ 
                  XX = XJ * FSQRTF (I) 
                  DO 46 II = 1, NQ 
                     I2 = I2 + 1 
   46       FH (I2) = FH (I2) * XX * FSQRTF (II) 
   48    EJ = - EJ 
   49 END DO 
!                                                                       
!   IF LEQU.NE.0 , PREPARE FOR A SECOND CALL OF PRORES TO CALCULATE THE 
!   PROVISIONAL X-DIRECTION RESPONSE MATRICES.                          
!                                                                       
   29 IF (LA.GT.1) GOTO 30 
      LA = LA + NGL * NGL 
      LF = LF + NGL * NGQ 
      HHX = HY 
      HHY = HX 
      IF (LEQU.NE.0) GOTO 6 
   30 IF (NG.GT.1) GOTO 35 
      N2 = NGL * NGQ * (LEQU + 1) 
      DO 32 N = 1, N2 
   32 F (N) = FH (N) 
      GOTO 50 
!                                                                       
!   CALCULATE THE "EXACT" RESPONSE MATRICES.                            
!                                                                       
   35 LY = 1 
      LX = LY + NGL * NGL * LEQU 
      LFY = 1 
      LFX = LFY + NGL * NGQ * LEQU 
!TEST                                                                   
!              WRITE(6,5555) LY,LX                                      
!5555  FORMAT(1X,'AFTER W BEFORE MAAT IN RESPMG: LY=',I3,';LX=',I3)     
!      NGGL = NG*NG*NL*NL*(1+LEQU)                                      
!      NGGQ = NG*NG*NQ*NQ*NL*(1+LEQU)                                   
!      DO 6666 I=1,NGGL                                                 
!6666          WRITE(6,1777) A(I),T(I),B(I),BA(I),W(I),WA(I)            
!              WRITE(6,1888) (FH(I),I=1,NGGQ)                           
!TEST                                                                   
!      ---------------------------------------------------------        
   40 CALL IPINV (H1, H2, W (LX), W (LY), NGL, NG, NL, - 4.D0) 
      CALL MAF (NGQ, NGL, FH (LFY), FH (LFX), W (LY), H1, H2, F (LFY) ) 
      CALL MAB (NGL, + 1.D0, ONE, B (LY), A (LX), T (LX), W (LY),       &
      H1, H2)                                                           
      CALL IPINV (H1, H2, WA (LX), WA (LY), NGL, NG, NL, + 4.D0) 
      CALL MAB (NGL, - 1.D0, ONE, BA (LY), A (LX), T (LX), WA (LY),     &
      H1, H2)                                                           
      CALL MAAT (NGL, ONE, T (LX), A (LX), B (LY), BA (LY), W (LX),     &
      WA (LX) )                                                         
!      ---------------------------------------------------------        
      LY = LX 
      LX = 1 
      LFY = LFX 
      LFX = 1 
      IF (LY.NE.1) GOTO 40 
!                                                                       
   50 RETURN 
!              WRITE(6,2222) LY                                         
!2222  FORMAT(1X,'AFTER MAAT IN RESPMG: LY   = ',I3)                    
!      DO 1111 I=1,NGGL                                                 
!1111          WRITE(6,1777) A(I),T(I),B(I),BA(I),W(I),WA(I)            
!              WRITE(6,1888) (FH(I),F(I),I=1,NGGQ)                      
!1777  FORMAT(1X,6(F10.6,1X))                                           
!1888  FORMAT(1X,8(F10.6,1X))                                           
      END SUBROUTINE RESPMG                         
!PRORES                      26 MAI 1982                             RES
!                                                                       
!                            **********                                 
!                            * PRORES *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE PRORES (HX, HY, DD0, DIF, EIG, VECT, A, T, B, F, H, H2,&
      JJ)                                                               
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MEMR = 200, MEMPSX = 10) 
      DIMENSION R (MEMR), PSX (MEMPSX), PSH (MEMPSX), DIF ( * ),        &
      EIG ( * ), VECT ( * ), A ( * ), T ( * ), B ( * ), F ( * ),        &
      H ( * ), H2 ( * )                                                 
      COMMON / CRMG / NG, NL, NQ, NQQ, NGL, NGQ, LQX, IADJNT 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / DINV / IDIAG 
      COMMON / WARN / IWARN 
!                                                                       
!   PRORES CALCULATES PROVISIONAL MULTIGROUP RECTANGULAR RESPONSE MATRIX
!                                                                       
      MCHECK = MAX (NL, NQ) 
      IF (MCHECK.GT.MEMPSX) THEN 
         WRITE (LUNFIL(2), 12345) MEMPSX, MCHECK 
         IWARN = 1 
12345 FORMAT (/1X,'LEGENDRE MOMENTS GREATER THAN',I2,' EXCEEDS MEMORY IN&
     & PRORES',/                                                        &
     & 1X,'CHANGE MEMPSX TO',I3)                                        
      ENDIF 
      NNG = NG * NG 
      DI = 1.D0 / DD0 
      AD = 0.25D0 * HX * DI 
      ADD = AD * AD 
      HYX = HY / HX 
      R1 = HYX * HYX 
      MG = 2 * NG 
      NGNG = NG * NG 
      NGL2 = NG * (LQX + 2) 
      LQ = 1 + MAX ( (MG * MG), (NG * NGL2) ) 
      LW = LQ + MG * MG 
      LC = LW + MG * NG * LQX 
!      TOTAL STORAGE    IN H = NG*NG*(6+2*LQX+MAX(4,LQX+2))             
!      TOTAL STORAGE    IN H2= MAX(NG*NG*3+2*NG , NG*NG + NG*NG*NL*NL)  
!                                                                       
!ALCULATE DD = (D/DD0 + E)                            -> H2             
      IF (IDIAG.NE.0) THEN 
         IAD = 6 * NG + NGNG + 1 
         IM1 = 1 + NGNG 
         IM2 = IM1 + NGNG 
         IM3 = IM2 + NGNG 
         IM4 = IM3 + NG 
         DDX = 1.D0 / DD0 
         CALL EQUATE (H2 (IM1), DIF (IAD), NGNG) 
         CALL MINV (H2 (IM1), NG, DET, H2 (IM3), H2 (IM4) ) 
         CALL MULTA (H2 (IM1), DDX, NGNG) 
         CALL SETID (H2 (IM2), NG) 
         CALL SUMAB (H2 (IM1), H2 (IM2), H2, NGNG) 
      ENDIF 
!----EZM                                                                
!                                                                       
!   THE ROOTS ARE CALCULATED IN TREQ AND STORED IN R.                   
!                                                                       
      EPS2 = 1.0D-8 
!      ------------------------                                         
      CALL TREQ (AD, R, JJ, EPS2) 
!      ------------------------                                         
!                                                                       
!  SETTING ALL ARRAYS EQUAL TO ZERO.                                    
!                                                                       
      N2 = NGL * NGL 
      DO 9 I = 1, N2 
         A (I) = 0.D0 
         T (I) = 0.D0 
    9 B (I) = 0.D0 
      N2 = NGL * NGQ 
      DO 10 I = 1, N2 
   10 F (I) = 0.D0 
!                                                                       
!  THE BIG SUMMATION LOOP OVER ALL ROOTS.                               
!                                                                       
      EM = 1.D0 
      JJJ = JJ + JJ 
      IF (JJJ.GT.MEMR) THEN 
         WRITE (LUNFIL(2), 12346) MEMR, JJJ 
         IWARN = 1 
12346 FORMAT (/1X,'MEMORY OF R IN PRORES EXCEEDED,  CHANGE MEMR FROM',I4&
     &,' TO',I4)                                                        
      ENDIF 
      DO 100 J0 = 1, JJJ 
         RX = R (J0) 
         RR = RX * RX 
         QQ = 1.D0 / (RR + ADD) 
         XNORM = 0.25D0 / (1.D0 + AD * QQ) 
         FA = 0.5D0 * EM * RX * DSQRT (QQ) 
         FB = FA * DI 
!                                                                       
!      FIND THE SCALAR PRODUCTS AND STORE THEM IN PSX.                  
!                                                                       
         IF (J0.GT.JJ) GOTO 20 
!      ----------------------------------                               
         CALL SKALDP (LQX, RX, SI, CO, PSX, PSH) 
         GOTO 30 
   20    CALL SKALDP (LQX, RX, SI, CO, PSH, PSX) 
!      ----------------------------------                               
!                                                                       
!  SLAB CALCULATES THE Y-FUNCTIONS IN THE +1 AND -1 POINTS AND THE      
!  SCALAR PRODUCTS OF THE Y-FUNCTIONS.                                  
!                                                                       
   30    R0 = R1 * RR 
!      -----------------------------------------------------------      
         CALL SLAB (NG, LQX, HY, R0, DIF, EIG, EIG (1 + NG), VECT, VECT &
         (1 + NNG), H, H (LQ), H (LW), H (LC), H2 (1 + NGNG) )          
!      -----------------------------------------------------------      
!                                                                       
!   CALCULATION OF THE RESPONSE MATRIX ELEMENTS.                        
!                                                                       
         XF = HY * XNORM 
         K1 = 0 
         IF (IDIAG.NE.0) GOTO 51 
         DO 50 N = 1, NG 
            N1 = N + NG 
            N2 = N1 + NG 
            XB = XNORM * (FA + DIF (N) * FB) * HYX 
            L1 = - 1 
            DO 48 M = 1, NG 
               M1 = (M - 1) * NGL2 
               NM = M1 + N 
               NM1 = M1 + N1 
               NM2 = M1 + N2 
!                                                                       
               DO 45 I = 1, NL 
                  K = K1 + I 
                  X1 = XNORM * PSX (I) 
                  X3 = XB * H (NM2) 
                  DO 40 J = 1, NL 
                     KL = (L1 + J) * NGL + K 
                     X2 = X1 * PSX (J) 
                     A (KL) = A (KL) + X2 * H (NM1) 
                     T (KL) = T (KL) + X2 * H (NM) 
   40             B (KL) = B (KL) + X3 * PSX (J) 
!                                                                       
                  X4 = XF * PSX (I) 
                  NK = NQQ * (N - 1 + NG * (I - 1 + NL * (M - 1) ) ) 
                  NM3 = N + NGL2 * (M - 1) + MG 
                  DO 43 KK = 1, NQ 
                     DO 42 LL = 1, NQ 
   42                F (NK + LL) = F (NK + LL) + X4 * PSX (LL) * H (NM3) 
                     NK = NK + NQ 
   43             NM3 = NM3 + NG 
   45          NM2 = NM2 + NG 
               L1 = L1 + NL 
   48       END DO 
            K1 = K1 + NL 
   50    END DO 
!                                                                       
!----EZM                                                                
         GOTO 69 
!ALCULATE B(G,L,G',L')  USING DD(G,G') = E(G,G)                         
   51    DO 60 N = 1, NG 
            N1 = N + NG 
            N2 = N1 + NG 
            XB = XNORM * FA * HYX 
            L1 = - 1 
            DO 58 M = 1, NG 
               M1 = (M - 1) * NGL2 
               NM = M1 + N 
               NM1 = M1 + N1 
               NM2 = M1 + N2 
               DO 55 I = 1, NL 
                  K = K1 + I 
                  X1 = XNORM * PSX (I) 
                  X3 = XB * H (NM2) 
                  DO 54 J = 1, NL 
                     KL = (L1 + J) * NGL + K 
                     X2 = X1 * PSX (J) 
                     A (KL) = A (KL) + X2 * H (NM1) 
                     T (KL) = T (KL) + X2 * H (NM) 
   54             B (KL) = B (KL) + X3 * PSX (J) 
                  X4 = XF * PSX (I) 
                  NK = NQQ * (N - 1 + NG * (I - 1 + NL * (M - 1) ) ) 
                  NM3 = N + NGL2 * (M - 1) + MG 
                  DO 53 KK = 1, NQ 
                     DO 52 LL = 1, NQ 
   52                F (NK + LL) = F (NK + LL) + X4 * PSX (LL) * H (NM3) 
                     NK = NK + NQ 
   53             NM3 = NM3 + NG 
   55          NM2 = NM2 + NG 
               L1 = L1 + NL 
   58       END DO 
            K1 = K1 + NL 
   60    END DO 
!----EZM                                                                
!                                                                       
   69    EM = - EM 
         IF (J0.EQ.JJ) EM = 1.D0 
  100 END DO 
!                                                                       
!----EZM                                                                
!ALCULATE B(L,L') = DD * B(L,L')                                        
      IF (IDIAG.EQ.0) RETURN 
      I1 = 0 
      DO 64 M = 1, NG 
         DO 64 J = 1, NL 
            DO 64 N = 1, NG 
               I11 = (J - 1) * NGL + (M - 1) * NGL * NL 
               DO 62 I = 1, NL 
                  I1 = I1 + 1 
                  I11 = I11 + 1 
                  SOM = 0.D0 
                  DO 61 IG = 1, NG 
                     IDD = (IG - 1) * NG + N 
                     DD = H2 (IDD) 
!                      I2=I1 + NL*(IG-1)                                
                     I2 = I11 + NL * (IG - 1) 
                     SOM = SOM + DD * B (I2) 
   61             END DO 
                  H2 (I1 + NGNG) = SOM 
   62          END DO 
   64 CONTINUE 
      KAK = NGNG * NL * NL 
      DO 65 K = 1, KAK 
         I = NGNG + K 
         B (K) = H2 (I) 
   65 END DO 
!----EZM                                                                
      RETURN 
      END SUBROUTINE PRORES                         
!HETSMG                      26 MAI 1982                             RES
!                                                                       
!                            **********                                 
!                            * HETSMG *                                 
!                            **********                                 
!                                                                       
!  MULTIGROUP RESPONSE MATRICES ARE IMPROVED BY ADDING                  
!  1/N**K - TERMS TO RESPMG SOLUTION.                                   
!                                                                       
!      -------------------------------------------------                
      SUBROUTINE HETSMG (LGL, LGQ, HX, HY, DD0, DIF, JJ, A, B, F, H) 
!      -------------------------------------------------                
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MEMG = 10) 
      DIMENSION DIF ( * ), A (LGL, * ), B (LGL, * ), F (LGQ, * ),       &
      ALF (MEMG), GAM (MEMG), SUM (4, 2), H ( * )                       
      COMMON / CRMG / NG, NL, NQ, NQQ, NGL, NGQ, LQX, IADJNT 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / DINV / IDIAG 
      COMMON / WARN / IWARN 
      DATA PI / 3.141592654D0 / 
!                                                                       
      MCHECK = MAX (NL, LQX) 
      IF (MCHECK.GT.MEMG) THEN 
         WRITE (LUNFIL(2), 12345) MEMG, MCHECK 
         IWARN = 1 
12345 FORMAT (/1X,'LEGENDRE MOMENTS GREATER THAN',I2,' EXCEEDS MEMORY IN&
     & HETSMG',/                                                        &
     & 1X,'CHANGE MEMG TO',I3)                                          
      ENDIF 
!                                                                       
      NGNG = NG * NG 
      A0 = 0.25D0 * HX / DD0 
      DO 10 I = 1, LQX 
         X = I * (I - 1) 
         ALF (I) = 2.D0 * A0 + X 
   10 GAM (I) = 0.5D0 * HX / HY * X 
!                                                                       
!  CALCULATE THE SUMS   SUM <PI*N>**-II  WHERE  II=4,5,6  AND  N=JJ-1,JJ
!  AND            SUM <PI*(N-0.5)>**-II  WHERE  II=4,5,6  AND  N=JJ,JJ+1
!                                                                       
      DO 20 I = 1, 4 
         II = I + 3 
         I1 = II - 1 
         EJ = JJ - 0.5D0 
         SS = DBLE (I1) * PI**II 
         XX = EJ + DBLE (II) / EJ / 24.D0 
         SUM (I, 1) = 1.D0 / SS / XX**I1 
         XX = JJ + DBLE (II) / DBLE (JJ) / 24.D0 
   20 SUM (I, 2) = 1.D0 / SS / XX**I1 
!                                                                       
!  IMPROVE MATRICES A AND B                                             
!                                                                       
      NN = 0 
      MM = 0 
      IF (IDIAG.NE.0) GOTO 60 
      DO 50 N = 1, NG 
         AG = 0.25D0 * HX / DIF (N) 
         XBB = 0.5D0 * AG * (1.D0 + DIF (N) / DD0) 
         IE = 1 
         DO 40 I1 = 1, NL 
            I = NN + I1 
            JE = 1 
            DO 30 J1 = 1, NL 
               J = NN + J1 
               M = 1 
               IF (JE.LT.0) M = 2 
               IF (IE * JE) 26, 26, 24 
   24          XA = AG * ALF (I1) * ALF (J1) 
               A (I, J) = A (I, J) + XA * (SUM (2, M) - AG * SUM (3, M) &
               )                                                        
   26          XB = DBLE (IE) * XBB * ALF (J1) 
               X2 = GAM (I1) + AG 
               B (I, J) = B (I, J) + XB * (SUM (1, M) - X2 * SUM (2, M) &
               )                                                        
   30       JE = - JE 
!                                                                       
!  IMPROVE MATRIX F                                                     
!                                                                       
            IL = 1 
            DO 38 L = 1, NQ 
               IK = 1 
               DO 36 K = 1, NQ 
                  IF (IK * IE) 36, 36, 34 
   34             XF = DBLE (IL) * HX * AG * ALF (K) * ALF (I1) 
                  X2 = GAM (L) + AG 
                  INGQ = MM + K + (L - 1) * NQ 
                  F (INGQ, I) = F (INGQ, I) + XF * (SUM (3, M) - X2 *   &
                  SUM (4, M) )                                          
   36          IK = - IK 
   38       IL = - IL 
   40    IE = - IE 
         MM = MM + NQQ 
   50 NN = NN + NL 
      RETURN 
!----EZM                                                                
!ALCULATE DD = DINV + E/DD0                            -> H             
   60 IAD = 6 * NG + NGNG + 1 
      IM1 = 1 + NGNG 
      DDX = 1.D0 / DD0 
      CALL SETID (H (IM1), NG) 
      CALL MULTA (H (IM1), DDX, NGNG) 
      CALL SUMAB (H (IM1), DIF (IAD), H, NGNG) 
      NJJ = 0 
      NNI = 0 
      DO 160 MMN = 1, NG 
         MM = 0 
         NII = 0 
         DO 150 N = 1, NG 
            IAD = 6 * NG + NGNG + (MMN - 1) * NG + N 
            IAJ = (MMN - 1) * NG + N 
            AG = 0.25D0 * HX * DIF (IAD) 
            XBB = 0.5D0 * 0.25D0 * HX * H (IAJ) 
            IE = 1 
            DO 140 I1 = 1, NL 
               I = NII + I1 
               II = NNI + I1 
               JE = 1 
               DO 130 J1 = 1, NL 
                  J = NJJ + J1 
                  M = 1 
                  IF (JE.LT.0) M = 2 
                  IF (IE * JE) 126, 126, 124 
  124             XA = AG * ALF (I1) * ALF (J1) 
                  A (I, J) = A (I, J) + XA * (SUM (2, M) - AG * SUM (3, &
                  M) )                                                  
  126             XB = DBLE (IE) * XBB * ALF (J1) 
                  X2 = GAM (I1) + AG 
                  B (I, J) = B (I, J) + XB * (SUM (1, M) - X2 * SUM (2, &
                  M) )                                                  
  130          JE = - JE 
!                                                                       
!  IMPROVE MATRIX F                                                     
!                                                                       
               IL = 1 
               DO 138 L = 1, NQ 
                  IK = 1 
                  DO 136 K = 1, NQ 
                     IF (IK * IE) 136, 136, 134 
  134                XF = DBLE (IL) * HX * AG * ALF (K) * ALF (I1) 
                     X2 = GAM (L) + AG 
                     INGQ = MM + K + (L - 1) * NQ 
                     F (INGQ, II) = F (INGQ, II) + XF * (SUM (3, M)     &
                     - X2 * SUM (4, M) )                                
  136             IK = - IK 
  138          IL = - IL 
  140       IE = - IE 
            MM = MM + NQQ 
  150    NII = NII + NL 
         NNI = NNI + NL 
  160 NJJ = NJJ + NL 
!                                                                       
      RETURN 
      END SUBROUTINE HETSMG                         
!BALPRO                      31 OCT 1983                             RES
!                                                                       
!                            **********                                 
!                            * BALPRO *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE BALPRO (LG, LGL, LGQ, DD0, XSEC, XK, A, T, B, F, H, H2) 
!                                                                       
!  BALPRO CHECKS THE NEUTRON CONSERVATION OF THE PROVISIONAL RESPONSE MA
!  CONSERVATION IS FORCED BY CORRECTING THE SIDE-TRANSMISSION MATRIX.   
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION XSEC (LG, * ), A (LGL, * ), T (LGL, * ), B (LGL,        &
      * ), F (LGQ, * ), H ( * ), H2 (LGL, * )                           
      COMMON / CRMG / NG, NL, NQ, NQQ, NGL, NGQ, LQX, IADJNT 
      COMMON / DINV / IDIAG 
!                                                                       
      NGNG = NG * NG 
      IA2 = 2 
      IA3 = 3 
      IA4 = 4 
      IF (IDIAG.NE.0) GOTO 30 
      DO 20 N1 = 1, NGL, NL 
         N = (N1 - 1) / NL + 1 
         DDD = 4.D0 * XSEC (N, 1) / (DD0 + XSEC (N, 1) ) 
         SIGREM = 0.D0 
         DO 8 I = 1, NG 
            IF (IADJNT.NE.1) SIGREM = SIGREM + XSEC (I, N + 4) 
    8    IF (IADJNT.EQ.1) SIGREM = SIGREM + XSEC (N, I + 4) 
         ABN = XSEC (N, IA2) + SIGREM 
         DO 18 JM = 1, NL, 2 
            DO 18 M1 = JM, NGL, NL 
               SUM = 0.D0 
               I1 = 1 
               DO 10 I = 1, NG 
!                                                                       
!      ORDINARY SOLUTION :                                              
                  IA = I 
                  NA = N 
                  IF (IADJNT.NE.1) GOTO 9 
!      ADJOINT SOLUTION :                                               
                  IA = N 
                  NA = I 
!                                                                       
    9             SUM = SUM + (XSEC (NA, IA4) / XK * XSEC (IA, IA3)     &
                  + XSEC (NA, IA + 4) ) * F (I1, M1)                    
   10          I1 = I1 + NQQ 
               NF = (N - 1) * NQQ + 1 
               ERR = A (N1, M1) + T (N1, M1) + DDD * B (N1, M1) + ABN * &
               F (NF, M1) - SUM                                         
               IF (M1.EQ. (N - 1) * NL + 1) ERR = ERR - 2.D0 
               B (N1, M1) = B (N1, M1) - ERR / DDD 
   18    CONTINUE 
   20 END DO 
      RETURN 
!----EZM                                                                
!ALCULATE DD = 4/(DD0*DINV + E)                         -> H            
   30 IM1 = 1 + NGNG 
      IM2 = IM1 + NGNG 
      CALL EQUATE (H (IM2), XSEC (1, 7 + NG), NGNG) 
      CALL MULTA (H (IM2), DD0, NGNG) 
      CALL SETID (H (IM1), NG) 
      CALL SUMAB (H (IM2), H (IM1), H, NGNG) 
      CALL MINV (H, NG, DET, H (IM1), H (IM2) ) 
      CALL MULTA (H, 4.D0, NGNG) 
      CALL SETVAL (H2, 0.D0, NGL * NGL) 
!ALCULATE H2(N1,M1) = DD*B(N1,M1) - ERR                  -> H2          
      DO 40 N1 = 1, NGL, NL 
         N = (N1 - 1) / NL + 1 
         SIGREM = 0.D0 
         DO 48 I = 1, NG 
            IF (IADJNT.NE.1) SIGREM = SIGREM + XSEC (I, N + 4) 
   48    IF (IADJNT.EQ.1) SIGREM = SIGREM + XSEC (N, I + 4) 
         ABN = XSEC (N, IA2) + SIGREM 
         DO 58 JM = 1, NL, 2 
            DO 58 M1 = JM, NGL, NL 
               SUM = 0.D0 
               I1 = 1 
               DO 60 I = 1, NG 
!                                                                       
!      ORDINARY SOLUTION :                                              
                  IA = I 
                  NA = N 
                  IF (IADJNT.NE.1) GOTO 59 
!      ADJOINT SOLUTION :                                               
                  IA = N 
                  NA = I 
!                                                                       
   59             SUM = SUM + (XSEC (NA, IA4) / XK * XSEC (IA, IA3)     &
                  + XSEC (NA, IA + 4) ) * F (I1, M1)                    
   60          I1 = I1 + NQQ 
               NF = (N - 1) * NQQ + 1 
               SOM = 0.D0 
               DO 70 II = 1, NG 
                  IADD = (II - 1) * NG + N 
                  NDD = (II - 1) * NL + 1 
                  DDD = H (IADD) 
                  SOM = SOM + DDD * B (NDD, M1) 
   70          END DO 
               H2 (N1, M1) = SOM 
               ERR = A (N1, M1) + T (N1, M1) + H2 (N1, M1) + ABN * F (  &
               NF, M1) - SUM                                            
!               ERR=A(N1,M1)+T(N1,M1)+DDD*B(N1,M1)+ABN*F(NF,M1)-SUM     
               IF (M1.EQ. (N - 1) * NL + 1) ERR = ERR - 2.D0 
               H2 (N1, M1) = H2 (N1, M1) - ERR 
!               B(N1,M1)=B(N1,M1) - ERR/DDD                             
   58    CONTINUE 
   40 END DO 
!ALCULATE B(N1,M1) = (1/DD)*H2(N1,M1)                      -> B         
      CALL MINV (H, NG, DET, H (IM1), H (IM2) ) 
      DO 140 N1 = 1, NGL, NL 
         N = (N1 - 1) / NL + 1 
         DO 158 JM = 1, NL, 2 
            DO 158 M1 = JM, NGL, NL 
               SOM = 0.D0 
               DO 170 II = 1, NG 
                  IADD = (II - 1) * NG + N 
                  NDD = (II - 1) * NL + 1 
                  DDD = H (IADD) 
                  SOM = SOM + DDD * H2 (NDD, M1) 
  170          END DO 
               B (N1, M1) = SOM 
  158    CONTINUE 
  140 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE BALPRO                         
!SLAB                        26 MAI 1982                             RES
!                                                                       
!                            ********                                   
!                            * SLAB *                                   
!                            ********                                   
!                                                                       
!   SLAB CALCULATES MULTIGROUP SLAB RESPONSE MATRICES.                  
!   ALBEDO (A), TRANSMISSION (T) AND FLUX (F) MATRICES ARE STORED IN ARR
!          H = TRANSPOSE ( TT; AA; F0; F1; ...; FL )                    
!   WHERE  T=0.5*TT,  A=0.5*AA - I   (I = IDENTITY MATRIX)              
!          FL = L:TH MOMEMT OF FLUX MATRIX  (LEGENDRE POLYNOMIALS USED).
!                                                                       
!   DIMENSIONING OF ARRAYS (NG2=NG*NG) : DIF(NG/NG2) , C(2*NG2),        
!      VR(NG2) , VI(NG2) , Q(4*NG2) , W(2*NL*NG2) , H(MAX((4*NG2),((NL+2
!      ER(NG)  , EI(NG)                                                 
!                                                                       
      SUBROUTINE SLAB (NG, ML, HY, R0, DIF, ER, EI, VR, VI, H, Q, W, C, &
      HELP)                                                             
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MEMR = 10, MEMS = 12) 
      DIMENSION DIF ( * ), C ( * ), VR ( * ), VI ( * ), Q ( * ),        &
      W (ML, * ), H ( * ), ER ( * ), EI ( * ), REPS (MEMR), REMS (MEMR),&
      YMPS (MEMR), YMMS (MEMR), SUM (MEMS), HELP ( * )                  
      COMMON / DINV / IDIAG 
      COMMON / LUNSCM / LUNFIL (4) 
      COMMON / WARN / IWARN 
!                                                                       
      NL = ML 
      MG = 2 * NG 
      NNG = NG * NG 
      NNG2 = 2 * NNG 
      H4 = 4.D0 / HY 
      HYY = HY * HY 
      NCHECK = NL + 2 
      IF (NL.GT.MEMR.OR.NCHECK.GE.MEMS) THEN 
         WRITE (LUNFIL(2), 12345) MEMR, NL, NCHECK 
         IWARN = 1 
12345 FORMAT (/1X,'LEGENDRE MOMENTS GREATER THAN',I2,' EXCEEDS ALLOCATED&
     & MEMORY IN SLAB',/                                                &
     & 1X,'CHANGE MEMR TO',I3,' AND MEMS TO',I3)                        
      ENDIF 
!                                                                       
!----EZM                                                                
      IF (IDIAG.NE.0) THEN 
         IPLACE = 6 * NG + NG * NG + 1 
         CALL MINV (DIF (IPLACE), NG, DET, HELP, HELP (1 + NG) ) 
      ENDIF 
!----EZM                                                                
!                                                                       
!  SWEEP OVER THE NG EIGENVALUES.                                       
!                                                                       
      NEX = 0 
      J = 1 
    2 K1 = (J - 1) * NG 
      I1 = 2 * K1 
      ERR = ER (J) * HYY + R0 
      EII = EI (J) * HYY 
      ICOM = 3 
! ASSUME IMAGINARY PARTS LESS THAN 1D-10 AS ZERO ?????????              
      IF (DABS (EII) .GT.1.D-10) GOTO 4 
      ICOM = 2 
      IF (ERR.LT.0.D0) GOTO 4 
      ICOM = 1 
!                                                                       
!  CALCULATION OF THE SQUARE ROOTS OF THE EIGENVALUES.                  
!                                                                       
    4 IF (ICOM - 2) 6, 6, 8 
    6 ALF = DSQRT (DABS (ERR) ) 
      YMP = 0.D0 
      YMM = 0.D0 
      DIMP = 0.D0 
      DIMM = 0.D0 
      DO 7 LP = 1, NL 
         YMPS (LP) = 0.D0 
    7 YMMS (LP) = 0.D0 
      GOTO 14 
    8 IF (ERR.LE.0.D0) GOTO 10 
      X = (EII * EII) / (ERR * ERR) 
      IF (X.GT.0.01D0) GOTO 10 
      BET = DSQRT (0.25D0 * ERR * X * (1.D0 - 0.25D0 * X * (1.D0 -      &
      0.5D0 * X * (1.D0 - 0.625D0 * X) ) ) )                            
      GOTO 12 
   10 BET = DSQRT (0.5D0 * (DSQRT (ERR * ERR + EII * EII) - ERR) ) 
   12 ALF = 0.5D0 * EII / BET 
!                                                                       
!  WAVE EQUATION SOLUTIONS AND THEIR DERIVATIVES AT THE POINT +1.       
!  THE INTEGRALS OF THE WAVE EQUATION SOLUTIONS.                        
!                                                                       
   14 IF (ICOM - 2) 16, 18, 20 
!                                                                       
!  ICOM=1.  REAL AND POSITIVE EIGENVALUES.                              
!                                                                       
   16 CALL SKALEX (NL, ALF, EX, REPS, REMS) 
      IF (ALF.GE.9.D0) NEX = NEX + 1 
      REP = 1.D0 
      REM = EX 
      DREP = ALF 
      DREM = - ALF * REM 
      GOTO 22 
!                                                                       
!  ICOM=2.  REAL AND NEGATIVE EIGENVALUES.                              
!                                                                       
   18 CALL SKALDP (NL, ALF, SI, CO, REPS, REMS) 
      REP = CO 
      REM = SI 
      DREP = - ALF * SI 
      DREM = ALF * CO 
      GOTO 22 
!                                                                       
!  ICOM=3.  COMPLEX EIGENVALUES.                                        
!                                                                       
   20 CALL SKALHE (NL, ALF, BET, CO, SI, EX, REPS, REMS, YMPS, YMMS) 
      REP = CO 
      YMP = SI 
      REM = EX * CO 
      YMM = - EX * SI 
!EZM                                                                    
      Q1 = - ALF * CO + BET * SI 
      Q2 = - ALF * SI - BET * CO 
      DREP = - Q1 
      DIMP = - Q2 
!EZM   DREP= ALF*REP - BET*YMP                                          
!EZM   DIMP= ALF*YMP + BET*REP                                          
!EZM                                                                    
      DREM = - ALF * REM + BET * YMM 
      DIMM = - ALF * YMM - BET * REM 
   22 NIX = - 1 
   21 NIX = NIX + 1 
!                                                                       
!  CALCULATION OF THE H MATRIX USED FOR CALCULATION OF THE COEFFICIENTS 
!  OF THE GENERAL SOLUTION, THE Q MATRIX USED FOR CALCULATION OF THE TRA
!  MISSION AND ALBEDO MATRICES AND THE W-L MATRICES USED FOR CALCULATION
!  OF THE P-L FLUX MOMENTS MATRICES.                                    
!                                                                       
      KK1 = K1 
      DO 42 I = 1, NG 
         I1 = I1 + 1 
         K1 = K1 + 1 
         I2 = I1 + NNG2 
         K2 = K1 + NNG 
         X = VR (K1) 
         Y = VI (K1) 
         Q (I1) = X * REP - Y * YMP 
         Q (I2) = X * REM - Y * YMM 
         IF (IDIAG.NE.0) THEN 
            SUM1 = 0.D0 
            SUM2 = 0.D0 
            DO 25 IJ = 1, NG 
               IADD = (IJ - 1) * NG + I + 6 * NG + NG * NG 
               IADV = KK1 + IJ 
               D = H4 * DIF (IADD) 
               XX = VR (IADV) 
               YY = VI (IADV) 
               SUM1 = SUM1 + D * (XX * DREP - YY * DIMP) 
               SUM2 = SUM2 + D * (XX * DREM - YY * DIMM) 
   25       END DO 
         ELSE 
            D = H4 * DIF (I) 
            SUM1 = D * (X * DREP - Y * DIMP) 
            SUM2 = D * (X * DREM - Y * DIMM) 
         ENDIF 
         H (I1) = Q (I1) + SUM1 
         H (I2) = Q (I2) + SUM2 
!      H(I1)=Q(I1) + D*(X*DREP - Y*DIMP)                                
!      H(I2)=Q(I2) + D*(X*DREM - Y*DIMM)                                
         IF (NIX) 23, 23, 31 
   23    DO 28 LP = 1, NL 
            W (LP, K1) = X * REPS (LP) - Y * YMPS (LP) 
   28    W (LP, K2) = X * REMS (LP) - Y * YMMS (LP) 
   31    IF (ICOM - 2) 42, 42, 32 
   32    I3 = I1 + MG 
         I4 = I2 + MG 
         K3 = K1 + NG 
         K4 = K2 + NG 
         Q (I3) = X * YMP + Y * REP 
         Q (I4) = X * YMM + Y * REM 
         IF (IDIAG.NE.0) THEN 
            SUM3 = 0.D0 
            SUM4 = 0.D0 
            DO 24 IJ = 1, NG 
               IADD = (IJ - 1) * NG + I + 6 * NG + NG * NG 
               IADV = KK1 + IJ 
               D = H4 * DIF (IADD) 
               XX = VR (IADV) 
               YY = VI (IADV) 
               SUM3 = SUM3 + D * (XX * DIMP + YY * DREP) 
               SUM4 = SUM4 + D * (XX * DIMM + YY * DREM) 
   24       END DO 
         ELSE 
            D = H4 * DIF (I) 
            SUM3 = D * (X * DIMP + Y * DREP) 
            SUM4 = D * (X * DIMM + Y * DREM) 
         ENDIF 
         H (I3) = Q (I3) + SUM3 
         H (I4) = Q (I4) + SUM4 
!      H(I3)=Q(I3) + D*(X*DIMP + Y*DREP)                                
!      H(I4)=Q(I4) + D*(X*DIMM + Y*DREM)                                
         IF (NIX) 33, 33, 42 
   33    DO 36 LP = 1, NL 
            W (LP, K3) = X * YMPS (LP) + Y * REPS (LP) 
   36    W (LP, K4) = X * YMMS (LP) + Y * REMS (LP) 
   42 END DO 
      H4 = - H4 
      IF (NIX) 43, 43, 54 
!                                                                       
!  THE WAVE EQUATION SOLUTIONS AND THEIR DERIVATIVES IN THE POINT  -1.  
!                                                                       
   43 IF (ICOM - 2) 44, 46, 48 
   44 REP = EX 
      REM = 1.D0 
      DREP = - DREM 
      DREM = - ALF 
      GOTO 50 
   46 REM = - SI 
      DREP = - DREP 
      GOTO 50 
   48 REP = REM 
      REM = CO 
      YMP = YMM 
      YMM = SI 
      DREP = - DREM 
      DIMP = - DIMM 
!EZM   DREM=-ALF*REM + BET*YMM                                          
!EZM   DIMM=-ALF*YMM - BET*REM                                          
      DREM = Q1 
      DIMM = Q2 
   50 K1 = K1 - NG 
      GOTO 21 
   54 J = J + 1 
      IF (ICOM.EQ.3) J = J + 1 
      IF (J.LE.NG) GOTO 2 
!  END OF THE SWEEP OVER THE EIGENVALUES.                               
!                                                                       
!  CALCULATION OF THE COEFFICIENTS OF THE GENERAL SOLUTION.             
!                                                                       
      IF (NG.NE.2.OR.NEX.NE.2) GOTO 55 
      DET = 4.D0 / (H (11) * H (16) - H (12) * H (15) ) 
      C (1) = 0.D0 
      C (2) = 0.D0 
      C (5) = 0.D0 
      C (6) = 0.D0 
      C (3) = DET * H (16) 
      C (4) = - DET * H (12) 
      C (7) = - DET * H (15) 
      C (8) = DET * H (11) 
      GOTO 57 
   55 MGS = 1 
      IF (NG.EQ.NEX) MGS = NG + 1 
      IS = 1 
      DO 56 J = 1, NG 
         DO 56 I = MGS, MG 
            C (IS) = 0.D0 
            IF ( (I - NG) .EQ.J) C (IS) = 4.D0 
   56 IS = IS + 1 
      IF (NG.NE.NEX) GOTO 256 
      NA = NG 
      IS = 0 
      DO 156 I = 1, NG 
         IT = MG * (NG + I - 1) + NG 
         DO 156 J = 1, NG 
            IS = IS + 1 
            IT = IT + 1 
  156 H (IS) = H (IT) 
      GOTO 356 
  256 NA = MG 
  356 EPS = 1.D-8 
!                                                                       
!   CALL THE STANDARD SUBROUTINE GELG THAT SOLVES NG SETS OF LINEAR EQUA
!   OF ORDER NA SIMULTANEOUSLY.                                         
!                                                                       
      CALL GELG (C, H, NA, NG, EPS, IER) 
      IF (IER.EQ.1) THEN 
         WRITE (LUNFIL(2), 290) IER, HY, R0, (DIF (II), II = 1, NG) 
!        WRITE(6,291) 'EIGENVALUES  - RE', (ER(II),II=1,NG)             
!        WRITE(6,291) 'EIGENVALUES  - IM', (EI(II),II=1,NG)             
!        WRITE(6,291) 'EIGENVECT.   - RE', (VR(II),II=1,NNG)            
!        WRITE(6,291) 'EIGENVECT.   - IM', (VI(II),II=1,NNG)            
         WRITE (LUNFIL(2), 291) (ER (II), II = 1, NG) 
         WRITE (LUNFIL(2), 292) (EI (II), II = 1, NG) 
         WRITE (LUNFIL(2), 293) (VR (II), II = 1, NNG) 
         WRITE (LUNFIL(2), 294) (VI (II), II = 1, NNG) 
  290 FORMAT    (//1X,75(1H*)/' ERROR IN GELG, IER=',I5//               &
     &          'HY =',F6.2,4X,'R0 =',F7.3,4X,'DIFF =',6F6.3)           
  291 FORMAT    (1X,'EIGENVALUES - RE   ',4X,10F10.6,/(22X,10F10.6)) 
  292 FORMAT    (1X,'EIGENVALUES - IM   ',4X,10F10.6,/(22X,10F10.6)) 
  293 FORMAT    (1X,'EIGENVECT.  - RE   ',4X,10F10.6,/(22X,10F10.6)) 
  294 FORMAT    (1X,'EIGENVECT.  - IM   ',4X,10F10.6,/(22X,10F10.6)) 
         STOP 'LABAN TRAGICALLY DECEASED DUE TO ERROR IN GELG' 
      ENDIF 
!                                                                       
  300 IF (NG.NE.NEX) GOTO 57 
      IS = NG * NG + 1 
      DO 157 I = 1, NG 
         IT = MG * (NG - I) + MG + 1 
         DO 157 J = 1, NG 
            IS = IS - 1 
            IT = IT - 1 
            C (IT) = C (IS) 
            IU = IT - NG 
  157 C (IU) = 0.D0 
!                                                                       
!  THE TRANSMISSION, ALBEDO AND P-0,P-1, .... P-L FLUX MATRICES ARE     
!  CALCULATED AND STORED IN THE ARRAY H.                                
!                                                                       
   57 NX = NL + 2 
      DO 72 I = 1, NG 
         IS = 0 
         IF0 = 0 
         DO 72 J = 1, NG 
            DO 58 K = 1, NX 
   58       SUM (K) = 0.D0 
            ING = 0 
            DO 68 K = 1, MG 
               IS = IS + 1 
               CC = C (IS) 
               IG = ING + I 
               IT = 2 * ING + I 
               IA = IT + NG 
               SUM (1) = SUM (1) + Q (IT) * CC 
               SUM (2) = SUM (2) + Q (IA) * CC 
               DO 60 LP = 1, NL 
   60          SUM (LP + 2) = SUM (LP + 2) + W (LP, IG) * CC 
   68       ING = ING + NG 
            IF = IF0 + I 
            DO 70 K = 1, NX 
               H (IF) = SUM (K) 
   70       IF = IF + NG 
   72 IF0 = IF0 + NG * NX 
!                                                                       
!----EZM                                                                
      IF (IDIAG.NE.0) THEN 
         IPLACE = 6 * NG + NG * NG + 1 
         CALL MINV (DIF (IPLACE), NG, DET, HELP, HELP (1 + NG) ) 
      ENDIF 
!----EZM                                                                
      RETURN 
      END SUBROUTINE SLAB                           
!TREQ                        26 MAI 1982                             RES
!                                                                       
!                            ********                                   
!                            * TREQ *                                   
!                            ********                                   
!                                                                       
      SUBROUTINE TREQ (P, R, N, EPS) 
!                                                                       
!         TREQ SOLVES THE TRANSCENDENTAL EQUATIONS AND STORES THE ROOTS 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION R ( * ) 
      DATA PI, PI2 / 3.141592654D0, 1.570796327D0 / 
!                                                                       
      JMAX = 8 
      A1 = 1.D0 / 3.D0 + 1.D0 / P 
      A2 = 0.2D0 + 4.D0 / 3.D0 / P + 2.D0 / (P * P) 
      S = P + 1.D0 
      B1 = P / 3.D0 / S 
      B2 = P * (3.D0 * P - 2.D0) / 15.D0 / (S * S) 
      IF (P.LE.1.D0) THEN 
         Q1 = DSQRT (P * (1.D0 - P * (1.D0 / 3.D0 - P * (4.D0 / 45.D0 - &
         0.016931D0 * P) ) ) )                                          
!      Q1=DSQRT(P*(1.D0-P*(0.333333D0-P*(0.088889D0-0.016931D0*P))))    
      ENDIF 
!                                                                       
      PIS = DBLE (N - 1) * PI 
      NN = 0 
    2 I = N 
    4 J = 0 
      I2 = I + NN 
!                                                                       
      IF (I2.GT.1) GOTO 101 
      IF (P.GT.1.D0) GOTO 103 
      X1 = Q1 
      GOTO 11 
!                                                                       
  101 Z = P / PIS 
      IF (Z.GT.0.7D0) GOTO 102 
      ZZ = Z * Z 
      X1 = Z * (1.D0 - ZZ * (A1 - A2 * ZZ) ) 
      GOTO 11 
!                                                                       
  102 IF (Z.GT.1.3D0) GOTO 103 
      X1 = (1.2337D0 + 0.5708D0 * PIS + P) / (2.5708D0 + 2.D0 * PIS) 
      GOTO 11 
!                                                                       
  103 Z = (PIS + PI2) / S 
      ZZ = Z * Z 
      X1 = PI2 - Z * (1.D0 - ZZ * (B1 - B2 * ZZ) ) 
!                                                                       
   11 X11 = PIS + X1 
!                                                                       
      CO = DCOS (X1) 
      U = CO * CO 
      DTAN = DSIN (X1) / CO 
      X2 = X1 - (X11 * DTAN - P) / (DTAN + X11 / U) 
      IF (X2.GT.PI2) X2 = DATAN (P / (PIS + PI2) ) 
      IF (DABS (X2 - X1) .LT.EPS) GOTO 16 
!                                                                       
      J = J + 1 
      IF (J.GT.JMAX) GOTO 16 
      X1 = X2 
      GOTO 11 
   16 R (I2) = X2 + PIS 
      I = I - 1 
      PIS = PIS - PI 
      IF (I - 1) 18, 4, 4 
   18 IF (NN.EQ.N) RETURN 
      PIS = DBLE (N - 1) * PI + PI2 
      NN = N 
      GOTO 2 
!                                                                       
      END SUBROUTINE TREQ                           
!SKALDP                      26 MAI 1982                             RES
!                                                                       
!                            **********                                 
!                            * SKALDP *                                 
!                            **********                                 
      SUBROUTINE SKALDP (NL, ALF, SI, CO, REPS, REMS) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MEMP = 20) 
      DIMENSION REPS ( * ), REMS ( * ), PS (MEMP), MS (MEMP) 
      DOUBLEPRECISION PS, MS, A, AI, X, SII, COO 
      COMMON / WARN / IWARN 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      MHECK = NL * 2 
      IF (MHECK.GT.MEMP) THEN 
         WRITE (LUNFIL(2), 12345) NL, MHECK 
         IWARN = 1 
12345 FORMAT (/1X,'LEGENDRE MOMENTS GREATER THAN',I2,' EXCEEDS ALLOCATED&
     & MEMORY IN SKALDP',/                                              &
     & 1X,'CHANGE MEMP TO',I3)                                          
      ENDIF 
!         SCALAR PRODUCTS FOR NEGATIVE EIGENVALUES.                     
      A = ALF 
      AI = 1.D+0 / A 
!      SII=DSIN(A)                                                      
!      COO=DCOS(A)                                                      
      SII = DSIN (ALF) 
      COO = DCOS (ALF) 
      PS (1) = 2.D+0 * AI * SII 
      MS (1) = 0.D+0 
      PS (3) = 0.D+0 
      MS (3) = AI * PS (1) - 2.D+0 * AI * COO 
      IF (NL - 3) 12, 8, 8 
    8 NL2 = 2 * NL 
      DO 10 I = 5, NL2, 2 
         X = I - 2 
         X = X * AI 
         PS (I) = PS (I - 4) - X * MS (I - 2) 
         MS (I) = MS (I - 4) + X * PS (I - 2) 
         IF (DABS (PS (I) ) .GT.1.D+10) PS (I) = 0.0D+0 
         IF (DABS (MS (I) ) .GT.1.D+10) MS (I) = 0.0D+0 
   10 END DO 
   12 DO 14 I = 1, NL 
         ID = 2 * I - 1 
         REPS (I) = PS (ID) 
   14 REMS (I) = MS (ID) 
      SI = SII 
      CO = COO 
!                                                                       
      RETURN 
      END SUBROUTINE SKALDP                         
!SKALEX                      26 MAI 1982                             RES
!                                                                       
!                            **********                                 
!                            * SKALEX *                                 
!                            **********                                 
!                                                                       
!         SCALAR PRODUCTS FOR POSITIVE EIGENVALUES.                     
!                                                                       
      SUBROUTINE SKALEX (NL, ALF, EX, REPS, REMS) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MEMP = 20) 
      DIMENSION REPS ( * ), REMS ( * ), PS (MEMP) 
      DOUBLEPRECISION PS, A, AI, X, EXX 
      COMMON / WARN / IWARN 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      MHECK = NL * 2 
      IF (MHECK.GT.MEMP) THEN 
         WRITE (LUNFIL(2), 12345) NL, MHECK 
         IWARN = 1 
12345 FORMAT (/1X,'LEGENDRE MOMENTS GREATER THAN',I2,' EXCEEDS ALLOCATED&
     & MEMORY IN SKALEX',/                                              &
     & 1X,'CHANGE MEMP TO',I3)                                          
      ENDIF 
      A = ALF 
      AI = 1.D+0 / A 
      EXX = 0.D+0 
!      IF (A.LT.16.D+0)  EXX=DEXP(-2.D+0*A)                             
      IF (A.LT.16.D+0) EXX = DEXP ( - 2.D+0 * ALF) 
      PS (1) = AI * (1.D+0 - EXX) 
      PS (3) = AI * (1.D+0 + EXX - PS (1) ) 
      IF (NL - 3) 22, 18, 18 
   18 NL2 = 2 * NL 
      DO 20 I = 5, NL2, 2 
         X = I - 2 
   20 PS (I) = PS (I - 4) - X * AI * PS (I - 2) 
   22 EI = 1.0 
      DO 24 I = 1, NL 
         ID = 2 * I - 1 
         REPS (I) = PS (ID) 
         REMS (I) = EI * REPS (I) 
   24 EI = - EI 
      EX = EXX 
!                                                                       
      RETURN 
      END SUBROUTINE SKALEX                         
!SKALHE                      26 MAI 1982                             RES
!                                                                       
!                            **********                                 
!                            * SKALHE *                                 
!                            **********                                 
!                                                                       
!         SCALAR PRODUCTS FOR COMPLEX EIGENVALUES.                      
!                                                                       
      SUBROUTINE SKALHE (NL, ALF, BET, CO, SI, EX, REPS, REMS, YMPS,    &
      YMMS)                                                             
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MEMP = 20) 
      DIMENSION REPS ( * ), REMS ( * ), YMPS ( * ), YMMS ( * ), PS (    &
      MEMP), PY (MEMP)                                                  
      DOUBLEPRECISION PS, PY, A, B, X, AP, C1, C2, C3, C4, SII, COO,    &
      EXX                                                               
      COMMON / WARN / IWARN 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
      MHECK = NL * 2 
      IF (MHECK.GT.MEMP) THEN 
         WRITE (LUNFIL(2), 12345) NL, MHECK 
         IWARN = 1 
12345 FORMAT (/1X,'LEGENDRE MOMENTS GREATER THAN',I2,' EXCEEDS ALLOCATED&
     & MEMORY IN SKALHE',/                                              &
     & 1X,'CHANGE MEMP TO',I3)                                          
      ENDIF 
      A = ALF 
      B = BET 
      AP = 1.D+0 / (A * A + B * B) 
!      COO=DCOS(B)                                                      
!      SII=DSIN(B)                                                      
      COO = DCOS (BET) 
      SII = DSIN (BET) 
      EXX = 0.D+0 
!      IF (A.LT.16.D+0)  EXX=DEXP(-2.D+0*A)                             
      IF (A.LT.16.D+0) EXX = DEXP ( - 2.D+0 * ALF) 
      C1 = A * COO + B * SII 
      C2 = A * COO - B * SII 
      C3 = A * SII - B * COO 
      C4 = A * SII + B * COO 
      PS (1) = AP * (C1 - EXX * C2) 
      PY (1) = AP * (C3 + EXX * C4) 
      PS (3) = AP * (C1 + EXX * C2 - A * PS (1) - B * PY (1) ) 
      PY (3) = AP * (C3 - EXX * C4 - A * PY (1) + B * PS (1) ) 
      IF (NL - 3) 32, 28, 28 
   28 NL2 = 2 * NL 
      DO 30 I = 5, NL2, 2 
         X = I - 2 
         X = X * AP 
         I1 = I - 2 
         PS (I) = PS (I - 4) - X * (A * PS (I1) + B * PY (I1) ) 
   30 PY (I) = PY (I - 4) - X * (A * PY (I1) - B * PS (I1) ) 
   32 EI = 1.0 
      DO 34 I = 1, NL 
         ID = 2 * I - 1 
         REPS (I) = PS (ID) 
         REMS (I) = EI * REPS (I) 
         YMPS (I) = PY (ID) 
         YMMS (I) = EI * YMPS (I) 
   34 EI = - EI 
      EX = EXX 
      SI = SII 
      CO = COO 
!                                                                       
      RETURN 
      END SUBROUTINE SKALHE                         
!ONECAL                      26 MAI 1982                             RES
!                                                                       
!                            **********                                 
!                            * ONECAL *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE ONECAL (NG, NL, ONE) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION ONE ( * ) 
!                                                                       
!         ONECAL FINDS THE MATRIX ONE (STORED COLUMNWISE) THAT HAS THE E
!         WHERE THE REFLECTION MATRIX A HAS A NON-ZERO ENTRY AND THE ENT
!         WHERE A HAS A ZERO ENTRY.                                     
!                                                                       
      IJ = 0 
      DO 20 IG = 1, NG 
         EI = 1.D0 
         DO 20 IL = 1, NL 
            DO 10 JG = 1, NG 
               EJ = 1.D0 
               DO 10 JL = 1, NL 
                  IJ = IJ + 1 
                  ONE (IJ) = EI * EJ 
   10       EJ = - EJ 
   20 EI = - EI 
!                                                                       
      RETURN 
      END SUBROUTINE ONECAL                         
!MAF                         26 MAI 1982                             RES
!                                                                       
!                            *******                                    
!                            * MAF *                                    
!                            *******                                    
!                                                                       
      SUBROUTINE MAF (NGQ, NGL, FY, FX, W, H1, H2, F) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION F (NGQ, * ), H2 (NGQ, * ), FY (NGQ, * ), FX (NGQ,       &
      * ), W (NGL, * ), H1 (NGL, * )                                    
!                                                                       
!   MAF CALCULATES THE MATRIX  F = (FY - 2*FX*W) * H1 .                 
!                                                                       
      DO 10 I = 1, NGQ 
         DO 10 J = 1, NGL 
            SUM = 0.D0 
            DO 8 K = 1, NGL 
    8       SUM = SUM + FX (I, K) * W (K, J) 
            H2 (I, J) = FY (I, J) - 2.D0 * SUM 
   10 CONTINUE 
!                                                                       
      DO 20 I = 1, NGQ 
         DO 20 J = 1, NGL 
            SUM = 0.D0 
            DO 18 K = 1, NGL 
   18       SUM = SUM + H2 (I, K) * H1 (K, J) 
   20 F (I, J) = SUM 
!                                                                       
      RETURN 
      END SUBROUTINE MAF                            
!MAB                         26 MAI 1982                             RES
!                                                                       
!                            *******                                    
!                            * MAB *                                    
!                            *******                                    
!                                                                       
      SUBROUTINE MAB (N, X, ONE, B, A, T, W, H1, H2) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION B (N, * ), H2 (N, * ), A (N, * ), T (N, * ), W (N,      &
      * ), H1 (N, * ), ONE (N, * )                                      
!                                                                       
!   MAB CALCULATES THE MATRIX  B = (B - (A + X*T)*W) * H1 .             
!                                                                       
      Y = 0.D0 
      IF (X.LT.0.D0) Y = 1.D0 
      DO 10 I = 1, N 
         DO 10 J = 1, N 
            H2 (I, J) = 0.D0 
            IF (ONE (I, J) .EQ.Y) GOTO 10 
            SUM = 0.D0 
            DO 8 K = 1, N 
    8       SUM = SUM + (A (I, K) + X * T (I, K) ) * W (K, J) 
            H2 (I, J) = B (I, J) - SUM 
   10 CONTINUE 
      DO 20 I = 1, N 
         DO 20 J = 1, N 
            SUM = 0.D0 
            IF (ONE (I, J) .EQ.Y) GOTO 20 
            DO 18 K = 1, N 
   18       SUM = SUM + H2 (I, K) * H1 (K, J) 
   20 B (I, J) = SUM 
!                                                                       
      RETURN 
      END SUBROUTINE MAB                            
!MAAT                        26 MAI 1982                             RES
!                                                                       
!                            ********                                   
!                            * MAAT *                                   
!                            ********                                   
!                                                                       
      SUBROUTINE MAAT (N, ONE, T, A, B, BA, W, WA) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION T (N, * ), A (N, * ), B (N, * ), W (N, * ), BA (N,      &
      * ), WA (N, * ), ONE (N, * )                                      
!                                                                       
!   MAAT CALCULATES THE MATRICES  T = T - 2*(B*W + BA*WA) ,             
!                                 A = A - 2*(B*W - BA*WA) ,             
!                                                                       
      DO 10 I = 1, N 
         DO 10 J = 1, N 
            S3 = 0.D0 
            S4 = 0.D0 
            IF (ONE (I, J) ) 6, 9, 6 
    6       S1 = 0.D0 
            S2 = 0.D0 
            DO 8 K = 1, N 
               S1 = S1 + B (I, K) * W (K, J) 
    8       S2 = S2 + BA (I, K) * WA (K, J) 
            S3 = T (I, J) - 2.D0 * (S1 + S2) 
            S4 = A (I, J) - 2.D0 * (S1 - S2) 
    9       T (I, J) = S3 
   10 A (I, J) = S4 
!                                                                       
      RETURN 
      END SUBROUTINE MAAT                           
!IPINV                       26 MAI 1982                             RES
!                                                                       
!                            *********                                  
!                            * IPINV *                                  
!                            *********                                  
!                                                                       
      SUBROUTINE IPINV (H1, H2, WX, WY, N, NG, NL, X) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (MEMK = 50) 
      DIMENSION H1 ( * ), H2 ( * ), WX (N, * ), WY (N, * ), KN (MEMK),  &
      KM (MEMK)                                                         
      COMMON / WARN / IWARN 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
!   IPINV CALCULATES THE MATRIX  H1 = INVERS (I + X*WX*WY) .            
!                                                                       
      IH1 = 0 
      IH2 = 0 
      J = 0 
      DO 10 JG = 1, NG 
         JE = 1 
         DO 10 JL = 1, NL 
            J = J + 1 
            I = 0 
            DO 9 IG = 1, NG 
               IE = 1 
               DO 9 IL = 1, NL 
                  I = I + 1 
                  IF (IE * JE) 9, 9, 1 
    1             SUM = 0.D0 
                  DO 2 K = 1, N 
    2             SUM = SUM + WX (I, K) * WY (K, J) 
                  SUM = X * SUM 
                  IF (I.EQ.J) SUM = 1.D0 + SUM 
                  IF (JE) 3, 3, 4 
    3             IH2 = IH2 + 1 
                  H2 (IH2) = SUM 
                  GOTO 9 
    4             IH1 = IH1 + 1 
                  H1 (IH1) = SUM 
    9       IE = - IE 
   10 JE = - JE 
!                                                                       
!   CALL THE STANDARD SUBROUTINE MINV FOR INVERSION OF H1.              
!                                                                       
      NH1 = NG * ( (NL + 1) / 2) 
      NH2 = N - NH1 
      MHECK = MAX (NH1, NH2) 
      IF (MHECK.GT.MEMK) THEN 
         WRITE (LUNFIL(2), 12345) NL, MHECK 
         IWARN = 1 
12345 FORMAT (/1X,'LEGENDRE MOMENTS GREATER THAN',I2,' EXCEEDS ALLOCATED&
     & MEMORY IN IPINV',/                                               &
     & 1X,'CHANGE MEMK TO',I3)                                          
      ENDIF 
      IF (IWARN.EQ.1) THEN 
         WRITE (LUNFIL(2), 12346) 
         STOP 'DIMENSION EXCEEDED' 
12346 FORMAT (/1X,'A DIMENSION IS EXCEEDED IN ONE OF THE RESPONSE MATRIX&
     & SUBROUTINES')                                                    
      ENDIF 
!      ----------------------------                                     
      CALL MINV (H1, NH1, DET, KN, KM) 
!      ----------------------------                                     
      IF (NL.EQ.1) RETURN 
!      ----------------------------                                     
      CALL MINV (H2, NH2, DET, KN, KM) 
!      ----------------------------                                     
!                                                                       
      IONE = 1 
      IF ( (NL / 2 * 2) .EQ.NL) IONE = - 1 
      I = N * N + 1 
      DO 20 JG = 1, NG 
         JE = IONE 
         DO 20 JL = 1, NL 
            DO 19 IG = 1, NG 
               IE = IONE 
               DO 19 IL = 1, NL 
                  I = I - 1 
                  IF (IE * JE) 15, 15, 16 
   15             H1 (I) = 0.D0 
                  GOTO 19 
   16             IF (JE) 17, 17, 18 
   17             H1 (I) = H2 (IH2) 
                  IH2 = IH2 - 1 
                  GOTO 19 
   18             H1 (I) = H1 (IH1) 
                  IH1 = IH1 - 1 
   19       IE = - IE 
   20 JE = - JE 
!                                                                       
      RETURN 
      END SUBROUTINE IPINV                          
!EIGENV                      26 MAI 1982                             RES
!                                                                       
!                            **********                                 
!                            * EIGENV *                                 
!                            **********                                 
!                                                                       
      SUBROUTINE EIGENV (NG, A, ER, EI, VR, VI, IND) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      COMMON / LUNSCM / LUNFIL (4) 
      DIMENSION A ( * ), ER ( * ), EI ( * ), VR ( * ), VI ( * ),        &
      IND ( * )                                                         
!                                                                       
!         EIGENV FINDS THE EIGENVALUES (REAL - ER , IMAGINARY - EI) AND 
!         EIGENVECTORS (VR,VI) OF THE MATRIX A.  THE ORDER OF A IS NG.  
!                                                                       
      NNG = NG * NG 
      IF (NG.LT.3) GOTO 3 
      TU = 56.D0 
!      ---------------------------------------                          
      CALL EIGENP (NG, NG, A, TU, ER, EI, VR, VI, IND) 
!      ---------------------------------------                          
      DO 40 I = 1, NG 
         IF (IND (I) .LE.1) THEN 
            WRITE (LUNFIL(2), 44) IND (I), I 
   44 FORMAT(//1X,'EIGENP UNSUCCESSFUL, INDIC=',I1,' FOR EIGENVALUE',I3,&
     &/1X,'THE MEANING OF INDIC IS AS FOLLOWS:',                        &
     &/1X,'VALUE OF INDIC(I)    EIGENVALUE I    EIGENVECTOR I',         &
     &/1X,'       0               NOT FOUND       NOT FOUND',           &
     &/1X,'       1               FOUND           NOT FOUND',           &
     &/1X,'       2               FOUND           FOUND'/)              
         ENDIF 
   40 END DO 
      RETURN 
!                                                                       
    3 DO 2 I = 1, NNG 
         VR (I) = 1.D0 
    2 VI (I) = 0.D0 
      EI (1) = 0.D0 
      IF (NG - 1) 1, 1, 4 
    1 ER (1) = A (1) 
      GOTO 30 
    4 A1 = A (1) 
      A2 = A (2) 
      A3 = A (3) 
      A4 = A (4) 
      A23 = A2 * A3 
      EI (2) = 0.D0 
      IF (A23) 10, 6, 10 
    6 ER (1) = A1 
      ER (2) = A4 
      IF (A3) 8, 7, 8 
    7 VR (2) = A2 / (A1 - A4) 
      VR (3) = 0.D0 
      GOTO 30 
    8 VR (2) = 0.D0 
      VR (3) = A3 / (A4 - A1) 
      GOTO 30 
   10 Y = 0.5D0 * (A1 + A4) 
      Z = 0.5D0 * (A1 - A4) 
      Z2 = Z * Z 
      X = Z2 + A23 
      IF (X) 20, 12, 12 
!         REAL EIGENVALUES                                              
   12 IF (Z2 - 100.D0 * DABS (A23) ) 13, 13, 16 
   13 X = DSQRT (X) 
      ER (1) = Y + X 
      ER (2) = Y - X 
      IF (X) 15, 15, 14 
   14 VR (2) = ( - Z + X) / A3 
      VR (4) = ( - Z - X) / A3 
      GOTO 30 
   15 VR (1) = 0.D0 
      VR (3) = 0.D0 
      GOTO 30 
   16 X = A23 / Z2 
      S = 1.D0 - 0.25D0 * X * (1.D0 - 0.5D0 * X * (1.D0 - 0.625D0 * X) ) 
      R = 0.5D0 * Z * X * S 
      ER (1) = A1 + R 
      ER (2) = A4 - R 
      VR (2) = S * A2 / (A1 - A4) 
      VR (4) = (A4 - A1) / (A3 * S) 
      GOTO 30 
!         COMPLEX EIGENVALUES.                                          
   20 X = DSQRT ( - X) 
      ER (1) = Y 
      ER (2) = Y 
      EI (1) = X 
      EI (2) = - X 
      VR (2) = - Z / A3 
      VR (4) = VR (2) 
      VI (2) = X / A3 
      VI (4) = - VI (2) 
!                                                                       
   30 RETURN 
      END SUBROUTINE EIGENV                         
!FMOVT                       1  NOV 1982                             LAB
!                                                                       
!                            **********                                 
!                            * FMOVT  *                                 
!                            **********                                 
!                                                                       
!  PURPOSE    : TO MOVE DATA FROM ONE ARRAY TO ANOTHER.                 
!                                                                       
      SUBROUTINE FMOVT (A, B, NANT) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ), B ( * ) 
!                                                                       
      DO 10 I = 1, NANT 
   10 A (I) = B (I) 
!                                                                       
      RETURN 
      END SUBROUTINE FMOVT                          
!*                    *************                                     
!*                    *  G E L G  *                                     
!*                    *************                                     
!*                                                                      
!                                                                       
!        SUBROUTINE GELG                                                
!                                                                       
!                                                                       
!        PURPOSE                                                        
!           TO SOLVE A GENERAL SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS. 
!                                                                       
!        USAGE                                                          
!           CALL GELG(R,A,M,N,EPS,IER)                                  
!                                                                       
!        DESCRIPTION OF PARAMETERS                                      
!           R      - THE M BY N MATRIX OF RIGHT HAND SIDES.  (DESTROYED)
!                    ON RETURN R CONTAINS THE SOLUTION OF THE EQUATIONS.
!           A      - THE M BY M COEFFICIENT MATRIX.  (DESTROYED)        
!           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             
!           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             
!           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        
!                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        
!           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         
!                    IER=0  - NO ERROR,                                 
!                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR     
!                             PIVOT ELEMENT AT ANY ELIMINATION STEP     
!                             EQUAL TO 0,                               
!                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  
!                             CANCE INDICATED AT ELIMINATION STEP K+1,  
!                             WHERE PIVOT ELEMENT WAS LESS THAN OR      
!                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES 
!                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.  
!                                                                       
!        REMARKS                                                        
!           INPUT MATRICES R AND A ARE ASSUMED TO BE STORED COLUMNWISE  
!           IN M*N RESP. M*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN    
!           SOLUTION MATRIX R IS STORED COLUMNWISE TOO.                 
!           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS 
!           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS  
!           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -    
!           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL  
!           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE 
!           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS     
!           GIVEN IN CASE M=1.                                          
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH         
!           COMPLETE PIVOTING.                                          
!                                                                       
!     ..................................................................
!                                                                       
!                                                                       
!                                                                       
      SUBROUTINE GELG (R, A, M, N, EPS, IER) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ), R ( * ) 
!                                                                       
      IF (M) 23, 23, 1 
!                                                                       
!     SEARCH FOR GREATEST ELEMENT IN MATRIX A                           
    1 IER = 0 
      PIV = 0.D0 
      MM = M * M 
      NM = N * M 
      DO 3 L = 1, MM 
         TB = DABS (A (L) ) 
         IF (TB - PIV) 3, 3, 2 
    2    PIV = TB 
         I = L 
    3 END DO 
      TOL = EPS * PIV 
!     A(I) IS PIVOT ELEMENT. PIV CONTAINS THE ABSOLUTE VALUE OF A(I).   
!                                                                       
!                                                                       
!     START ELIMINATION LOOP                                            
      LST = 1 
      DO 17 K = 1, M 
!                                                                       
!     TEST ON SINGULARITY                                               
         IF (PIV) 23, 23, 4 
    4    IF (IER) 7, 5, 7 
    5    IF (PIV - TOL) 6, 6, 7 
    6    IER = K - 1 
    7    PIVI = 1.D0 / A (I) 
         J = (I - 1) / M 
         I = I - J * M - K 
         J = J + 1 - K 
!     I+K IS ROW-INDEX, J+K COLUMN-INDEX OF PIVOT ELEMENT               
!                                                                       
!     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      
         DO 8 L = K, NM, M 
            LL = L + I 
            TB = PIVI * R (LL) 
            R (LL) = R (L) 
    8    R (L) = TB 
!                                                                       
!     IS ELIMINATION TERMINATED                                         
         IF (K - M) 9, 18, 18 
!                                                                       
!     COLUMN INTERCHANGE IN MATRIX A                                    
    9    LEND = LST + M - K 
         IF (J) 12, 12, 10 
   10    II = J * M 
         DO 11 L = LST, LEND 
            TB = A (L) 
            LL = L + II 
            A (L) = A (LL) 
   11    A (LL) = TB 
!                                                                       
!     ROW INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A               
   12    DO 13 L = LST, MM, M 
            LL = L + I 
            TB = PIVI * A (LL) 
            A (LL) = A (L) 
   13    A (L) = TB 
!                                                                       
!     SAVE COLUMN INTERCHANGE INFORMATION                               
         A (LST) = J 
!                                                                       
!     ELEMENT REDUCTION AND NEXT PIVOT SEARCH                           
         PIV = 0.0D0 
         LST = LST + 1 
         J = 0 
         DO 16 II = LST, LEND 
            PIVI = - A (II) 
            IST = II + M 
            J = J + 1 
            DO 15 L = IST, MM, M 
               LL = L - J 
               A (L) = A (L) + PIVI * A (LL) 
               TB = DABS (A (L) ) 
               IF (TB - PIV) 15, 15, 14 
   14          PIV = TB 
               I = L 
   15       END DO 
            DO 16 L = K, NM, M 
               LL = L + J 
   16    R (LL) = R (LL) + PIVI * R (L) 
   17 LST = LST + M 
!     END OF ELIMINATION LOOP                                           
!                                                                       
!                                                                       
!     BACK SUBSTITUTION AND BACK INTERCHANGE                            
   18 IF (M - 1) 23, 22, 19 
   19 IST = MM + M 
      LST = M + 1 
      DO 21 I = 2, M 
         II = LST - I 
         IST = IST - LST 
         L = IST - M 
         L = A (L) + .5D0 
         DO 21 J = II, NM, M 
            TB = R (J) 
            LL = J 
            DO 20 K = IST, MM, M 
               LL = LL + 1 
   20       TB = TB - A (K) * R (LL) 
            K = J + L 
            R (J) = R (K) 
   21 R (K) = TB 
   22 RETURN 
!                                                                       
!                                                                       
!     ERROR RETURN                                                      
   23 IER = - 1 
      RETURN 
      END SUBROUTINE GELG                           
!MINV                        26 MAI 1982                             MAT
!                                                                       
!                            ********                                   
!                            * MINV *                                   
!                            ********                                   
!                                                                       
!  PURPOSE   : MINV INVERTS THE MATRIX A OF ORDER N BY THE STANDARD GAUS
!              METHOD. THE DETERMINANT IS ALSO CALCULATED = D. M AND L A
!              VECTORS OF LENGTH N.                                     
!-----------------------------------------------------------------------
      SUBROUTINE MINV (A, N, D, L, M) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ), L ( * ), M ( * ) 
!                                                                       
      D = 1.D0 
      NK = - N 
      DO 80 K = 1, N 
         NK = NK + N 
         L (K) = K 
         M (K) = K 
         KK = NK + K 
         BIGA = A (KK) 
         DO 20 J = K, N 
            IZ = N * (J - 1) 
            DO 20 I = K, N 
               IJ = IZ + I 
               IF (DABS (BIGA) - DABS (A (IJ) ) ) 15, 20, 20 
   15          BIGA = A (IJ) 
               L (K) = I 
               M (K) = J 
   20    CONTINUE 
         J = L (K) 
         IF (J - K) 35, 35, 25 
   25    KI = K - N 
         DO 30 I = 1, N 
            KI = KI + N 
            HOLD = - A (KI) 
            JI = KI - K + J 
            A (KI) = A (JI) 
   30    A (JI) = HOLD 
   35    I = M (K) 
         IF (I - K) 45, 45, 38 
   38    JP = N * (I - 1) 
         DO 40 J = 1, N 
            JK = NK + J 
            JI = JP + J 
            HOLD = - A (JK) 
            A (JK) = A (JI) 
   40    A (JI) = HOLD 
   45    IF (BIGA) 48, 46, 48 
   46    D = 0.D0 
         RETURN 
   48    DO 55 I = 1, N 
            IF (I - K) 50, 55, 50 
   50       IK = NK + I 
            A (IK) = A (IK) / ( - BIGA) 
   55    END DO 
         DO 65 I = 1, N 
            IK = NK + I 
            HOLD = A (IK) 
            IJ = I - N 
            DO 65 J = 1, N 
               IJ = IJ + N 
               IF (I - K) 60, 65, 60 
   60          IF (J - K) 62, 65, 62 
   62          KJ = IJ - I + K 
               A (IJ) = HOLD * A (KJ) + A (IJ) 
   65    CONTINUE 
         KJ = K - N 
         DO 75 J = 1, N 
            KJ = KJ + N 
            IF (J - K) 70, 75, 70 
   70       A (KJ) = A (KJ) / BIGA 
   75    END DO 
         D = D * BIGA 
         A (KK) = 1.D0 / BIGA 
   80 END DO 
      K = N 
  100 K = K - 1 
      IF (K) 150, 150, 105 
  105 I = L (K) 
      IF (I - K) 120, 120, 108 
  108 JQ = N * (K - 1) 
      JR = N * (I - 1) 
      DO 110 J = 1, N 
         JK = JQ + J 
         HOLD = A (JK) 
         JI = JR + J 
         A (JK) = - A (JI) 
  110 A (JI) = HOLD 
  120 J = M (K) 
      IF (J - K) 100, 100, 125 
  125 KI = K - N 
      DO 130 I = 1, N 
         KI = KI + N 
         HOLD = A (KI) 
         JI = KI - K + J 
         A (KI) = - A (JI) 
  130 A (JI) = HOLD 
      GOTO 100 
  150 RETURN 
      END SUBROUTINE MINV                           
!*                    *************                                     
!*                    * M U L T A *                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE MULTA (A, V, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ) 
!*                                                                      
!*    A=A*VALUE                                                         
!*                                                                      
      DO 1 I = 1, N 
    1 A (I) = A (I) * V 
      RETURN 
      END SUBROUTINE MULTA                          
!*                    *************                                     
!*                    *S E T V A L*                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE SETVAL (A, V, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ) 
!*                                                                      
!*    A=VALUE                                                           
!:                                                                      
      DO 1 I = 1, N 
    1 A (I) = V 
      RETURN 
      END SUBROUTINE SETVAL                         
!*                    *************                                     
!*                    * S E T I D *                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE SETID (A, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (N, * ) 
!*                                                                      
!*    A=E (IDENTITY MATRIX)                                             
!*                                                                      
      CALL SETVAL (A, 0.D0, N * N) 
      DO 1 I = 1, N 
    1 A (I, I) = 1.D0 
      RETURN 
      END SUBROUTINE SETID                          
!*                    *************                                     
!*                    * S U M A B *                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE SUMAB (A, B, C, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ), B ( * ), C ( * ) 
!*                                                                      
!*    C=A+B                                                             
!*                                                                      
      DO 1 I = 1, N 
    1 C (I) = A (I) + B (I) 
      RETURN 
      END SUBROUTINE SUMAB                          
!*                    *************                                     
!*                    * D I F A B *                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE DIFAB (A, B, C, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ), B ( * ), C ( * ) 
!*                                                                      
!*    C=A-B                                                             
!*$                                                                     
      DO 1 I = 1, N 
    1 C (I) = A (I) - B (I) 
      RETURN 
      END SUBROUTINE DIFAB                          
!*                    *************                                     
!*                    *M U L T A B*                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE MULTAB (A, B, C, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (N, * ), B (N, * ), C (N, * ) 
!*                                                                      
!*      C=A*B                                                           
!*                                                                      
      DO 1 I = 1, N 
         DO 1 J = 1, N 
            S = 0.D0 
            DO 2 K = 1, N 
    2       S = S + A (I, K) * B (K, J) 
    1 C (I, J) = S 
      RETURN 
      END SUBROUTINE MULTAB                         
!*                    *************                                     
!*                    *E Q U A T E*                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE EQUATE (A, B, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ), B ( * ) 
!*                                                                      
!*   A=B                                                                
!*                                                                      
      DO 1 I = 1, N 
    1 A (I) = B (I) 
      RETURN 
      END SUBROUTINE EQUATE                         
!*                    *************                                     
!*                    * O U T F L *                                     
!*                    *************                                     
!*                                                                      
      SUBROUTINE OUTFL (A, N, NAME) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (N, * ) 
      COMMON / LUNSCM / LUNFIL (4) 
      CHARACTER NAME * 6 
!*                                                                      
!*    OUTPUT OF A IN FLOATING POINT FORM                                
!*                                                                      
      WRITE (LUNFIL(2), 1) NAME 
    1 FORMAT(1X,'<',A6,'>') 
      DO 2 I = 1, N 
    2 WRITE (LUNFIL(2), 3) (A (I, J), J = 1, N) 
    3 FORMAT(1X,1P9D14.6) 
      RETURN 
      END SUBROUTINE OUTFL                          
!*                    ***************                                   
!*                    * O U T F L V *                                   
!*                    ***************                                   
!*                                                                      
      SUBROUTINE OUTFLV (A, N, NAME) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      COMMON / LUNSCM / LUNFIL (4) 
      DIMENSION A (N) 
      CHARACTER NAME * 6 
!*                                                                      
!*    OUTPUT OF A IN FLOATING POINT FORM                                
!*                                                                      
      WRITE (LUNFIL(2), 1) NAME 
    1 FORMAT(1X,'<',A6,'>') 
      WRITE (LUNFIL(2), 3) (A (J), J = 1, N) 
    3 FORMAT(1X,1P9D14.6) 
      RETURN 
      END SUBROUTINE OUTFLV                         
!-----------------------------------------------------------------------
!*                                                                      
!*                    ************                                      
!*                    *  EIGENP  *                                      
!*                    ************                                      
!*                                                                      
!EIGENP            EIGENP/EIGENP 28.8.70 G.B.                           
      SUBROUTINE EIGENP (N, NM, A, T, EVR, EVI, VECR, VECI, INDIC) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!     DOUBLE PRECISION D1,D2,D3,PRFACT                                  
!     INTEGER I,IVEC,J,K,K1,KON,L,L1,M,N,NM                             
!     REAL ENORM,EPS,EX,R,R1,T                                          
      DIMENSION A (NM, * ), VECR (NM, * ), VECI (NM, * ), EVR (NM),     &
      EVI (NM), INDIC (NM)                                              
      PARAMETER (MAXW = 100) 
      COMMON / LUNSCM / LUNFIL (4) 
      DIMENSION IWORK (MAXW), LOCAL (MAXW), PRFACT (MAXW), SUBDIA (MAXW)&
      , WORK1 (MAXW), WORK2 (MAXW), WORK (MAXW)                         
!                                                                       
! THIS SUBROUTINE FINDS ALL THE EIGENVALUES AND THE                     
! EIGENVECTORS OF A REAL GENERAL MATRIX OF ORDER N.                     
!                                                                       
! FIRST IN THE SUBROUTINE SCALE THE MATRIX IS SCALED SO THAT            
! THE CORRESPONDING ROWS AND COLUMNS ARE APPROXIMATELY                  
! BALANCED AND THEN THE MATRIX IS NORMALISED SO THAT THE                
! VALUE OF THE EUCLIDIAN NORM OF THE MATRIX IS EQUAL TO ONE.            
!                                                                       
! THE EIGENVALUES ARE COMPUTED BY THE QR DOUBLE-STEP METHOD             
! IN THE SUBROUTINE HESQR.                                              
! THE EIGENVECTORS ARE COMPUTED BY INVERSE ITERATION IN                 
! THE SUBROUTINE REALVE,FOR THE REAL EIGENVALUES,OR IN THE              
! SUBROUTINE COMPVE,FOR THE COMPLEX EIGENVALUES.                        
!                                                                       
! THE ELEMENTS OF THE MATRIX ARE TO BE STORED IN THE FIRST N            
! ROWS AND COLUMNS OF THE TWO DIMENSIONAL ARRAY A. THE                  
! ORIGINAL MATRIX IS DESTROYED BY THE SUBROUTINE.                       
! N IS THE ORDER OF THE MATRIX.                                         
! NM DEFINES THE FIRST DIMENSION OF THE TWO DIMENSIONAL                 
! ARRAYS A,VECR,VECI AND THE DIMENSION OF THE ONE                       
! DIMENSIONAL ARRAYS EVR,EVI AND INDIC. THEREFORE THE                   
! CALLING PROGRAM SHOULD CONTAIN THE FOLLOWING DECLARATION              
!     DIMENSION A(NM,NN),VECR(NM,NN),VECI(NM,NN),                       
!    1EVR(NM),EVI(NM),INDIC(NM)                                         
! WHERE NM AND NN ARE ANY NUMBERS EQUAL TO OR GREATER THAN N            
! NM AND NN ARE OF COURSE BOUNDED BY THE SIZE OF THE STORE.             
!                                                                       
! THE REAL PARAMETER T MUST BE SET EQUAL TO THE NUMBER OF               
! BINARY DIGITS IN THE MANTISSA OF A DOUBLE PRECISION                   
! FLOATING-POINT NUMBER.                                                
!                                                                       
! THE REAL PARTS OF THE N COMPUTED EIGENVALUES WILL BE FOUND            
! IN THE FIRST N PLACES OF THE ARRAY EVR AND THE IMAGINARY              
! PARTS IN THE FIRST N PLACES OF THE ARRAY EVI.                         
! THE REAL COMPONENTS OF THE NORMALISED EIGENVECTOR I                   
! (I_1,2,....,N) CORRESPONDING TO THE EIGENVALUE STORED IN              
! EVR(I) AND EVI(I) WILL BE FOUND IN THE FIRST N PLACES OF              
! THE COLUMN I OF THE TWO DIMENSIONAL ARRAY VECR AND THE                
! IMAGINARY COMPONENTS IN THE FIRST N PLACES OF THE COLUMN I            
! OF THE TWO DIMENSIONAL ARRAY VECI.                                    
!                                                                       
! THE REAL EIGENVECTOR IS NORMALISED SO THAT THE SUM OF THE             
! SQUARES OF THE COMPONENTS IS EQUAL TO ONE.                            
! THE COMPLEX EIGENVECTOR IS NORMALISED SO THAT THE                     
! COMPONENT WITH THE LARGEST VALUE IN MODULUS HAS ITS REAL              
! PART EQUAL TO ONE AND THE IMAGINARY PART EQUAL TO ZERO.               
! THE ARRAY INDIC INDICATES THE SUCCESS OF THE SUBROUTINE               
! EIGENP AS FOLLOWS                                                     
!     VALUE OF INDIC(I)    EIGENVALUE I    EIGENVECTOR I                
!            0               NOT FOUND       NOT FOUND                  
!            1               FOUND           NOT FOUND                  
!            2               FOUND           FOUND                      
!                                                                       
!                                                                       
      IF (NM.GT.MAXW) THEN 
         WRITE (LUNFIL(2), 1234) NM 
         STOP 'DIMENSION EXCEEDED' 
 1234 FORMAT(/1X,'DIMENSION EXCEEDED:INCREASE MAXW IN EIGENP TO MORE THA&
     &N',I6)                                                            
      ENDIF 
      IF (N.NE.1) GOTO 1 
      EVR (1) = A (1, 1) 
      EVI (1) = 0.0D0 
      VECR (1, 1) = 1.0D0 
      VECI (1, 1) = 0.0D0 
      INDIC (1) = 2 
      GOTO 25 
!                                                                       
    1 CALL SKALA (N, NM, A, VECI, PRFACT, ENORM) 
! THE COMPUTATION OF THE EIGENVALUES OF THE NORMALISED                  
! MATRIX.                                                               
      EX = DEXP ( - T * DLOG (2.0D0) ) 
      CALL HESQR (N, NM, A, VECI, EVR, EVI, SUBDIA, INDIC, EPS, EX) 
!                                                                       
! THE POSSIBLE DECOMPOSITION OF THE UPPER-HESSENBERG MATRIX             
! INTO THE SUBMATRICES OF LOWER ORDER IS INDICATED IN THE               
! ARRAY LOCAL. THE DECOMPOSITION OCCURS WHEN SOME                       
! SUBDIAGONAL ELEMENTS ARE IN MODULUS LESS THAN A SMALL                 
! POSITIVE NUMBER EPS DEFINED IN THE SUBROUTINE HESQR . THE             
! AMOUNT OF WORK IN THE EIGENVECTOR PROBLEM MAY BE                      
! DIMINISHED IN THIS WAY.                                               
      J = N 
      I = 1 
      LOCAL (1) = 1 
      IF (J.EQ.1) GOTO 4 
    2 IF (DABS (SUBDIA (J - 1) ) .GT.EPS) GOTO 3 
      I = I + 1 
      LOCAL (I) = 0 
    3 J = J - 1 
      LOCAL (I) = LOCAL (I) + 1 
      IF (J.NE.1) GOTO 2 
!                                                                       
! THE EIGENVECTOR PROBLEM.                                              
    4 K = 1 
      KON = 0 
      L = LOCAL (1) 
      M = N 
      DO 10 I = 1, N 
         IVEC = N - I + 1 
         IF (I.LE.L) GOTO 5 
         K = K + 1 
         M = N - L 
         L = L + LOCAL (K) 
    5    IF (INDIC (IVEC) .EQ.0) GOTO 10 
         IF (EVI (IVEC) .NE.0.0D0) GOTO 8 
!                                                                       
! TRANSFER OF AN UPPER-HESSENBERG MATRIX OF THE ORDER M FROM            
! THE ARRAYS VECI AND SUBDIA INTO THE ARRAY A.                          
         DO 7 K1 = 1, M 
            DO 6 L1 = K1, M 
    6       A (K1, L1) = VECI (K1, L1) 
            IF (K1.EQ.1) GOTO 7 
            A (K1, K1 - 1) = SUBDIA (K1 - 1) 
    7    END DO 
!                                                                       
! THE COMPUTATION OF THE REAL EIGENVECTOR IVEC OF THE UPPER-            
! HESSENBERG MATRIX CORRESPONDING TO THE REAL EIGENVALUE                
! EVR(IVEC).                                                            
         CALL REALVE (N, NM, M, IVEC, A, VECR, EVR, EVI, IWORK, WORK,   &
         INDIC, EPS, EX)                                                
         GOTO 10 
!                                                                       
! THE COMPUTATION OF THE COMPLEX EIGENVECTOR IVEC OF THE                
! UPPER-HESSENBERG MATRIX CORRESPONDING TO THE COMPLEX                  
! EIGENVALUE EVR(IVEC) / I*EVI(IVEC). IF THE VALUE OF KON IS            
! NOT EQUAL TO ZERO THEN THIS COMPLEX EIGENVECTOR HAS                   
! ALREADY BEEN FOUND FROM ITS CONJUGATE.                                
    8    IF (KON.NE.0) GOTO 9 
         KON = 1 
         CALL COMPVE (N, NM, M, IVEC, A, VECR, VECI, EVR, EVI, INDIC,   &
         IWORK, SUBDIA, WORK1, WORK2, WORK, EPS, EX)                    
         GOTO 10 
    9    KON = 0 
   10 END DO 
!                                                                       
! THE RECONSTRUCTION OF THE MATRIX USED IN THE REDUCTION OF             
! MATRIX A TO AN UPPER-HESSENBERG FORM BY HOUSEHOLDER METHOD            
      DO 12 I = 1, N 
         DO 11 J = I, N 
            A (I, J) = 0.0D0 
   11    A (J, I) = 0.0D0 
   12 A (I, I) = 1.0D0 
      IF (N.LE.2) GOTO 15 
      M = N - 2 
      DO 14 K = 1, M 
         L = K + 1 
         DO 14 J = 2, N 
            D1 = 0.0D0 
            DO 13 I = L, N 
               D2 = VECI (I, K) 
   13       D1 = D1 + D2 * A (J, I) 
            DO 14 I = L, N 
   14 A (J, I) = A (J, I) - VECI (I, K) * D1 
!                                                                       
! THE COMPUTATION OF THE EIGENVECTORS OF THE ORIGINAL NON-              
! SCALED MATRIX.                                                        
   15 KON = 1 
      DO 24 I = 1, N 
         L = 0 
         IF (EVI (I) .EQ.0.0D0) GOTO 16 
         L = 1 
         IF (KON.EQ.0) GOTO 16 
         KON = 0 
         GOTO 24 
   16    DO 18 J = 1, N 
            D1 = 0.0D0 
            D2 = 0.0D0 
            DO 17 K = 1, N 
               D3 = A (J, K) 
               D1 = D1 + D3 * VECR (K, I) 
               IF (L.EQ.0) GOTO 17 
               D2 = D2 + D3 * VECR (K, I - 1) 
   17       END DO 
            WORK (J) = D1 / PRFACT (J) 
            IF (L.EQ.0) GOTO 18 
            SUBDIA (J) = D2 / PRFACT (J) 
   18    END DO 
!                                                                       
! THE NORMALISATION OF THE EIGENVECTORS AND THE COMPUTATION             
! OF THE EIGENVALUES OF THE ORIGINAL NON-NORMALISED MATRIX.             
         IF (L.EQ.1) GOTO 21 
         D1 = 0.0D0 
         DO 19 M = 1, N 
   19    D1 = D1 + WORK (M) * WORK (M) 
         D1 = DSQRT (D1) 
         DO 20 M = 1, N 
            VECI (M, I) = 0.0D0 
   20    VECR (M, I) = WORK (M) / D1 
         EVR (I) = EVR (I) * ENORM 
         GOTO 24 
!                                                                       
   21    KON = 1 
         EVR (I) = EVR (I) * ENORM 
         EVR (I - 1) = EVR (I) 
         EVI (I) = EVI (I) * ENORM 
         EVI (I - 1) = - EVI (I) 
         R = 0.0D0 
         DO 22 J = 1, N 
            R1 = WORK (J) * WORK (J) + SUBDIA (J) * SUBDIA (J) 
            IF (R.GE.R1) GOTO 22 
            R = R1 
            L = J 
   22    END DO 
         D3 = WORK (L) 
         R1 = SUBDIA (L) 
         DO 23 J = 1, N 
            D1 = WORK (J) 
            D2 = SUBDIA (J) 
            VECR (J, I) = (D1 * D3 + D2 * R1) / R 
            VECI (J, I) = (D2 * D3 - D1 * R1) / R 
            VECR (J, I - 1) = VECR (J, I) 
   23    VECI (J, I - 1) = - VECI (J, I) 
   24 END DO 
!                                                                       
   25 RETURN 
      END SUBROUTINE EIGENP                         
!*                                                                      
!*                    ************                                      
!*                    *  SCALA   *                                      
!*                    ************                                      
!*                                                                      
!SCALE             SCALE/EIGENP 28.8.70.G.B.                            
      SUBROUTINE SKALA (N, NM, A, H, PRFACT, ENORM) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!     DOUBLE PRECISION COLUMN,FACTOR,FNORM,PRFACT,Q,ROW                 
!     INTEGER I,J,ITER,N,NCOUNT,NM                                      
!     REAL BOUND1,BOUND2,ENORM                                          
      DIMENSION A (NM, * ), H (NM, * ), PRFACT (NM) 
      COMMON / LUNSCM / LUNFIL (4) 
!                                                                       
! THIS SUBROUTINE STORES THE MATRIX OF THE ORDER N FROM THE             
! ARRAY A INTO THE ARRAY H. AFTERWARD THE MATRIX IN THE                 
! ARRAY A IS SCALED SO THAT THE QUOTIENT OF THE ABSOLUTE SUM            
! OF THE OFF-DIAGONAL ELEMENTS OF COLUMN I AND THE ABSOLUTE             
! SUM OF THE OFF-DIAGONAL ELEMENTS OF ROW I LIES WITHIN THE             
! VALUES OF BOUND1 AND BOUND2.                                          
! THE COMPONENT I OF THE EIGENVECTOR OBTAINED BY USING THE              
! SCALED MATRIX MUST BE DIVIDED BY THE VALUE FOUND IN THE               
! PRFACT(I) OF THE ARRAY PRFACT. IN THIS WAY THE EIGENVECTOR            
! OF THE NON-SCALED MATRIX IS OBTAINED.                                 
!                                                                       
! AFTER THE MATRIX IS SCALED IT IS NORMALISED SO THAT THE               
! VALUE OF THE EUCLIDIAN NORM IS EQUAL TO ONE.                          
! IF THE PROCESS OF SCALING WAS NOT SUCCESSFUL THE ORIGINAL             
! MATRIX FROM THE ARRAY H WOULD BE STORED BACK INTO A AND               
! THE EIGENPROBLEM WOULD BE SOLVED BY USING THIS MATRIX.                
! NM DEFINES THE FIRST DIMENSION OF THE ARRAYS A AND H. NM              
! MUST BE GREATER OR EQUAL TO N.                                        
! THE EIGENVALUES OF THE NORMALISED MATRIX MUST BE                      
! MULTIPLIED BY THE SCALAR ENORM IN ORDER THAT THEY BECOME              
! THE EIGENVALUES OF THE NON-NORMALISED MATRIX.                         
!                                                                       
      DO 2 I = 1, N 
         DO 1 J = 1, N 
    1    H (I, J) = A (I, J) 
    2 PRFACT (I) = 1.0D0 
      BOUND1 = 0.75D0 
      BOUND2 = 1.33D0 
      ITER = 0 
    3 NCOUNT = 0 
      DO 8 I = 1, N 
         COLUMN = 0.0D0 
         ROW = 0.0D0 
         DO 4 J = 1, N 
            IF (I.EQ.J) GOTO 4 
            COLUMN = COLUMN + DABS (A (J, I) ) 
            ROW = ROW + DABS (A (I, J) ) 
    4    END DO 
         IF (COLUMN.EQ.0.0D0) GOTO 5 
         IF (ROW.EQ.0.0D0) GOTO 5 
         Q = COLUMN / ROW 
         IF (Q.LT.BOUND1) GOTO 6 
         IF (Q.GT.BOUND2) GOTO 6 
    5    NCOUNT = NCOUNT + 1 
         GOTO 8 
    6    FACTOR = DSQRT (Q) 
         DO 7 J = 1, N 
            IF (I.EQ.J) GOTO 7 
            A (I, J) = A (I, J) * FACTOR 
            A (J, I) = A (J, I) / FACTOR 
    7    END DO 
         PRFACT (I) = PRFACT (I) * FACTOR 
    8 END DO 
      ITER = ITER + 1 
      IF (ITER.GT.30) GOTO 11 
      IF (NCOUNT.LT.N) GOTO 3 
!                                                                       
      FNORM = 0.0D0 
      DO 9 I = 1, N 
         DO 9 J = 1, N 
            Q = A (I, J) 
    9 FNORM = FNORM + Q * Q 
      FNORM = DSQRT (FNORM) 
      DO 10 I = 1, N 
         DO 10 J = 1, N 
   10 A (I, J) = A (I, J) / FNORM 
      ENORM = FNORM 
      GOTO 13 
!                                                                       
   11 DO 12 I = 1, N 
         PRFACT (I) = 1.0D0 
         DO 12 J = 1, N 
   12 A (I, J) = H (I, J) 
      ENORM = 1.0D0 
!                                                                       
   13 RETURN 
      END SUBROUTINE SKALA                          
!*                                                                      
!*                    ************                                      
!*                    *  HESQR   *                                      
!*                    ************                                      
!*                                                                      
!HESQR             HESQR/EIGENP28.8.70.G.B.                             
      SUBROUTINE HESQR (N, NM, A, H, EVR, EVI, SUBDIA, INDIC, EPS, EX) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!     DOUBLE PRECISION S,SR,SR2,X,Y,Z                                   
!     INTEGER I,J,K,L,M,MAXST,M1,N,NM,NS                                
!     REAL EPS,EX,R,SHIFT,T                                             
      DIMENSION A (NM, * ), H (NM, * ), EVR (NM), EVI (NM), SUBDIA (NM) 
      DIMENSION INDIC (NM) 
!                                                                       
! THIS SUBROUTINE FINDS ALL THE EIGENVALUES OF A REAL                   
! GENERAL MATRIX. THE ORIGINAL MATRIX A OF ORDER N IS                   
! REDUCED TO THE UPPER-HESSENGERG FORM H BY MEANS OF                    
! SIMILARITY TRANSFORMATIONS(HOUSEHOLDER METHOD). THE MATRIX            
! H IS PRESERVED IN THE UPPER HALF OF THE ARRAY H AND IN THE            
! ARRAY SUBDIA. THE SPECIAL VECTORS USED IN THE DEFINITION              
! THE LOWER PART OF THE ARRAY H.                                        
! OF THE HOUSEHOLDER TRANSFORMATION MATRICES ARE STORED IN              
! NM IS THE FIRST DIMENSION OF THE ARRAYS A AND H. NM MUST              
! BE EQUAL TO OR GREATER THAN N.                                        
! THE REAL PARTS OF THE N EIGENVALUES WILL BE FOUND IN THE              
! FIRST N PLACES OF THE ARRAY EVR,AND                                   
! THE IMAGINARY PARTS IN THE FIRST N PLACES OF THE ARRAY EVI            
! THE ARRAY INDIC INDICATES THE SUCCESS OF THE ROUTINE AS               
! FOLLOWS                                                               
!         VALUE OF INDIC(I)     EIGENVALUE I                            
!                0                NOT FOUND                             
!                1                  FOUND                               
! EPS IS A SMALL POSITIVE NUMBER THAT NUMERICALLY REPRESENTS            
! ZERO IN THE PROGRAM. EPS _ (EUCLIDIAN NORM OF H)*EX ,WHERE            
! EX _ 2**(-T). T IS THE NUMBER OF BINARY DIGITS IN THE                 
! MANTISSA OF A FLOATING POINT NUMBER.                                  
!                                                                       
!                                                                       
!                                                                       
! REDUCTION OF THE MATRIX A TO AN UPPER-HESSENBERG FORM H.              
! THERE ARE N-2 STEPS.                                                  
      IF (N - 2) 14, 1, 2 
    1 SUBDIA (1) = A (2, 1) 
      GOTO 14 
    2 M = N - 2 
      DO 12 K = 1, M 
         L = K + 1 
         S = 0.0D0 
         DO 3 I = L, N 
            H (I, K) = A (I, K) 
    3    S = S + DABS (A (I, K) ) 
         IF (S.NE.DABS (A (K + 1, K) ) ) GOTO 4 
         SUBDIA (K) = A (K + 1, K) 
         H (K + 1, K) = 0.0D0 
         GOTO 12 
    4    SR2 = 0.0D0 
         DO 5 I = L, N 
            SR = A (I, K) 
            SR = SR / S 
            A (I, K) = SR 
    5    SR2 = SR2 + SR * SR 
         SR = DSQRT (SR2) 
         IF (A (L, K) .LT.0.0D0) GOTO 6 
         SR = - SR 
    6    SR2 = SR2 - SR * A (L, K) 
         A (L, K) = A (L, K) - SR 
         H (L, K) = H (L, K) - SR * S 
         SUBDIA (K) = SR * S 
         X = S * DSQRT (SR2) 
         DO 7 I = L, N 
            H (I, K) = H (I, K) / X 
    7    SUBDIA (I) = A (I, K) / SR2 
! PREMULTIPLICATION BY THE MATRIX PR.                                   
         DO 9 J = L, N 
            SR = 0.0D0 
            DO 8 I = L, N 
    8       SR = SR + A (I, K) * A (I, J) 
            DO 9 I = L, N 
    9    A (I, J) = A (I, J) - SUBDIA (I) * SR 
! POSTMULTIPLICATION BY THE MATRIX PR.                                  
         DO 11 J = 1, N 
            SR = 0.0D0 
            DO 10 I = L, N 
   10       SR = SR + A (J, I) * A (I, K) 
            DO 11 I = L, N 
   11    A (J, I) = A (J, I) - SUBDIA (I) * SR 
   12 END DO 
      DO 13 K = 1, M 
   13 A (K + 1, K) = SUBDIA (K) 
! TRANSFER OF THE UPPER HALF OF THE MATRIX A INTO THE                   
! ARRAY H AND THE CALCULATION OF THE SMALL POSITIVE NUMBER              
! EPS.                                                                  
      SUBDIA (N - 1) = A (N, N - 1) 
   14 EPS = 0.0D0 
      DO 15 K = 1, N 
         INDIC (K) = 0 
         IF (K.NE.N) EPS = EPS + SUBDIA (K) * SUBDIA (K) 
         DO 15 I = K, N 
            H (K, I) = A (K, I) 
   15 EPS = EPS + A (K, I) * A (K, I) 
      EPS = EX * DSQRT (EPS) 
!                                                                       
! THE QR ITERATIVE PROCESS. THE UPPER-HESSENBERG MATRIX H IS            
! REDUCED TO THE UPPER/MODIFIED TRIANGULAR FORM.                        
!                                                                       
! DETERMINATION OF THE SHIFT OF ORIGIN FOR THE FIRST STEP OF            
! THE QR ITERATIVE PROCESS.                                             
      SHIFT = A (N, N - 1) 
      IF (N.LE.2) SHIFT = 0.0D0 
      IF (A (N, N) .NE.0.0D0) SHIFT = 0.0D0 
      IF (A (N - 1, N) .NE.0.0D0) SHIFT = 0.0D0 
      IF (A (N - 1, N - 1) .NE.0.0D0) SHIFT = 0.0D0 
      M = N 
      NS = 0 
      MAXST = N * 10 
!                                                                       
! TESTING IF THE UPPER HALF OF THE MATRIX IS EQUAL TO ZERO.             
! IF IT IS EQUAL TO ZERO THE QR PROCESS IS NOT NECESSARY.               
      DO 16 I = 2, N 
         DO 16 K = I, N 
            IF (A (I - 1, K) .NE.0.0D0) GOTO 18 
   16 CONTINUE 
      DO 17 I = 1, N 
         INDIC (I) = 1 
         EVR (I) = A (I, I) 
   17 EVI (I) = 0.0D0 
      GOTO 37 
!                                                                       
! START THE MAIN LOOP OF THE QR PROCESS.                                
   18 K = M - 1 
      M1 = K 
      I = K 
! FIND ANY DECOMPOSITIONS OF THE MATRIX.                                
! JUMP TO 34 IF THE LAST SUBMATRIX OF THE DECOMPOSITION IS              
! OF THE ORDER ONE.                                                     
! JUMP TO 35 IF THE LAST SUBMATRIX OF THE DECOMPOSITION IS              
! OF THE ORDER TWO.                                                     
      IF (K) 37, 34, 19 
   19 IF (DABS (A (M, K) ) .LE.EPS) GOTO 34 
      IF (M - 2.EQ.0) GOTO 35 
   20 I = I - 1 
      IF (DABS (A (K, I) ) .LE.EPS) GOTO 21 
      K = I 
      IF (K.GT.1) GOTO 20 
   21 IF (K.EQ.M1) GOTO 35 
! TRANSFORMATION OF THE MATRIX OF THE ORDER GREATER THAN TWO            
      S = A (M, M) + A (M1, M1) + SHIFT 
      SR = A (M, M) * A (M1, M1) - A (M, M1) * A (M1, M) + 0.25D0 *     &
      SHIFT * SHIFT                                                     
      A (K + 2, K) = 0.0D0 
! CALCULATE X1,Y1,Z1,FOR THE SUBMATRIX OBTAINED BY THE                  
! DECOMPOSITION.                                                        
      X = A (K, K) * (A (K, K) - S) + A (K, K + 1) * A (K + 1, K)       &
      + SR                                                              
      Y = A (K + 1, K) * (A (K, K) + A (K + 1, K + 1) - S) 
      R = DABS (X) + DABS (Y) 
      IF (R.EQ.0.0D0) SHIFT = A (M, M - 1) 
      IF (R.EQ.0.0D0) GOTO 21 
      Z = A (K + 2, K + 1) * A (K + 1, K) 
      SHIFT = 0.0D0 
      NS = NS + 1 
!                                                                       
! THE LOOP FOR ONE STEP OF THE QR PROCESS.                              
      DO 33 I = K, M1 
         IF (I.EQ.K) GOTO 22 
! CALCULATE XR,YR,ZR.                                                   
         X = A (I, I - 1) 
         Y = A (I + 1, I - 1) 
         Z = 0.0D0 
         IF (I + 2.GT.M) GOTO 22 
         Z = A (I + 2, I - 1) 
   22    SR2 = DABS (X) + DABS (Y) + DABS (Z) 
         IF (SR2.EQ.0.0D0) GOTO 23 
         X = X / SR2 
         Y = Y / SR2 
         Z = Z / SR2 
   23    S = DSQRT (X * X + Y * Y + Z * Z) 
         IF (X.LT.0.0D0) GOTO 24 
         S = - S 
   24    IF (I.EQ.K) GOTO 25 
         A (I, I - 1) = S * SR2 
   25    IF (SR2.NE.0.0D0) GOTO 26 
         IF (I + 3.GT.M) GOTO 33 
         GOTO 32 
   26    SR = 1.0D0 - X / S 
         S = X - S 
         X = Y / S 
         Y = Z / S 
! PREMULTIPLICATION BY THE MATRIX PR.                                   
         DO 28 J = I, M 
            S = A (I, J) + A (I + 1, J) * X 
            IF (I + 2.GT.M) GOTO 27 
            S = S + A (I + 2, J) * Y 
   27       S = S * SR 
            A (I, J) = A (I, J) - S 
            A (I + 1, J) = A (I + 1, J) - S * X 
            IF (I + 2.GT.M) GOTO 28 
            A (I + 2, J) = A (I + 2, J) - S * Y 
   28    END DO 
! POSTMULTIPLICATION BY THE MATRIX PR.                                  
         L = I + 2 
         IF (I.LT.M1) GOTO 29 
         L = M 
   29    DO 31 J = K, L 
            S = A (J, I) + A (J, I + 1) * X 
            IF (I + 2.GT.M) GOTO 30 
            S = S + A (J, I + 2) * Y 
   30       S = S * SR 
            A (J, I) = A (J, I) - S 
            A (J, I + 1) = A (J, I + 1) - S * X 
            IF (I + 2.GT.M) GOTO 31 
            A (J, I + 2) = A (J, I + 2) - S * Y 
   31    END DO 
         IF (I + 3.GT.M) GOTO 33 
         S = - A (I + 3, I + 2) * Y * SR 
   32    A (I + 3, I) = S 
         A (I + 3, I + 1) = S * X 
         A (I + 3, I + 2) = S * Y + A (I + 3, I + 2) 
   33 END DO 
!                                                                       
      IF (NS.GT.MAXST) GOTO 37 
      GOTO 18 
!                                                                       
! COMPUTE THE LAST EIGENVALUE.                                          
   34 EVR (M) = A (M, M) 
      EVI (M) = 0.0D0 
      INDIC (M) = 1 
      M = K 
      GOTO 18 
!                                                                       
! COMPUTE THE EIGENVALUE OF THE LAST 2*2 MATRIX OBTAINED BY             
! THE DECOMPOSITION.                                                    
   35 R = 0.5D0 * (A (K, K) + A (M, M) ) 
      S = 0.5D0 * (A (M, M) - A (K, K) ) 
      S = S * S + A (K, M) * A (M, K) 
      INDIC (K) = 1 
      INDIC (M) = 1 
      IF (S.LT.0.0D0) GOTO 36 
      T = DSQRT (S) 
      EVR (K) = R - T 
      EVR (M) = R + T 
      EVI (K) = 0.0D0 
      EVI (M) = 0.0D0 
      M = M - 2 
      GOTO 18 
   36 T = DSQRT ( - S) 
      EVR (K) = R 
      EVI (K) = T 
      EVR (M) = R 
      EVI (M) = - T 
      M = M - 2 
      GOTO 18 
!                                                                       
   37 RETURN 
      END SUBROUTINE HESQR                          
!*                                                                      
!*                    ************                                      
!*                    *  REALVE  *                                      
!*                    ************                                      
!*                                                                      
!REALVE            REALVE/EIGENP 28.8.70.G.B.                           
      SUBROUTINE REALVE (N, NM, M, IVEC, A, VECR, EVR, EVI, IWORK, WORK,&
      INDIC, EPS, EX)                                                   
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!     DOUBLE PRECISION S,SR                                             
!     INTEGER I,IVEC,ITER,J,K,L,M,N,NM,NS                               
!     REAL BOUND,EPS,EVALUE,EX,PREVIS,R,R1,T                            
      DIMENSION A (NM, * ), VECR (NM, * ), EVR (NM) 
      DIMENSION EVI (NM), IWORK (NM), WORK (NM), INDIC (NM) 
!                                                                       
! THIS SUBROUTINE FINDS THE REAL EIGENVECTOR OF THE REAL                
! UPPER-HESSENBERG MATRIX IN THE ARRAY A,CORRESPONDING TO               
! THE REAL EIGENVALUE STORED IN EVR(IVEC). THE INVERSE                  
! ITERATION METHOD IS USED.                                             
! NOTE THE MATRIX IN A IS DESTROYED BY THE SUBROUTINE.                  
! N IS THE ORDER OF THE UPPER-HESSENBERG MATRIX.                        
! NM DEFINES THE FIRST DIMENSION OF THE TWO DIMENSIONAL                 
! ARRAYS A AND VECR. NM MUST BE EQUAL TO OR GREATER THAN N.             
! M IS THE ORDER OF THE SUBMATRIX OBTAINED BY A SUITABLE                
! DECOMPOSITION OF THE UPPER-HESSENBERG MATRIX IF SOME                  
! SUBDIAGONAL ELEMENTS ARE EQUAL TO ZERO. THE VALUE OF M IS             
! CHOSEN SO THAT THE LAST N-M COMPONENTS OF THE EIGENVECTOR             
! ARE ZERO.                                                             
! IVEC GIVES THE POSITION OF THE EIGENVALUE IN THE ARRAY EVR            
! FOR WHICH THE CORRESPONDING EIGENVECTOR IS COMPUTED.                  
! THE ARRAY EVI WOULD CONTAIN THE IMAGINARY PARTS OF THE N              
! EIGENVALUES IF THEY EXISTED.                                          
!                                                                       
! THE M COMPONENTS OF THE COMPUTED REAL EIGENVECTOR WILL BE             
! FOUND IN THE FIRST M PLACES OF THE COLUMN IVEC OF THE TWO             
! DIMENSIONAL ARRAY VECR.                                               
!                                                                       
! IWORK AND WORK ARE THE WORKING STORES USED DURING THE                 
! GAUSSIAN ELIMINATION AND BACKSUBSTITUTION PROCESS.                    
! THE ARRAY INDIC INDICATES THE SUCCESS OF THE ROUTINE AS               
! FOLLOWS                                                               
!     VALUE OF INDIC(I)     EIGENVECTOR I                               
!            1                NOT FOUND                                 
!            2                  FOUND                                   
! EPS IS A SMALL POSITIVE NUMBER THAT NUMERICALLY REPRESENTS            
! ZERO IN THE PROGRAM. EPS _ (EUCLIDIAN NORM OF A)*EX,WHERE             
! EX _ 2**(-T). T IS THE NUMBER OF BINARY DIGITS IN THE                 
! MANTISSA OF A FLOATING POINT NUMBER.                                  
      VECR (1, IVEC) = 1.0D0 
      IF (M.EQ.1) GOTO 24 
! SMALL PERTURBATION OF EQUAL EIGENVALUES TO OBTAIN A FULL              
! SET OF EIGENVECTORS.                                                  
      EVALUE = EVR (IVEC) 
      IF (IVEC.EQ.M) GOTO 2 
      K = IVEC + 1 
      R = 0.0D0 
      DO 1 I = K, M 
         IF (EVALUE.NE.EVR (I) ) GOTO 1 
         IF (EVI (I) .NE.0.0D0) GOTO 1 
         R = R + 3.0D0 
    1 END DO 
      EVALUE = EVALUE+R * EX 
    2 DO 3 K = 1, M 
    3 A (K, K) = A (K, K) - EVALUE 
!                                                                       
! GAUSSIAN ELIMINATION OF THE UPPER-HESSENBERG MATRIX A. ALL            
! ROW INTERCHANGES ARE INDICATED IN THE ARRAY IWORK.ALL THE             
! MULTIPLIERS ARE STORED AS THE SUBDIAGONAL ELEMENTS OF A.              
      K = M - 1 
      DO 8 I = 1, K 
         L = I + 1 
         IWORK (I) = 0 
         IF (A (I + 1, I) .NE.0.0D0) GOTO 4 
         IF (A (I, I) .NE.0.0D0) GOTO 8 
         A (I, I) = EPS 
         GOTO 8 
    4    IF (DABS (A (I, I) ) .GE.DABS (A (I + 1, I) ) ) GOTO 6 
         IWORK (I) = 1 
         DO 5 J = I, M 
            R = A (I, J) 
            A (I, J) = A (I + 1, J) 
    5    A (I + 1, J) = R 
    6    R = - A (I + 1, I) / A (I, I) 
         A (I + 1, I) = R 
         DO 7 J = L, M 
    7    A (I + 1, J) = A (I + 1, J) + R * A (I, J) 
    8 END DO 
      IF (A (M, M) .NE.0.0D0) GOTO 9 
      A (M, M) = EPS 
!                                                                       
! THE VECTOR (1,1,...,1) IS STORED IN THE PLACE OF THE RIGHT            
! HAND SIDE COLUMN VECTOR.                                              
    9 DO 11 I = 1, N 
         IF (I.GT.M) GOTO 10 
         WORK (I) = 1.0D0 
         GOTO 11 
   10    WORK (I) = 0.0D0 
   11 END DO 
!                                                                       
! THE INVERSE ITERATION IS PERFORMED ON THE MATRIX UNTIL THE            
! INFINITE NORM OF THE RIGHT-HAND SIDE VECTOR IS GREATER                
! THAN THE BOUND DEFINED AS 0.01/(N*EX).                                
      BOUND = 0.01D0 / (EX * DBLE (N) ) 
      NS = 0 
      ITER = 1 
!                                                                       
! THE BACKSUBSTITUTION.                                                 
   12 R = 0.0D0 
      DO 15 I = 1, M 
         J = M - I + 1 
         S = WORK (J) 
         IF (J.EQ.M) GOTO 14 
         L = J + 1 
         DO 13 K = L, M 
            SR = WORK (K) 
   13    S = S - SR * A (J, K) 
   14    WORK (J) = S / A (J, J) 
         T = DABS (WORK (J) ) 
         IF (R.GE.T) GOTO 15 
         R = T 
   15 END DO 
!                                                                       
! THE COMPUTATION OF THE RIGHT-HAND SIDE VECTOR FOR THE NEW             
! ITERATION STEP.                                                       
      DO 16 I = 1, M 
   16 WORK (I) = WORK (I) / R 
!                                                                       
! THE COMPUTATION OF THE RESIDUALS AND COMPARISON OF THE                
! RESIDUALS OF THE TWO SUCCESSIVE STEPS OF THE INVERSE                  
! ITERATION.IF THE INFINITE NORM OF THE RESIDUAL VECTOR IS              
! GREATER THAN THE INFINITE NORM OF THE PREVIOUS RESIDUAL               
! VECTOR THE COMPUTED EIGENVECTOR OF THE PREVIOUS STEP IS               
! TAKEN AS THE FINAL EIGENVECTOR.                                       
      R1 = 0.0D0 
      DO 18 I = 1, M 
         T = 0.0D0 
         DO 17 J = I, M 
   17    T = T + A (I, J) * WORK (J) 
         T = DABS (T) 
         IF (R1.GE.T) GOTO 18 
         R1 = T 
   18 END DO 
      IF (ITER.EQ.1) GOTO 19 
      IF (PREVIS.LE.R1) GOTO 24 
   19 DO 20 I = 1, M 
   20 VECR (I, IVEC) = WORK (I) 
      PREVIS = R1 
      IF (NS.EQ.1) GOTO 24 
      IF (ITER.GT.6) GOTO 25 
      ITER = ITER + 1 
      IF (R.LT.BOUND) GOTO 21 
      NS = 1 
!                                                                       
! GAUSSIAN ELIMINATION OF THE RIGHT-HAND SIDE VECTOR.                   
   21 K = M - 1 
      DO 23 I = 1, K 
         R = WORK (I + 1) 
         IF (IWORK (I) .EQ.0) GOTO 22 
         WORK (I + 1) = WORK (I) + WORK (I + 1) * A (I + 1, I) 
         WORK (I) = R 
         GOTO 23 
   22    WORK (I + 1) = WORK (I + 1) + WORK (I) * A (I + 1, I) 
   23 END DO 
      GOTO 12 
!                                                                       
   24 INDIC (IVEC) = 2 
   25 IF (M.EQ.N) GOTO 27 
      J = M + 1 
      DO 26 I = J, N 
   26 VECR (I, IVEC) = 0.0D0 
   27 RETURN 
      END SUBROUTINE REALVE                         
!*                                                                      
!*                    ************                                      
!*                    *  COMPVE  *                                      
!*                    ************                                      
!*                                                                      
!COMPVE            COMPVE/EIGENP 28.8.70.G.B.                           
      SUBROUTINE COMPVE (N, NM, M, IVEC, A, VECR, H, EVR, EVI, INDIC,   &
      IWORK, SUBDIA, WORK1, WORK2, WORK, EPS, EX)                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!     DOUBLE PRECISION D,D1                                             
!     INTEGER I,I1,I2,ITER,IVEC,J,K,L,M,N,NM,NS                         
!     REAL B,BOUND,EPS,ETA,EX,FKSI,PREVIS,R,S,U,V                       
      DIMENSION A (NM, * ), VECR (NM, * ), H (NM, * ), EVR (NM),        &
      EVI (NM), INDIC (NM), IWORK (NM), SUBDIA (NM), WORK1 (NM),        &
      WORK2 (NM), WORK (NM)                                             
!                                                                       
! THIS SUBROUTINE FINDS THE COMPLEX EIGENVECTOR OF THE REAL             
! UPPER-HESSENBERG MATRIX OF ORDER N CORRESPONDING TO THE               
! COMPLEX EIGENVALUE WITH THE REAL PART IN EVR(IVEC) AND THE            
! CORRESPONDING IMAGINARY PART IN EVI(IVEC). THE INVERSE                
! ITERATION METHOD IS USED MODIFIED TO AVOID THE USE OF                 
! COMPLEX ARITHMETIC.                                                   
! THE MATRIX ON WHICH THE INVERSE ITERATION IS PERFORMED IS             
! BUILT UP IN THE ARRAY A BY USING THE UPPER-HESSENBERG                 
! MATRIX PRESERVED IN THE UPPER HALF OF THE ARRAY H AND IN              
! THE ARRAY SUBDIA.                                                     
! NM DEFINES THE FIRST DIMENSION OF THE TWO DIMENSIONAL                 
! ARRAYS A,VECR AND H. NM MUST BE EQUAL TO OR GREATER                   
! THAN N.                                                               
! M IS THE ORDER OF THE SUBMATRIX OBTAINED BY A SUITABLE                
! DECOMPOSITION OF THE UPPER-HESSENBERG MATRIX IF SOME                  
! SUBDIAGONAL ELEMENTS ARE EQUAL TO ZERO. THE VALUE OF M IS             
! CHOSEN SO THAT THE LAST N-M COMPONENTS OF THE COMPLEX                 
! EIGENVECTOR ARE ZERO.                                                 
!                                                                       
! THE REAL PARTS OF THE FIRST M COMPONENTS OF THE COMPUTED              
! COMPLEX EIGENVECTOR WILL BE FOUND IN THE FIRST M PLACES OF            
! THE COLUMN WHOSE TOP ELEMENT IS VECR(1,IVEC) AND THE                  
! CORRESPONDING IMAGINARY PARTS OF THE FIRST M COMPONENTS OF            
! THE COMPLEX EIGENVECTOR WILL BE FOUND IN THE FIRST M                  
! PLACES OF THE COLUMN WHOSE TOP ELEMENT IS VECR(1,IVEC-1).             
!                                                                       
! THE ARRAY INDIC INDICATES THE SUCCESS OF THE ROUTINE AS               
! FOLLOWS                                                               
!     VALUE OF INDIC(I)     EIGENVECTOR I                               
!            1                NOT FOUND                                 
!            2                  FOUND                                   
! THE ARRAYS IWORK,WORK1,WORK2 AND WORK ARE THE WORKING                 
! STORES USED DURING THE INVERSE ITERATION PROCESS.                     
! EPS IS A SMALL POSITIVE NUMBER THAT NUMERICALLY REPRESENTS            
! ZERO IN THE PROGRAM. EPS _ (EUCLIDIAN NORM OF H)*EX, WHERE            
! EX _ 2**(-T). T IS THE NUMBER OF BINARY DIGITS IN THE                 
! MANTISSA OF A FLOATING POINT NUMBER.                                  
!                                                                       
      FKSI = EVR (IVEC) 
      ETA = EVI (IVEC) 
! THE MODIFICATION OF THE EIGENVALUE (FKSI / I*ETA) IF MORE             
! EIGENVALUES ARE EQUAL.                                                
      IF (IVEC.EQ.M) GOTO 2 
      K = IVEC + 1 
      R = 0.0D0 
      DO 1 I = K, M 
         IF (FKSI.NE.EVR (I) ) GOTO 1 
         IF (DABS (ETA) .NE.DABS (EVI (I) ) ) GOTO 1 
         R = R + 3.0D0 
    1 END DO 
      R = R * EX 
      FKSI = FKSI + R 
      ETA = ETA + R 
!                                                                       
! THE MATRIX  ((H-FKSI*I)*(H-FKSI*I) / (ETA*ETA)*I)  IS                 
! STORED INTO THE ARRAY A.                                              
    2 R = FKSI * FKSI + ETA * ETA 
      S = 2.0D0 * FKSI 
      L = M - 1 
      DO 5 I = 1, M 
         DO 4 J = I, M 
            D = 0.0D0 
            A (J, I) = 0.0D0 
            DO 3 K = I, J 
    3       D = D+H (I, K) * H (K, J) 
    4    A (I, J) = D-S * H (I, J) 
    5 A (I, I) = A (I, I) + R 
      DO 9 I = 1, L 
         R = SUBDIA (I) 
         A (I + 1, I) = - S * R 
         I1 = I + 1 
         DO 6 J = 1, I1 
    6    A (J, I) = A (J, I) + R * H (J, I + 1) 
         IF (I.EQ.1) GOTO 7 
         A (I + 1, I - 1) = R * SUBDIA (I - 1) 
    7    DO 8 J = I, M 
    8    A (I + 1, J) = A (I + 1, J) + R * H (I, J) 
    9 END DO 
!                                                                       
! THE GAUSSIAN ELIMINATION OF THE MATRIX                                
! ((H-FKSI*I)*(H-FKSI*I) / (ETA*ETA)*I) IN THE ARRAY A. THE             
! ROW INTERCHANGES THAT OCCUR ARE INDICATED IN THE ARRAY                
! IWORK. ALL THE MULTIPLIERS ARE STORED IN THE FIRST AND IN             
! THE SECOND SUBDIAGONAL OF THE ARRAY A.                                
      K = M - 1 
      DO 18 I = 1, K 
         I1 = I + 1 
         I2 = I + 2 
         IWORK (I) = 0 
         IF (I.EQ.K) GOTO 10 
         IF (A (I + 2, I) .NE.0.0D0) GOTO 11 
   10    IF (A (I + 1, I) .NE.0.0D0) GOTO 11 
         IF (A (I, I) .NE.0.0D0) GOTO 18 
         A (I, I) = EPS 
         GOTO 18 
!                                                                       
   11    IF (I.EQ.K) GOTO 12 
         IF (DABS (A (I + 1, I) ) .GE.DABS (A (I + 2, I) ) ) GOTO 12 
         IF (DABS (A (I, I) ) .GE.DABS (A (I + 2, I) ) ) GOTO 16 
         L = I + 2 
         IWORK (I) = 2 
         GOTO 13 
   12    IF (DABS (A (I, I) ) .GE.DABS (A (I + 1, I) ) ) GOTO 15 
         L = I + 1 
         IWORK (I) = 1 
!                                                                       
   13    DO 14 J = I, M 
            R = A (I, J) 
            A (I, J) = A (L, J) 
   14    A (L, J) = R 
   15    IF (I.NE.K) GOTO 16 
         I2 = I1 
   16    DO 17 L = I1, I2 
            R = - A (L, I) / A (I, I) 
            A (L, I) = R 
            DO 17 J = I1, M 
   17    A (L, J) = A (L, J) + R * A (I, J) 
   18 END DO 
      IF (A (M, M) .NE.0.0D0) GOTO 19 
      A (M, M) = EPS 
!                                                                       
! THE VECTOR (1,1,...,1) IS STORED INTO THE RIGHT-HAND SIDE             
! VECTORS VECR( ,IVEC) AND VECR( ,IVEC-1) REPRESENTING THE              
! COMPLEX RIGHT-HAND SIDE VECTOR.                                       
   19 DO 21 I = 1, N 
         IF (I.GT.M) GOTO 20 
         VECR (I, IVEC) = 1.0D0 
         VECR (I, IVEC - 1) = 1.0D0 
         GOTO 21 
   20    VECR (I, IVEC) = 0.0D0 
         VECR (I, IVEC - 1) = 0.0D0 
   21 END DO 
!                                                                       
! THE INVERSE ITERATION IS PERFORMED ON THE MATRIX UNTIL THE            
! INFINITE NORM OF THE RIGHT-HAND SIDE VECTOR IS GREATER                
! THAN THE BOUND DEFINED AS 0.01/(N*EX).                                
      BOUND = 0.01D0 / (EX * DBLE (N) ) 
      NS = 0 
      ITER = 1 
      DO 22 I = 1, M 
   22 WORK (I) = H (I, I) - FKSI 
!                                                                       
! THE SEQUENCE OF THE COMPLEX VECTORS Z(S) _ P(S)/I*Q(S) AND            
! W(S/1)_ U(S/1)/I*V(S/1) IS GIVEN BY THE RELATIONS                     
! (A - (FKSI-I*ETA)*I)*W(S/1) _ Z(S) AND                                
!     C Z(S/1) _ W(S/1)/MAX(W(S/1)).                                    
! THE FINAL W(S) IS TAKEN AS THE COMPUTED EIGENVECTOR.                  
!                                                                       
! THE COMPUTATION OF THE RIGHT-HAND SIDE VECTOR                         
! (A-FKSI*I)*P(S)-ETA*Q(S). A IS AN UPPER-HESSENBERG MATRIX.            
   23 DO 27 I = 1, M 
         D = WORK (I) * VECR (I, IVEC) 
         IF (I.EQ.1) GOTO 24 
         D = D+SUBDIA (I - 1) * VECR (I - 1, IVEC) 
   24    L = I + 1 
         IF (L.GT.M) GOTO 26 
         DO 25 K = L, M 
   25    D = D+H (I, K) * VECR (K, IVEC) 
   26    VECR (I, IVEC - 1) = D-ETA * VECR (I, IVEC - 1) 
   27 END DO 
!                                                                       
! GAUSSIAN ELIMINATION OF THE RIGHT-HAND SIDE VECTOR.                   
      K = M - 1 
      DO 28 I = 1, K 
         L = I + IWORK (I) 
         R = VECR (L, IVEC - 1) 
         VECR (L, IVEC - 1) = VECR (I, IVEC - 1) 
         VECR (I, IVEC - 1) = R 
         VECR (I + 1, IVEC - 1) = VECR (I + 1, IVEC - 1) + A (I + 1, I) &
         * R                                                            
         IF (I.EQ.K) GOTO 28 
         VECR (I + 2, IVEC - 1) = VECR (I + 2, IVEC - 1) + A (I + 2, I) &
         * R                                                            
   28 END DO 
!                                                                       
! THE COMPUTATION OF THE REAL PART U(S/1) OF THE COMPLEX                
! VECTOR W(S/1). THE VECTOR U(S/1) IS OBTAINED AFTER THE                
! BACKSUBSTITUTION.                                                     
      DO 31 I = 1, M 
         J = M - I + 1 
         D = VECR (J, IVEC - 1) 
         IF (J.EQ.M) GOTO 30 
         L = J + 1 
         DO 29 K = L, M 
            D1 = A (J, K) 
   29    D = D-D1 * VECR (K, IVEC - 1) 
   30    VECR (J, IVEC - 1) = D / A (J, J) 
   31 END DO 
!                                                                       
! THE COMPUTATION OF THE IMAGINARY PART V(S/1) OF THE VECTOR            
! W(S/1), WHERE V(S/1) _ (P(S)-(A-FKSI*I)*U(S/1))/ETA.                  
      DO 35 I = 1, M 
         D = WORK (I) * VECR (I, IVEC - 1) 
         IF (I.EQ.1) GOTO 32 
         D = D+SUBDIA (I - 1) * VECR (I - 1, IVEC - 1) 
   32    L = I + 1 
         IF (L.GT.M) GOTO 34 
         DO 33 K = L, M 
   33    D = D+H (I, K) * VECR (K, IVEC - 1) 
   34    VECR (I, IVEC) = (VECR (I, IVEC) - D) / ETA 
   35 END DO 
!                                                                       
! THE COMPUTATION OF (INFIN. NORM OF W(S/1))**2.                        
      L = 1 
      S = 0.0D0 
      DO 36 I = 1, M 
         R = VECR (I, IVEC) * VECR (I, IVEC) + VECR (I, IVEC - 1)       &
         * VECR (I, IVEC - 1)                                           
         IF (R.LE.S) GOTO 36 
         S = R 
         L = I 
   36 END DO 
! THE COMPUTATION OF THE VECTOR Z(S/1),WHERE Z(S/1)=/W(S/1)/            
! (COMPONENT OF W(S/1) WITH THE LARGEST ABSOLUTE VALUE) .               
      U = VECR (L, IVEC - 1) 
      V = VECR (L, IVEC) 
      DO 37 I = 1, M 
         B = VECR (I, IVEC) 
         R = VECR (I, IVEC - 1) 
         VECR (I, IVEC) = (R * U + B * V) / S 
   37 VECR (I, IVEC - 1) = (B * U - R * V) / S 
! THE COMPUTATION OF THE RESIDUALS AND COMPARISON OF THE                
! RESIDUALS OF THE TWO SUCCESSIVE STEPS OF THE INVERSE                  
! GREATER THAN THE INFINITE NORM OF THE PREVIOUS RESIDUAL               
! ITERATION. IF THE INFINITE NORM OF THE RESIDUAL VECTOR IS             
! VECTOR THE COMPUTED VECTOR OF THE PREVIOUS STEP IS TAKEN              
! AS THE COMPUTED APPROXIMATION TO THE EIGENVECTOR.                     
      B = 0.0D0 
      DO 41 I = 1, M 
         R = WORK (I) * VECR (I, IVEC - 1) - ETA * VECR (I, IVEC) 
         U = WORK (I) * VECR (I, IVEC) + ETA * VECR (I, IVEC - 1) 
         IF (I.EQ.1) GOTO 38 
         R = R + SUBDIA (I - 1) * VECR (I - 1, IVEC - 1) 
         U = U + SUBDIA (I - 1) * VECR (I - 1, IVEC) 
   38    L = I + 1 
         IF (L.GT.M) GOTO 40 
         DO 39 J = L, M 
            R = R + H (I, J) * VECR (J, IVEC - 1) 
   39    U = U + H (I, J) * VECR (J, IVEC) 
   40    U = R * R + U * U 
         IF (B.GE.U) GOTO 41 
         B = U 
   41 END DO 
      IF (ITER.EQ.1) GOTO 42 
      IF (PREVIS.LE.B) GOTO 44 
   42 DO 43 I = 1, N 
         WORK1 (I) = VECR (I, IVEC) 
   43 WORK2 (I) = VECR (I, IVEC - 1) 
      PREVIS = B 
      IF (NS.EQ.1) GOTO 46 
      IF (ITER.GT.6) GOTO 47 
      ITER = ITER + 1 
      IF (BOUND.GT.DSQRT (S) ) GOTO 23 
      NS = 1 
      GOTO 23 
!                                                                       
   44 DO 45 I = 1, N 
         VECR (I, IVEC) = WORK1 (I) 
   45 VECR (I, IVEC - 1) = WORK2 (I) 
   46 INDIC (IVEC - 1) = 2 
      INDIC (IVEC) = 2 
   47 RETURN 
      END SUBROUTINE COMPVE                         
