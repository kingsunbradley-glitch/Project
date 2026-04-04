      SUBROUTINE FUSE(ELAB,AP,AT,ZP,ZT,SIGML,JDIM,IFUS,FLGRAZ)          00010011
C
C  Berechnet Fusionsquerschnitte, verschiedene Optionen, einschliess-
C  lich Bass-Modell,  modifiziert am 26.1.1980
C  Stand: 3.5.1994
C
C     PERFORMS SMOOTH CUT OFF OF UPPER PARTIAL WAVES (ENTRANCE)         00030000
C                                                                       00040000
C     IFUS=0  NO MODIFICATION OF SIGML                                  00050000
C     IFUS=1  FUSION GIVEN LCRIT AND DELTAL(DELTAL=0 IS SHARP CUTOFF)   00060000
C             WILL READ LCRIT,DELTAL                                    00070000
C     IFUS=2  READ RATIO SIGFUS TO SIG REACT                            00080000
C     IFUS=3  READ LCRIT                                                00090000
C     IFUS=4  BASS1 CUT OFF                                             00100000
C     IFUS=5  FUSION GIVEN SIGFUS  AND DELTAL                           00110001
C             (DELTAL=0 IS SHARP CUTOFF)                                00120001
C             WILL READ SIGFUS,DELTAL                                   00130001
C     IFUS>5  NO MODIFICATION                                           00140000
C                                                                       00150000
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIGML(JDIM)                                             00160011
      DIMENSION SIGFS (200) ,ARRAY(10)                                  00170004
C ------------------------------------------------------------------    00180011
C     MAYBE UNINITIALIZED: LMAX1,DSIG
C
C     I.Giese 2007
      SAVE
C
      IF(IFUS.LE.0 .OR. IFUS.GT.5) RETURN                               00190011
C --------------------------------------------------------------------  00200011
C                                                                       00210011
      IF(IFUS.EQ.2 .OR. IFUS .EQ. 3) THEN                               00220011
         IF(FLGRAZ .LT.1.) RETURN                                       00230011
         LGRAZ=FLGRAZ+0.5                                               00240011
         DO 30 I=1,JDIM                                                 00250011
         IF(SIGML(I))30,30,20                                           00260011
 20      LMAX=I-1                                                       00270011
 30      CONTINUE                                                       00280011
         WRITE(6,35) LMAX,FLGRAZ                                        00290011
 35      FORMAT(/' FUSE: VALUES FROM PARAP LMAX,LGRAZE',I10,F10.1)      00300011
         LMAX=MIN0(JDIM-1,LMAX)                                         00310011
         LMAX1=LMAX+1                                                   00320011
      ENDIF                                                             00330011
C --------------------------------------------------------------------  00340011
C                                                                       00350011
      ECM=ELAB*AT/(AT+AP)                                               00360011
      U=AT*AP/(AT+AP)                                                   00370011
      PWAVE=SQRT(2.*931.16*U*ECM)                                       00380011
      WAVEL=19.732/PWAVE                                                00390011
      ARWAVE=3141.59*(WAVEL**2.)                                        00400011
      GOTO(50,80,160,180,60),IFUS                                       00410011
C----------------------------------------------------------------       00420011
C     FUSION GIVEN LCRIT,DELTAL(WELL ABOVE CB)                 IFUS=1   00430011
C                                                                       00440011
 50   CALL MYRD(ARRAY,2,23,5)                                           00450011
      FLCRIT=ARRAY(1)                                                   00460011
      DELTAL=ARRAY(2)                                                   00470011
      LUP =FLCRIT+5.*DELTAL+1.5                                         00480011
      LUP =MIN(JDIM,LUP )                                               00490011
      SIGFUS=0.                                                         00500011
      DO 55 L=1,LUP                                                     00510011
      AL=L-1                                                            00520011
      TL=1.                                                             00530011
      IF(DELTAL.LE.0)GOTO 53                                            00540011
      DUM=(AL-FLCRIT)/DELTAL                                            00550011
      TL=1./(1.+EXP( DUM))                                              00560011
 53   SIGML(L)=ARWAVE*TL*(2.*AL+1.)                                     00570011
      SIGFUS=SIGFUS+SIGML(L)                                            00580011
 55   CONTINUE                                                          00590011
      IF (LUP.LT.JDIM) THEN                                             00600011
          DO 57 L=LUP+1,JDIM                                            00610011
 57       SIGML(L)=0.                                                   00620011
      ENDIF                                                             00630011
      LMAX=LUP-1                                                        00640011
      WRITE(6,56) ELAB,FLCRIT,DELTAL,LMAX,SIGFUS                        00650011
 56   FORMAT(' FUSE OPTION 1',3X,'ELAB=',F8.2,3X,'FLCRIT,DELTAL',2F8.2, 00660011
     1       3X,'LMAX=',I3,3X,'SIGFUS(MB)=',E12.4/)                     00670011
      RETURN                                                            00680011
C ------------------------------------------------------------------    00690011
C                             GIVEN SIGFUS,DELTAL               IFUS=5  00700011
 60   CALL MYRD(ARRAY,2,23,5)                                           00710011
      SIGFUS=ARRAY(1)                                                   00720004
      DELTAL=ARRAY(2)                                                   00730004
C                                                                       00740011
      FLCRIT=0.                                                         00750011
      DO 61 L=1,JDIM                                                    00760011
 61   SIGML(L)=0.                                                       00770011
      DO 62 I=1,3                                                       00780011
 62   FLCRIT=SQRT(SIGFUS/ARWAVE-FLCRIT)                                 00790011
      IREP=0                                                            00800001
 63   SIGF=0.                                                           00810011
      LUP =FLCRIT+5.*DELTAL+1.5                                         00820001
      LUP =MIN(JDIM,LUP )                                               00830011
      DO 65 L=1,LUP                                                     00840001
      AL=L-1                                                            00850001
      TL=1.                                                             00860001
      IF(DELTAL.GT.0.01) THEN                                           00870011
          DUM=(AL-FLCRIT)/DELTAL                                        00880011
          IF(DUM.GE.150.) TL=0.                                         00890011
          IF(ABS(DUM).LT.150.) TL=1./(1.+EXP( DUM))                     00900011
      ENDIF                                                             00910011
      SIGML(L)=ARWAVE*TL*(2.*AL+1.)                                     00920011
      SIGF=SIGF+SIGML(L)                                                00930001
 65   CONTINUE                                                          00940001
      IF(DELTAL.LE.0.) GOTO 68                                          00950003
      IF(IREP.GT.1) GOTO 74                                             00960003
      DSIG=SIGFUS-SIGF                                                  00970001
      DUM=0.01*SIGFUS                                                   00980001
      IF(ABS(DSIG).LT.DUM)GOTO 74                                       00990001
      IREP=IREP+1                                                       01000003
      DUM=     SIGF /(FLCRIT*FLCRIT)                                    01010001
      DLCRIT=DSIG/(2.*DUM*FLCRIT)                                       01020001
      FLCRIT=FLCRIT+DLCRIT                                              01030001
      GOTO 63                                                           01040011
 68   SIGF=SIGFUS-DSIG                                                  01050003
      IF(DSIG.GT.0.) GOTO 69                                            01060003
      SIGML(LUP)=SIGML(LUP)+DSIG                                        01070003
      GOTO 74                                                           01080003
 69   LUP=LUP+1                                                         01090003
      SIGML(LUP)=DSIG                                                   01100003
      FLCRIT=LUP-1                                                      01110003
 74   LMAX=LUP-1                                                        01120011
      WRITE(6,75) ELAB,SIGFUS,DELTAL,SIGF,FLCRIT,LMAX                   01130011
 75   FORMAT(/' FUSE OPTION 5  ELAB(MEV)=',F8.2/                        01140011
     1       ' ENTERED SIGFUS,DELTAL:',E12.4,F8.3,5X,                   01150012
     2       ' OUT SIGFUS,FLCRIT:',E12.4,F8.2,3X,'LMAX=',I3)            01160011
      RETURN                                                            01170011
C ------------------------------------------------------------------    01180011
C            GIVEN RATIO( 0 T0 1) COMPARED TO PARAP         IFUS=2      01190011
 80   CALL MYRD(ARRAY,1,23,5)                                           01200009
      RATIO=ARRAY(1)                                                    01210004
 81   FORMAT( 9F8.3)                                                    01220000
 90   IF(RATIO.LE.0. .OR. RATIO .GE.1.) RETURN                          01230011
      FLCRIT=SQRT(RATIO)*FLGRAZ                                         01240000
      GOTO 195                                                          01250011
C ------------------------------------------------------------------    01260011
C            SMALLER LCRIT        COMPARED TO PARAP         IFUS=3      01270011
 160  CALL MYRD(ARRAY,1,23,5)                                           01280011
      FLCRIT=ARRAY(1)                                                   01290004
      GOTO 195                                                          01300011
C-------------------------------------------------------------------    01310011
C     BASS MODEL  NPA231(1974)45  EQUATION NR IN COL 72-80  IFUS=4      01320011
C                                                                       01330000
 180  AS=17.                                                            01340000
      RZERO=1.07                                                        01350000
      FF=5./7.                                                          01360000
      D=1.35                                                            01370000
      AP3=AP**0.33333                                                   01380000
      AT3=AT**0.33333                                                   01390000
      SUMA3=AP3+AT3                                                     01400000
      R12=RZERO*SUMA3                                                   01410000
      DUM=AP3*AT3*(AP3+AT3)                                             01420000
      X=0.0792*ZP*ZT/DUM                                                01430000
      DUM=SUMA3*SUMA3*AP*AT*AT3*AP3                                     01440000
      Y=1.065*(AP+AT)/DUM                                               01450000
      VR=AS*AP3*AT3*D/R12                                               01460000
      F=VR                                                              01470000
      VCOUL=1.44*ZT*ZP/R12                                              01480000
      E1=VCOUL*(1.+0.5*(1.-X)/X-D/(X*R12))                              01490000
      E1LAB=E1*(AT+AP)/AT                                               01500000
      CFUE1=SQRT((1.-X)/(2.*Y))                                         01510000
      E2=VCOUL*(1.+0.5*(1.-X)/(X*FF*FF)-D/(X*R12))                      01520000
      E2LAB=E2*(AT+AP)/AT                                               01530000
      CFUE2=CFUE1/FF                                                    01540000
      PRINT 170                                                         01550000
 170  FORMAT(/ 2X,'BASS MODEL' )                                        01560000
      PRINT 181,E1,E1LAB,CFUE1,E2,E2LAB,CFUE2                           01570000
 181  FORMAT(1X ,'E1,E1LAB,LFUE1',3F10.3/1H0,'E2,E2LAB,LFUE2',3F10.3 )  01580000
C                                                                       01590000
C     FOR ECM.LT.E1 USE INVERTED PARABOLA(INVPAR) FOR LCRIT             01600000
C     FOR E1.LT.ECM.LT.E2   LCRIT IS LINEAR INTERPOLATION BETWEEN CFUE1 01610000
C     FOR ECM.GT.E2  LCRIT=CFUE2=CONSTANT                               01620000
      IF(ELAB-E1LAB)183,183,185                                         01630000
 183  CALL INVPAR(ELAB,F,RZERO,D,AP,ZP,AT,ZT,SIGFS ,FLCRIT)             01640002
      DO 184 L=1,200                                                    01650000
 184  SIGML(L)=0.                                                       01660000
      DO 187 L=1,199                                                    01670000
      SIGML(L)=SIGFS (L)                                                01680002
      IF(SIGFS (L+1).LE..1E-3) GOTO 188                                 01690002
 187  CONTINUE                                                          01700000
 188  LUP=L                                                             01710000
      LUP1=L+1                                                          01720000
      GOTO 218                                                          01730000
 185  IF(ELAB-E2LAB) 186,186,190                                        01740000
 186  ECM=ELAB*AT/(AT+AP)                                               01750000
      FLCRIT=CFUE1*CFUE1      +(CFUE2*CFUE2-CFUE1*CFUE1)*(ECM-E1)       01760000
     1 /(E2-E1)                                                         01770000
      FLCRIT=SQRT(FLCRIT)                                               01780000
      GOTO 195                                                          01790000
 190  FLCRIT=CFUE2                                                      01800000
 195  LCRIT=FLCRIT+0.5                                                  01810000
 200  FDL=FLGRAZ-FLCRIT                                                 01820000
      IDL=FDL                                                           01830000
      IF(FDL.GT.0.) GOTO 201                                            01840000
      IDL=-FDL+1.                                                       01850000
      IDL=-IDL                                                          01860000
 201  IF(IFUS.EQ.3) GOTO 202                                            01870000
      IF(FDL.LE.0.) GOTO 400                                            01880000
 202  CONTINUE                                                          01890000
C                                                                       01900000
C     SMOOTH CUT OFF FUSION PARTIAL X-SECTIONS                          01910000
C                                                                       01920000
      DO 205 L=1,LMAX1,10                                               01930000
      L1=MIN0(L+9,LMAX)                                                 01940000
      WRITE(6,225) (SIGML(L2),L2=L,L1)                                  01950000
 205  CONTINUE                                                          01960000
C                                                                       01970000
      LUP=LMAX-IDL                                                      01980000
      LUP1=LUP+1                                                        01990000
      LDN1=MAX0(1,1-IDL)                                                02000000
      DO 210 L=LDN1,LUP1                                                02010000
      AL=L                                                              02020000
      L1=AL+FDL                                                         02030000
      FL1=L1                                                            02040000
      X=AL+FDL-FL1                                                      02050000
      FL=L1+L1-1                                                        02060000
      FL1=L1+L1+1                                                       02070000
      DUM=(1.-X)*SIGML(L1)/FL+X*SIGML(L1+1)/FL1                         02080000
      FL=L+L-1                                                          02090000
      SIGFS (L)=DUM*FL                                                  02100002
 210  CONTINUE                                                          02110000
      DO 213 L=LDN1,LUP1                                                02120000
 213  SIGML(L)=SIGFS (L)                                                02130002
        LUP1=LUP1+1                                                     02140000
      IF(LUP1.GT.LMAX1) GOTO 217                                        02150000
 214  DO 215 L=LUP1,LMAX1                                               02160000
 215  SIGML(L)=0.                                                       02170000
 217  CONTINUE                                                          02180000
C                                                                       02190000
 218  CONTINUE                                                          02200000
      WRITE(6,220) LMAX,FLGRAZ,LUP,FLCRIT,ELAB,RATIO                    02210000
 220  FORMAT(//2X,'LMAX,LGRAZ,LUP,LCRIT',I8,F8.1,I8,F8.1/               02220000
     1 1H0,'ELAB,RATIO',2F10.3/2X,'FUSION PARTIAL X-SECTIONS')          02230000
      DO 230 L=1,LUP1,10                                                02240000
      L1=MIN0(L+9,LUP1)                                                 02250000
      WRITE(6,225) (SIGML(L2),L2=L,L1)                                  02260000
 225  FORMAT(1X,10F11.3)                                                02270008
      FLGRAZ=FLCRIT                                                     02280000
 230  CONTINUE                                                          02290000
 400  CONTINUE                                                          02300000
      RETURN                                                            02310000
      END                                                               02320000
