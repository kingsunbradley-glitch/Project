      SUBROUTINE   ROT (NSTEP,NROW,AMASS,ZEE,SHELK,ITEST)               00020000
C
C  Liest bzw. berechnet Yrastlinien (auch am Sattelpunkt),
C  liest Yrastlinien von logischer Einheit 'IYR' wenn IRAT=1;
C  modifiziert fuer Sierk-Barrieren bei schweren Elementen (8.2.1989)
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/YR/JDIM,EJAY(500),JUPYR,IYR,IRAST,RATIOS(6),RATIOF(6),     00050000
     1          INERF,FINERT,FINERF,NOROTF                              00060000
      COMMON/BARR/BFLDM(25,25),BFS(25,25),BAR0,BARFAC,IBF,IFISRT        00070000
      COMMON/DNS/YRSMO(500),CS,CE,UCRIT,RPAGAP,IASYM,IASYMS,ENHANC(400),00080000
     1           ENHANS(400)                                            00090000
      COMMON/FISRT/ROT0,A2MS,RKAPPA,X0,Y0,RKAPA4                        00100003
      DIMENSION ARRAY(500),SHELK(5)                                     00110009
      CHARACTER*1 AW(72)                                                00120008
C     MAYBE UNINITIALIZED: DUMF,DUMG,NOF
C
C     I.Giese 2007
      SAVE
C
      EPS=0.01                                                          00130000
      AMPAR1=0.
C---------------------------------------------------------------------- 00140000
C             'IRAST' NOT ZERO: YRAST FROM DISC LOG UNIT IYR            00150000
      IF(IRAST   .EQ.0) GOTO 150                                        00160000
      ITEST=0                                                           00170000
      J0=JDIM                                                           00180000
      JDIM1=JDIM-1                                                      00190000
      IRC=0                                                             00200000
      CALL DECODN(ARRAY,AW,IRC,5,72)                                    00210008
      NJ =0                                                             00220000
      IF(IRC.GT.0) NJ =ARRAY(1)+EPS                                     00230000
      NJ1=0                                                             00240000
      IF(IRC.GT.1) NJ1=ARRAY(2)+EPS                                     00250000
      NJ2=0                                                             00260000
      IF(IRC.GT.2) NJ2=ARRAY(3)+EPS                                     00270000
      NJ3=0                                                             00280000
      IF(IRC.GT.3) NJ3=ARRAY(4)+EPS                                     00290000
      NJ4=0                                                             00300000
      IF(IRC.GT.4) NJ4=ARRAY(5)+EPS                                     00310000
C                                                                       00320000
      IF(NJ .GT.0) CALL MYRD(EJAY   ,NJ,-1,IYR)                         00330000
      DO 40 K=NJ ,JDIM1                                                 00340000
 40   EJAY(K+1)=9999.                                                   00350000
C                                                                       00360000
      IF(NJ1.GT.0) CALL MYRD(EJAY(J0+1),NJ1,-1,IYR)                     00370000
      DO 41 K=NJ1,JDIM1                                                 00380000
      KK=K+1+J0                                                         00390000
 41   EJAY(  KK )=9999.                                                 00400000
      J0=J0+JDIM                                                        00410000
C                                                                       00420000
      IF(NJ2.GT.0) CALL MYRD(EJAY(J0+1),NJ2,-1,IYR)                     00430000
      DO 42 K=NJ2,JDIM1                                                 00440000
      KK=K+1+J0                                                         00450000
 42   EJAY(  KK )=9999.                                                 00460000
      J0=J0+JDIM                                                        00470000
C                                                                       00480000
      IF(NJ3.GT.0) CALL MYRD(EJAY(J0+1),NJ3,-1,IYR)                     00490000
      DO 43 K=NJ3,JDIM1                                                 00500000
      KK=K+1+J0                                                         00510000
 43   EJAY(  KK )=9999.                                                 00520000
      IF(NOF.EQ.1) GOTO 150                                             00530002
      J0=J0+JDIM                                                        00540000
C                                                                       00550000
      IF(NJ4.GT.0) CALL MYRD(EJAY(J0+1),NJ4,-1,IYR)                     00560000
      DO 44 K=NJ4,JDIM1                                                 00570000
      KK=K+1+J0                                                         00580000
 44   EJAY(  KK )=9999.                                                 00590000
C-----------------------------------------------------------------------00600000
C                               LIQUID DROP ROTOR (CPS 1974)            00610000
 150  IF(FINERT.LE.0.) FINERT=1.                                        00620000
      IF(FINERF.LE.0.) FINERF=1.                                        00630000
      DUM=0.                                                            00640000
      IF(JUPYR.LE.0) JUPYR=JDIM-1                                       00650000
      AN=AMASS-ZEE                                                      00660000
      EJAY(1)=0.                                                        00670000
      YRSMO(1)=0.                                                       00680000
      J0F=4*JDIM                                                        00690000
C----------------------------------------------------                   00700000
C     FISSION BARRIER FOR ZERO AM                                       00710000
      BF0=BFS(NSTEP,NROW)                                               00720000
      IF(IRAST.EQ.0) EJAY(J0F+1)=BF0                                    00730000
      SPT0=BFLDM(NSTEP,NROW)                                            00740000
      DM0=SPT0*BARFAC                                                   00750000
      DBF0=SPT0 * (BARFAC-1.)                                           00760000
      YRSMO(J0F+1)=DM0                                                  00770000
C                                                                       00780000
C------------------------------------------------------                 00790000
C                          YRAST OF GAMMA DAUGHTER AND SADDLE POINT     00800000
      ALIM=RATIOS(6)+16.                                                00810000
      IRAT =RATIOS(6)+1.01                                              00820000
      DIFF=0.                                                           00830000
      SHELK0=SHELK(1)                                                   00840000
      JUPYR1=MAX0(JUPYR,31)                                             00850006
C
C     I.Giese 2007
C     0.0D0
C
      CALL FISROT(AMASS,ZEE,AN,0.0D0,DELR,BFCPS,EROT,0.0D0)             00860011
      DO 160 J=2,JUPYR1                                                 00870006
      AA=J-1                                                            00880000
      CALL FISROT(AMASS,ZEE,AN,AA,DELR,SPT,EROT,DUM)                    00890006
      IF(IFISRT.EQ.2) THEN                                              00900011
         IZEE=ZEE+0.01                                                  00910011
         MASS=AMASS+0.01                                                00920011
         IF(IZEE.LE.102) THEN                                           00930011
             CALL BARFIT(IZEE,MASS,J-1,SPT,DELR,ELMAX)                  00940011
             SPT=SPT+DELR                                               00950011
         ELSE                                                           00960011
             SPT=SPT-BFCPS +SPT0                                        00970011
         ENDIF                                                          00980011
      ENDIF                                                             00990011
C                                                                       01000011
      DELR=DELR*(AA+1.)/(AA*FINERT)                                     01010011
      SPT=(SPT-SPT0)*(AA+1.)/(AA*FINERF) + SPT0                         01020005
      IF(NOROTF.EQ.2) THEN                                              01030011
         IF(J.EQ.2) DUMG=DELR*0.5                                       01040011
         DELR=AA*(AA+1.)*DUMG                                           01050011
         IF(J.EQ.2)DUMF=(SPT-BFLDM(NSTEP,NROW))*0.5                     01060011
         SPT=AA*(AA+1.)*DUMF + BFLDM(NSTEP,NROW)                        01070011
      ENDIF                                                             01080011
      DM=(SPT-DELR)*BARFAC                                              01090011
      DM=DMAX1(DM,AMPAR1)                                               01100000
      YRSMO(J)=DELR                                                     01110000
      IF(AA.GT.RATIOS(6).OR. IRAST.NE.0) GOTO 154                       01120000
      CALL IP2(RATIOS(2),RATIOS(4),RATIOS(6),AA,RATIOS(1),RATIOS(3),    01130000
     1 RATIOS(5),EJAY(J))                                               01140000
      IF( J.EQ.IRAT     ) SHELK0=(EJAY(J)-YRSMO(J)+SHELK(1))*DM0/DM     01150000
      GOTO 155                                                          01160000
 154  IF(IRAST.EQ.0) EJAY(J)=YRSMO(J)-SHELK(1)+SHELK0*DM /DM0           01170000
 155  YRSMO(J0F+J)=DELR+DM                                              01180000
      IF(NOROTF.EQ.1) YRSMO(J0F+J)=SPT+DBF0                             01190000
      IF(IRAST.EQ.0)EJAY(J0F+J)=YRSMO(J0F+J)-SHELK(1)+SHELK(5)*DM/DM0   01200000
 160  CONTINUE                                                          01210000
C                                                                       01220000
C-----------------------------------------------------------------      01230000
C                                    YRAST OF N,P,A DAUGHTERS           01240000
      DUM=0.                                                            01250000
      K=1                                                               01260000
 161  K=K+1                                                             01270000
      GOTO(163,163,166,169,190),K                                       01280000
 163  AF=AMASS-1                                                        01290000
      ZF=ZEE                                                            01300000
      GOTO 175                                                          01310000
 166  ZF=ZEE-1                                                          01320000
      GOTO 175                                                          01330000
 169  AF=AMASS-4.                                                       01340000
      ZF=ZEE-2.                                                         01350000
 175  J0=(K-1)*JDIM                                                     01360000
      SHELK0=SHELK(K)                                                   01370000
      AN=AF-ZF                                                          01380000
      MASS=AF+0.01                                                      01390006
      IZEE=ZF+0.01                                                      01400006
      AA=0.                                                             01410000
      CALL FISROT(AF,ZF,AN,AA,DELR,DM0,EROT,0.0D0)                      01420011
      BFCPS=DM0                                                         01430011
      IF(IFISRT.EQ.2)                                                   01440006
     1CALL BARFIT(IZEE,MASS,0,DM0,DELR,ELMAX)                           01450006
      EJAY(J0+1)=0.                                                     01460000
      YRSMO(J0+1)=0.                                                    01470000
      DIFF=0.                                                           01480000
C---------------------------------                                      01490000
      DO 185 J=2,JUPYR1                                                 01500006
      AA=J-1                                                            01510000
      CALL FISROT(AF,ZF,AN,AA,DELR,SPT,EROT,DUM)                        01520011
      IF(IFISRT.EQ.2) THEN                                              01530011
         IF(IZEE.LE.102) THEN                                           01540011
             CALL BARFIT(IZEE,MASS,J-1,SPT,DELR,ELMAX)                  01550011
             SPT=SPT+DELR                                               01560011
         ELSE                                                           01570011
             SPT=SPT-BFCPS +DM0                                         01580011
         ENDIF                                                          01590011
      ENDIF                                                             01600011
C                                                                       01610011
      DELR=DELR*(AA+1.)/(AA*FINERT)                                     01620011
      SPT=(SPT-DM0)*(AA+1.)/(AA*FINERF) + DM0                           01630005
      IF(NOROTF.EQ.2) THEN                                              01640011
         IF(J.EQ.2) DUMG=DELR*0.5                                       01650011
         DELR=AA*(AA+1.)*DUMG                                           01660011
         IF(J.EQ.2)DUMF=(SPT-DM0)*0.5                                   01670011
         SPT=AA*(AA+1.)*DUMF + DM0                                      01680011
      ENDIF                                                             01690011
      DM =SPT-DELR                                                      01700011
      DM=DMAX1(AMPAR1,DM)                                               01710000
      YRSMO(J0+J)=DELR                                                  01720000
      IF(AA.GT.RATIOS(6).OR. IRAST.NE.0) GOTO 178                       01730000
      CALL IP2(RATIOS(2),RATIOS(4),RATIOS(6),AA,RATIOS(1),RATIOS(3),    01740000
     1 RATIOS(5),EJAY(J0+J))                                            01750000
      IF( J.EQ.IRAT     ) SHELK0=(EJAY(J0+J)-YRSMO(J0+J) +SHELK(K))*    01760000
     1                            DM0/DM                                01770000
      GOTO 185                                                          01780000
 178  IF(IRAST.EQ.0)EJAY(J0+J)=YRSMO(J0+J)-SHELK(K)+SHELK0*DM/DM0       01790000
 185  CONTINUE                                                          01800000
C----------------------------                                           01810000
      GOTO 161                                                          01820000
 190  IF(IRAST.NE.0) GOTO 999                                           01830000
      JUP1=JUPYR1+1                                                     01840006
      DO 205 K=1,5                                                      01850000
      J0=(K-1)*JDIM                                                     01860000
      DO 200 J=JUP1,JDIM                                                01870000
      JJ=J0+J                                                           01880000
 200  EJAY(JJ)=9999.                                                    01890000
 205  CONTINUE                                                          01900000
C----------------------------------------------------------------       01910000
C                                      TEST OUTPUT                      01920000
      IF(ITEST.EQ.0) GOTO 999                                           01930000
      ITEST=0                                                           01940000
      J1=1                                                              01950000
      J11=0                                                             01960000
      J21=1                                                             01970000
      WRITE(16,210) J1,ITEST,ITEST                                      01980004
 210  FORMAT(' ROT OUTPUT  FORMAT:',3I4/' EJAY,EJAYF')                  01990000
 215  FORMAT(15F8.2)                                                    02000000
 220  FORMAT(6X)                                                        02010000
      DO 240 K=1,5,4                                                    02020000
      J2=J1+JUPYR1-1                                                    02030006
      WRITE(16,230) (EJAY(J),J=J1,J2)                                   02040004
 230  FORMAT(5F10.4)                                                    02050000
      WRITE(16,220)                                                     02060004
      WRITE(16,221)J11,J21                                              02070004
 221  FORMAT(2I6/)                                                      02080000
      J1=J1+4*JDIM                                                      02090000
 240  CONTINUE                                                          02100000
      J0=4*JDIM                                                         02110000
      DO 235 J=1,JUPYR1                                                 02120006
 235  ARRAY(J)=EJAY  (J0+J)-EJAY  (J)                                   02130000
      WRITE(16,236)                                                     02140004
 236  FORMAT(' BARRIER')                                                02150000
      WRITE(16,230) (ARRAY(J),J=1,JUPYR1)                               02160006
      WRITE(16,220)                                                     02170004
      WRITE(16,221)J11,J21                                              02180004
      WRITE(16,241)                                                     02190004
 241  FORMAT(' ROT OUTPUT  YRSMO,YRSMOF')                               02200000
      J1=1                                                              02210000
      DO 250 K=1,5,4                                                    02220000
      J2=J1+JUPYR1-1                                                    02230006
      WRITE(16,230) (YRSMO(J),J=J1,J2)                                  02240004
      WRITE(16,220)                                                     02250004
      WRITE(16,221)J11,J21                                              02260004
 232  FORMAT(1X,10F10.4)                                                02270000
      J1=J1+4*JDIM                                                      02280000
 250  CONTINUE                                                          02290000
      J0=4*JDIM                                                         02300000
      DO 255 J=1,JUPYR1                                                 02310006
 255  ARRAY(J)=YRSMO (J0+J)-YRSMO (J)                                   02320000
      WRITE(16,236)                                                     02330004
      WRITE(16,230) (ARRAY(J),J=1,JUPYR1)                               02340006
      WRITE(16,220)                                                     02350004
      WRITE(16,221)J11,J21                                              02360004
      WRITE(16,220)                                                     02370004
 999  RETURN                                                            02380000
      END                                                               02390000
