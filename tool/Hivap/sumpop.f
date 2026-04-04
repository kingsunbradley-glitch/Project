      SUBROUTINE SUMPOP(N,M,EXMAX,POP,TOT,ISIZE,N1,M1,EXMAX1,POP1,TOT1, 00020003
     1                  POP0,FLOST,IPRINT)                              00030011
C  Addiert Bevoelkerungswahrscheinlichkeiten aus verschiedenen
C  Mutterkernen, die ueber xn, pxn, alphaxn in den gleichen Restkern
C  (Tochterkern) zerfallen
C  Stand: 3.5.1994
C 
C     ADDS POP AND POP1 INTO POP0,NEW DIMENSIONS INTO N1,M1,EXMAX1      00040001
C                   LOST PARTS INTO FLOST                               00050001
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),POP1(1),POP0(1),SUME(120),SUMJ(400)              00060018
C
C     I.Giese 2007
      SAVE
C                                                                       00070000
      GLOST=0.                                                          00080017
C
C     I.Giese 2007
C
      SUM=0.0D0
      IF(TOT .LE.0. .OR. M .EQ.0) GOTO 300                              00090006
      IF(TOT1.LE.0. .OR. M1.EQ.0) GOTO 400                              00100006
C                                                                       00110005
C                            GET ENCOMPASSING DIMENSIONS M0,N0,EXMAX0   00120005
      N0=MAX0(N,N1)                                                     00130000
      EXMAX0=DMAX1(EXMAX,EXMAX1)                                        00140000
      EXMIN=EXMAX-M+1                                                   00150000
      EXMIN1=EXMAX1-M1+1                                                00160000
      EXMIN0=DMIN1(EXMIN,EXMIN1)                                        00170000
      M0=EXMAX0-EXMIN0+1.                                               00180000
      MN0=M0*N0                                                         00190000
      MN=M*N                                                            00200014
      MN1=M1*N1                                                         00210014
      GLOST=0.                                                          00220014
      IF(MN0.LE.ISIZE) GOTO 150                                         00230000
C                                  GET SUME(J) OF (POP+POP1)            00240005
      JM=-M                                                             00250005
      JM1=-M1                                                           00260005
      DO 40 J=1,N0                                                      00270000
      SUM=0.                                                            00280000
      IF(J.GT.N) GOTO 25                                                00290000
      JM=JM+M                                                           00300005
      DO 20 K=1,M                                                       00310000
 20   SUM=SUM+POP(JM+K)                                                 00320000
 25   IF(J.GT.N1) GOTO 35                                               00330000
      JM1=JM1+M1                                                        00340005
      DO 30 K=1,M1                                                      00350000
 30   SUM=SUM+POP1(JM1+K)                                               00360000
 35   SUME(J)=SUM                                                       00370001
 40   CONTINUE                                                          00380001
C                                GET SUMJ(KE) OF (POP+POP1)             00390005
      IK=EXMAX0-EXMAX   +0.5                                            00400010
      IK1=EXMAX0-EXMAX1 +0.5                                            00410010
      DO 70 K=1,M0                                                      00420004
      SUM=0.                                                            00430001
      JM=K-IK                                                           00440005
      IF(JM.LE.0) GOTO 55                                               00450005
      IF(JM.GT.M) GOTO 55                                               00460005
      DO 50 J=1,N                                                       00470005
      SUM=SUM+POP(JM)                                                   00480005
 50   JM=JM+M                                                           00490005
 55   JM1=K-IK1                                                         00500005
      IF(JM1.LE.0) GOTO 65                                              00510005
      IF(JM1.GT.M1) GOTO 65                                             00520005
      DO 60 J=1,N1                                                      00530005
      SUM=SUM+POP1(JM1)                                                 00540005
 60   JM1=JM1+M1                                                        00550005
 65   SUMJ(K)=SUM                                                       00560001
 70   CONTINUE                                                          00570001
C                                                                       00580001
      I=1                                                               00590001
 80   SUM=DMIN1(SUME(N0),SUMJ(I),SUMJ(M0))                              00600001
 85   IF(SUME(N0).GT.SUM) GOTO 90                                       00610001
      N0=N0-1                                                           00620001
      MN0=M0*N0                                                         00630012
      IF(MN0.LE.ISIZE) GOTO 150                                         00640001
      GOTO 80                                                           00650001
 90   IF(SUMJ(M0).GT.SUM) GOTO 100                                      00660001
      M0=M0-1                                                           00670001
      MN0=M0*N0                                                         00680012
      IF(MN0.LE.ISIZE) GOTO 150                                         00690001
      GOTO 80                                                           00700001
 100  IF(SUMJ(I).GT.SUM) GOTO 85                                        00710001
      DO 105 K=2,M0                                                     00720012
 105  SUMJ(K-1)=SUMJ(K)                                                 00730013
      M0=M0-1                                                           00740012
      EXMAX0=EXMAX0-1                                                   00750001
      MN0=M0*N0                                                         00760012
      IF(MN0.LE.ISIZE) GOTO 150                                         00770001
      GOTO 80                                                           00780001
C---------------------------------------------------------------------  00790010
C                                     SUM POP+POP1 INTO POP0            00800010
 150  M0=MAX0(M0,1)                                                     00810015
      N0=MAX0(N0,1)                                                     00820015
      DUM=EXMAX0-EXMAX+0.5                                              00830015
      IK=DUM                                                            00840008
      IK=-IK                                                            00850008
      DUM=EXMAX0-EXMAX1+0.5                                             00860008
      IK1=DUM                                                           00870008
      IK1=-IK1                                                          00880008
      SUM=0.                                                            00890001
      JM0=0                                                             00900005
      JM1=0                                                             00910005
      JM=0                                                              00920005
      DO 200 J0=1,N0                                                    00930001
      J=J0-1                                                            00940001
      DO 180 K0=1,M0                                                    00950001
      DUM1=0.                                                           00960001
      K1=K0+IK1                                                         00970001
      IF(K1.LE.0) GOTO 170                                              00980001
      IF(K1.GT.M1)GOTO 170                                              00990014
      JK1=JM1+K1                                                        01000001
      IF(JK1.GT.MN1  ) GOTO 170                                         01010014
      DUM1=POP1(JK1)                                                    01020001
 170  DUM=0.                                                            01030001
      K=K0+IK                                                           01040001
      IF(K.LE.0) GOTO 175                                               01050001
      IF(K.GT.M) GOTO 175                                               01060014
      JK=JM+K                                                           01070001
      IF(JK.GT.MN   ) GOTO 175                                          01080014
      DUM=POP(JK)                                                       01090001
 175  JK0=JM0+K0                                                        01100001
      IF(JK0.GT.ISIZE)GOTO 180                                          01110014
      POP0(JK0)=DUM+DUM1                                                01120001
      SUM=SUM+POP0(JK0)                                                 01130001
 180  CONTINUE                                                          01140001
      JM0=JM0+M0                                                        01150005
      JM1=JM1+M1                                                        01160005
      JM=JM+M                                                           01170005
 200  CONTINUE                                                          01180001
C-------------------------------------------------------------------    01190010
      GLOST=TOT1+TOT-SUM                                                01200014
      IF(IPRINT.LT.3)                                                   01210016
     1WRITE(6,210) GLOST,SUM,EXMAX0,M0,N0,TOT,EXMAX,M,N,TOT1,EXMAX1,M1, 01220014
     2             N1,IK,IK1                                            01230011
 210  FORMAT(' SUMPOP FLOST,SUM,EXMAX,M,N',E12.4,3(E12.4,F8.2,2I4),2I3) 01240008
      TOT1=SUM                                                          01250001
      N1=N0                                                             01260001
      M1=M0                                                             01270001
      EXMAX1=EXMAX0                                                     01280001
      TOT1=SUM                                                          01290002
 215  I1=M1*N1+1                                                        01300014
      DO 220 I=I1,ISIZE                                                 01310014
 220  POP0(I)=0.                                                        01320014
      FLOST=FLOST+GLOST                                                 01330014
      RETURN                                                            01340001
C                                                                       01350005
 300  MN1=M1*N1                                                         01360014
      DO 310 I=1,MN1                                                    01370014
 310  POP0(I)=POP1(I)                                                   01380005
      IF(IPRINT.LT.3)                                                   01390017
     1WRITE(6,210) GLOST,TOT1,EXMAX1,M1,N1                              01400014
      GOTO 215                                                          01410014
C                                                                       01420005
 400  MN=M*N                                                            01430005
      DO 410  I=1,MN                                                    01440005
 410  POP0(I)=POP(I)                                                    01450005
      N1=N                                                              01460005
      M1=M                                                              01470005
      EXMAX1=EXMAX                                                      01480005
      TOT1=TOT                                                          01490012
      IF(IPRINT.LT.3)                                                   01500017
     1WRITE(6,210) GLOST,TOT1,EXMAX1,M1,N1                              01510015
      GOTO 215                                                          01520014
      END                                                               01530001
