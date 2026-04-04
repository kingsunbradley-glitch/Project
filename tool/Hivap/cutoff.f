      SUBROUTINE CUTOFF(N1,M1,POP,SPECJ,TOTAL,NSTEP,AMX,CUT,PRCNT,      00020002
     1                  ICUT,TRIM,IPOPS,EXMAX1,ALEVEL,IPRINT)           00030008
C
C  Eliminiert vernachlaessigbar kleine Elemente in der E-J-Bevoelkerung
C  Stand: 3.5.1994
C  
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),SPECJ(1),TOTAL(1)                                00050000
C                                                                       00060000
C------------------------------------------------------------------     00070012
C     MAYBE UNINITIALIZED: SUM1
C
C     I.Giese 2007
      SAVE
C
      TEST=TOTAL(NSTEP)*CUT                                             00080000
      TRIMJ=0.                                                          00090012
      TRIME=0.                                                          00100012
      TRIMEL=0.                                                         00110016
C------------------------------------------------------------------     00120012
C                                        TRIM OFF J STRIPE              00130012
  200 SUM=0.                                                            00140000
      DO210J1=1,N1                                                      00150000
      INDEX=N1+1-J1                                                     00160000
      SUM=SUM+SPECJ(INDEX)                                              00170000
      IF(SUM-TEST)210,220,220                                           00180000
  210 CONTINUE                                                          00190000
  220 N11=N1-J1+1                                                       00200000
      IF(N1.EQ.N11)GOTO235                                              00210000
      K1=N11*M1+1                                                       00220000
      KLAST=N1*M1                                                       00230000
      DO 230  K=K1,KLAST                                                00240000
  230 POP(K)=0.                                                         00250000
      N1=N11                                                            00260000
      TRIMJ=SUM-SPECJ(N1)                                               00270006
C-------------------------------------------------------------------    00280012
C                                        TRIM OFF E STRIPE              00290012
C                                                                       00300000
  235 SUM=0.                                                            00310000
      DO 260 KE=1,M1                                                    00320000
      SUM1=0.                                                           00330000
      INDEX=KE-M1                                                       00340003
      DO 250  J1=1,N1                                                   00350000
      INDEX=INDEX+M1                                                    00360003
  250 SUM1=SUM1+POP(INDEX)                                              00370000
      SUM=SUM+SUM1                                                      00380000
      IF(SUM.GE.TEST) GO TO 270                                         00390000
  260 CONTINUE                                                          00400000
  270 M11=M1-KE+1                                                       00410000
      M11=MAX0(M11,1)                                                   00420009
      IF(M11.EQ.M1) GO TO 300                                           00430000
      TRIME=    SUM-SUM1                                                00440006
C----------------------  REARRANGE POP FOR NEW VALUE OF M1 ----------   00450012
      I=KE                                                              00460003
      IROW=-M11                                                         00470003
      ITROW=-M1                                                         00480003
      DO 290 J1=1,N1                                                    00490000
      IROW= IROW+M11                                                    00500003
      ITROW= ITROW+M1                                                   00510003
      DO 285 KE=1,M11                                                   00520000
      INDEX=IROW+KE                                                     00530000
      ITM=ITROW+I-1+KE                                                  00540000
      POP(INDEX)=POP(ITM)                                               00550000
  285 CONTINUE                                                          00560000
  290 CONTINUE                                                          00570000
      EXMAX1=EXMAX1-FLOAT(I-1)                                          00580000
      M1=M11                                                            00590000
C---------------------------------------------------------------------  00600012
C                                TRIM OFF INDIVIDUAL SMALL ELEMENTS     00610012
C                                                                       00620000
 300  ALEVEL=0.                                                         00630015
      IF(AMX.LE.0.) GOTO 350                                            00640015
      IAGAIN=0                                                          00650003
      TRIMEL=0.                                                         00660006
      IF(PRCNT.GE.1.) GOTO 350                                          00670006
      IF(PRCNT.LE..5) GOTO 350                                          00680006
      ALEVEL=AMX                                                        00690015
      FRACTN=PRCNT*(TOTAL(NSTEP)-TRIME-TRIMJ)                           00700006
C                                                                       00710000
 310  ALEVEL=ALEVEL*0.5                                                 00720000
      IAGAIN=IAGAIN+1                                                   00730000
      INDEX=-M1                                                         00740000
      SUM=0.                                                            00750000
C                                                                       00760000
      DO 330 J1=1,N1                                                    00770000
      INDEX=INDEX+M1                                                    00780000
      DO 320 KE=1,M1                                                    00790000
      DUM=POP(INDEX+KE)                                                 00800000
      IF(DUM.GT.ALEVEL)SUM=SUM+DUM                                      00810000
 320  CONTINUE                                                          00820000
 330  CONTINUE                                                          00830000
C                                                                       00840000
      IF(IAGAIN.GT.20) GOTO 350                                         00850000
      IF(SUM.LT.FRACTN) GOTO 310                                        00860000
C                                                                       00870010
      ICUT=0                                                            00880010
      IF(M1.LT.2) GOTO 350                                              00890010
      INDEX=-M1                                                         00900003
      DO 340 J1=1,N1                                                    00910000
      INDEX=INDEX+M1                                                    00920000
      DO 335 KE=1,M1                                                    00930000
      I=INDEX+KE                                                        00940000
      IF(POP(I).GT.ALEVEL) GOTO 335                                     00950000
      IF(POP(I).GT.0) ICUT=ICUT+1                                       00960003
      TRIMEL=TRIMEL+POP(I)                                              00970006
      POP(I)=0.                                                         00980000
 335  CONTINUE                                                          00990000
 340  CONTINUE                                                          01000000
 350  TRIM=TRIM+TRIMJ+TRIME+TRIMEL                                      01010006
      IF(IPRINT.LT.3)                                                   01020008
     1WRITE(6,351) NSTEP,TOTAL(NSTEP),TRIM,TRIMJ,TRIME,TRIMEL,ICUT,AMX, 01030008
     2             ALEVEL,IAGAIN ,CUT,PRCNT,M1,N1                       01040012
 351  FORMAT(' CUTOFF   NSTEP,TOTAL,TRIM,TRIMJ,TRIME,TRIMLEV',          01050006
     1 I4,5E12.4/10X,'ICUT,AMX,ALEVEL,IAGAIN',I4,2E12.4,I4,             01060010
     2 4X,'CUT,PRCNT',2E11.3,'  M1,N1=',2I4)                            01070012
C--------------------------------------------------------------------   01080012
C                     ADJUST M1 SO IT EXTENDS AS FAR DOWN AS POSSIBLE   01090012
C                           REARRANGE POP FOR NEW VALUES OF M1 AND N1   01100012
C                                                                       01110000
  400 IA=EXMAX1+1.                                                      01120000
      M11=MIN0(IA,IPOPS/N1)                                             01130000
      M11=MAX0(M11,1)                                                   01140009
      IF(M1.EQ.M11)GOTO 500                                             01150003
      IF(M1.GT.M11)GOTO440                                              01160003
C                                               M11>M1                  01170012
      IDIFF=M11-M1                                                      01180003
      MN=M11*N1                                                         01190003
      MN1=M1*N1                                                         01200003
 410  J=MN+1                                                            01210003
      K=MN1+1                                                           01220003
      DO 420 I=1,IDIFF                                                  01230003
      J=J-1                                                             01240003
 420  POP(J)=0.                                                         01250003
      IF(J.LE.(M1+1)) GOTO 430                                          01260003
      DO425 I=1,M1                                                      01270004
      J=J-1                                                             01280003
      K=K-1                                                             01290003
 425  POP(J)=POP(K)                                                     01300004
      MN=MN-M11                                                         01310003
      MN1=MN1-M1                                                        01320003
      GOTO 410                                                          01330003
 430  IF(IPRINT.LT.3)                                                   01340008
     1WRITE(6,435) M1,M11                                               01350008
 435  FORMAT(' CUTOFF OLD,NEW M1',2I4)                                  01360003
      M1=M11                                                            01370003
      RETURN                                                            01380003
C                                                M1>M11                 01390012
 440  IDIFF=M1-M11                                                      01400003
      SUM=0.                                                            01410003
      J=0                                                               01420003
      K=-M11                                                            01430003
      MN=M11*N1                                                         01440003
      MN1=M1*N1                                                         01450003
 450  J=J+M11                                                           01460003
      K=K+M11                                                           01470003
      DO 460 I=1,IDIFF                                                  01480003
      J=J+1                                                             01490003
      IF(J.GT.IPOPS) GOTO 460                                           01500003
      SUM=SUM+POP(J)                                                    01510003
 460  CONTINUE                                                          01520003
      IF(J.GE.MN1) GOTO 480                                             01530003
      DO 470 I=1,M11                                                    01540003
      J=J+1                                                             01550003
      K=K+1                                                             01560003
      IF(J.GT.IPOPS) GOTO 470                                           01570003
      POP(K)=POP(J)                                                     01580003
 470  CONTINUE                                                          01590003
      GOTO 450                                                          01600003
 480  M1=M11                                                            01610003
      TRIM=TRIM+SUM                                                     01620000
 500  CONTINUE                                                          01630003
      RETURN                                                            01640000
      END                                                               01650000
