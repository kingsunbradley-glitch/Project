      SUBROUTINE GRIDS(NROW,NSTEP,CUTPOP)                               00010034
C                                                                       00020048
C  Bereitet die E-J - Bevoelkerung vor der Hauptrechnung (in EVA) zu;
C  ruft GETPUT, SUMPOP, TOT, CUTOFF auf
C  Stand 3.5.1994
C
      INCLUDE 'common.f'                                                00030050
C     MAYBE UNINITIALIZED: N,M,EXMAX,SUMN,JYR
C
C     I.Giese 2007
      SAVE
C
      IPRINT=LPRINT                                                     00040040
      IF(NSTEP.EQ.1 .AND. LPRINT.LT.3) WRITE(6,5)                       00050036
 5    FORMAT(1H1,' NEW ROW')                                            00060034
      IF(LPRINT.LT.2)                                                   00070034
     1WRITE (6,10) NROW,NSTEP                                           00080034
 10   FORMAT(/' GRIDS  NROW,NSTEP',2I4)                                 00090034
C                                                                       00100034
C                                      GET P                            00110034
C                                                                       00120034
 250  MN=M2*N2                                                          00130034
      SUMN=0.                                                           00140034
      DO 260 I=1,MN                                                     00150034
 260  SUMN=SUMN+POPN(I)                                                 00160034
      IF(LPRINT.LT.2)                                                   00170034
     1WRITE(6,265) SUMN,SUMLOW                                          00180034
 265  FORMAT(' SUMN,SUMLOW',2E14.4)                                     00190034
      IDUM= STORE(NROW,NSTEP)                                           00200034
      IF(IDUM.EQ.0 .OR. IDUM.EQ.2) GOTO 275                             00210034
      IUNIT=10+IU1                                                      00220034
      IREAD=1                                                           00230034
      CALL GETPUT(POPP,IPOPS,IUNIT,IREAD,M,N,EXMAX,MMZ,MMA,AACN,ZZCN)   00240034
      MN=N*M                                                            00250034
      SUMP=0.                                                           00260034
      DO 270 I=1,MN                                                     00270034
 270  SUMP=SUMP+POPP(I)                                                 00280034
      I=AACN+0.001                                                      00290034
      II=ZZCN+0.001                                                     00300034
      IF(LPRINT.LT.2)                                                   00310034
     1WRITE(6,272) IUNIT,SUMP,M,N,EXMAX,MMA,ELEMNT(MMZ),I,ELEMNT(II)    00320034
 272  FORMAT(' GOT POPP FROM UNIT',I4,6X,'SUMP',E12.4,3X,'M,N,EXMAX',   00330034
     1        2I4,F8.1,4X,I4,A4,4X,'CN=',I4,A4)                         00340034
      GOTO 280                                                          00350034
 275  DO 276 I=1,IPOPS                                                  00360034
 276  POPP(I)=0.                                                        00370034
      SUMP=0.                                                           00380034
 280  DO 282 I=1,IPOPS                                                  00390034
 282  POPA(I)=0.                                                        00400034
      CALL SUMPOP(N ,M ,EXMAX ,POPP,SUMP,IPOPS,N2,M2,EXMAX2,POPN,       00410034
     1                 SUMN  ,POPA,FLOST(NSTEP),IPRINT)                 00420036
C                                                                       00430034
C                                      GET A                            00440034
C                                                                       00450034
 300  IDUM= STORE(NROW,NSTEP)                                           00460034
      IF(IDUM.EQ.0 .OR. IDUM.EQ.1) GOTO 350                             00470034
      IDIM=-1                                                           00480034
      IUNIT=12+IU2                                                      00490034
      K=4                                                               00500034
      IREAD=1                                                           00510034
      CALL GETPUT(POPN,IPOPS,IUNIT,IREAD,M,N,EXMAX,MMZ,MMA,AACN,ZZCN)   00520034
      MN=M*N                                                            00530034
      SUMA=0.                                                           00540034
      DO 310 I=1,MN                                                     00550034
 310  SUMA=SUMA+POPN(I)                                                 00560034
      IF(LPRINT.LT.2)                                                   00570034
     1WRITE(6,312) IUNIT,SUMA,M,N,EXMAX,MMA,ELEMNT(MMZ)                 00580034
 312  FORMAT(' GOT POPA FROM UNIT',I4,6X,'SUMA',E12.4,3X,'M,N,EXMAX',   00590034
     1        2I4,F8.1,4X,I4,A4)                                        00600034
      GOTO 360                                                          00610034
 350  DO 355 I=1,IPOPS                                                  00620034
 355  POPN(I)=0.                                                        00630034
      SUMA=0.                                                           00640034
 360  DO 362 I=1,IPOPS                                                  00650034
 362  POP(I)=0.                                                         00660034
      CALL SUMPOP(N2,M2,EXMAX2,POPA,SUMN,IPOPS,N,M ,EXMAX, POPN,        00670034
     1            SUMA       ,POP ,FLOST(NSTEP),IPRINT)                 00680036
      N1=N                                                              00690034
      M1=M                                                              00700034
      EXMAX1=EXMAX                                                      00710034
C                                                                       00720034
C                           TRIM NEW PARENT                             00730034
C                                                                       00740034
      CALL TOT(N1,M1,POP,SPECJ,AMX,AVGJ(NSTEP),TOTAL(NSTEP),SPECE,      00750034
     1         AVGE(NSTEP),EXMAX1)                                      00760034
      IF(TOTAL(NSTEP).GT.SUMLOW) GOTO 450                               00770034
      ISEQ(NSTEP)=0                                                     00780034
      RETURN                                                            00790034
 450  ALEVEL=0.                                                         00800034
      CUT1=CUT                                                          00810046
      PRCNT1=PRCNT                                                      00820046
      IF(NROW.EQ.1 .AND. NSTEP.EQ.1) CUT1=CUT*0.1                       00830046
      IF(NROW.EQ.1 .AND. NSTEP.EQ.1) PRCNT1=1.-(1.-PRCNT)*0.1           00840046
      CALL CUTOFF(N1,M1,POP,SPECJ,TOTAL,NSTEP,AMX,CUT1,PRCNT1,ICUT,     00850046
     1            TRIM(NSTEP),IPOPS,EXMAX1,ALEVEL,IPRINT)               00860036
      CALL TOT(N1,M1,POP,SPECJ,AMX,AVGJ(NSTEP),TOTAL(NSTEP),SPECE,      00870034
     1         AVGE(NSTEP),EXMAX1)                                      00880034
      IDUM=M1*N1                                                        00890043
      DO 455 I=1,IDUM                                                   00900043
 455  PENTRY(I)=POP(I)                                                  00910043
      CUTPOP=ALEVEL/2.                                                  00920034
      DUM=0.05*CUT*TOTAL(NSTEP)                                         00930048
      IF(ALEVEL.LE.0. .OR. CUTPOP.GT.DUM) CUTPOP=DUM                    00940048
      DUM=0.                                                            00950034
      MM1=MIN0(IEDIM,M1)                                                00960034
      DO 460 I=1,MM1                                                    00970034
 460  DUM=DMAX1(DUM,SPECE(I))                                           00980034
      IAMX=DUM*10000./TOTAL(NSTEP)  +0.5                                00990034
C                                                                       01000034
C                           SET LIMITS OF NEW DAUGHTERS                 01010034
C                                                       N               01020034
 505  K=2                                                               01030034
      N2=1                                                              01040034
      M2=1                                                              01050034
      EXMAX2=EXMAX1+Q(K)   -DEL(K)                                      01060034
      IF(EXMAX2.LE.0.)GOTO 518                                          01070034
      N22=MIN0(N1+5,JDIM)                                               01080034
      JZERO=JDIM*(K-1)+1                                                01090034
      DO 509 J1=1,N22                                                   01100034
      JYR=N22-J1+JZERO                                                  01110034
      IF(EJAY(JYR).LE.EXMAX2) GO TO 510                                 01120034
  509 CONTINUE                                                          01130034
  510 N2=JYR-JZERO+1                                                    01140034
C                                                                       01150034
  511 M22=EXMAX2+1.                                                     01160034
      M2=M1+KEN                                                         01170034
      M2=MIN0(M22,M2)                                                   01180034
      IF((M2*N2).LE.IPOPS) GOTO 518                                     01190034
      M2=IPOPS/N2                                                       01200034
 518  MN=M2*N2                                                          01210034
      DO  530  I=1,IPOPS                                                01220034
 530  POPN(I)=0.                                                        01230034
C                                                       P               01240034
C                                                                       01250034
 605  K=3                                                               01260034
      N3=1                                                              01270034
      M3=1                                                              01280034
      EXMAX3=EXMAX1+Q(K)   -DEL(K)                                      01290034
      IF(EXMAX3.LE.0.)GOTO 618                                          01300034
      N33=MIN0(N1+5,JDIM)                                               01310034
      JZERO=JDIM*(K-1)+1                                                01320034
      DO 609 J1=1,N33                                                   01330034
      JYR=N33-J1+JZERO                                                  01340034
      IF(EJAY(JYR).LE.EXMAX3) GO TO 610                                 01350034
  609 CONTINUE                                                          01360034
  610 N3=JYR-JZERO+1                                                    01370034
C                                                                       01380034
  611 M33=EXMAX3+1.                                                     01390034
      M3=M1+KEP                                                         01400034
      M3=MIN0(M33,M3)                                                   01410034
      IF((M3*N3).LE.IPOPS) GOTO 618                                     01420034
      M3=IPOPS/N3                                                       01430034
 618  MN=M3*N3                                                          01440034
      DO  630  I=1,IPOPS                                                01450034
 630  POPP(I)=0.                                                        01460034
C                                                                       01470034
C                                                        A              01480034
  705 K=4                                                               01490034
      N4=1                                                              01500034
      M4=1                                                              01510034
      EXMAX4=EXMAX1+Q(K)   -DEL(K)                                      01520034
      IF(EXMAX4.LE.0.) GOTO 718                                         01530034
      N44=MIN0(N1+5,JDIM)                                               01540034
      JZERO=JDIM*(K-1)+1                                                01550034
      DO 709 J1=1,N44                                                   01560034
      JYR=N44-J1+JZERO                                                  01570034
      IF(EJAY(JYR).LE.EXMAX4) GO TO 710                                 01580034
  709 CONTINUE                                                          01590034
  710 N4=JYR-JZERO+1                                                    01600034
C                                                                       01610034
  711 M44=EXMAX4+1.                                                     01620034
      M4=M1+KEP                                                         01630034
      M4=MIN0(M44,M4)                                                   01640034
      IF((M4*N4).LE.IPOPS) GOTO 718                                     01650034
      M4=IPOPS/N4                                                       01660034
 718  MN=M4*N4                                                          01670034
      DO  730  I=1,IPOPS                                                01680034
 730  POPA(I)=0.                                                        01690034
      IF(LPRINT.LT.2)                                                   01700034
     1WRITE(6,740) M1,N1,EXMAX1,M2,N2,EXMAX2,M3,N3,EXMAX3,M4,N4,EXMAX4  01710034
 740  FORMAT(' DIMENSIONS',4(2I4,F8.1))                                 01720034
      RETURN                                                            01730034
      END                                                               01740034
