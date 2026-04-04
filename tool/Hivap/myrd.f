      SUBROUTINE MYRD(ARRAY,IR  ,IPR,W)                                 00010010
C 
C  Unterprogramm zum Einlesen der Eingabedaten
C  Stand: 3.5.1994
C 
C     IRC1 NUMBER OF VARIABLES TO BE READ INTO 'ARRAY'                  00030008
C     IPR  LOG UNIT FOR PRINT OUT  (IF<6, NO PRINT OUT)                 00040005
C     W    LOG UNIT TO READ FROM                                        00050009
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ARRAY(36)                                               00060011
      INTEGER*4 W5,W                                                    00070009
      CHARACTER*1 AW(72),AST,COMM,BLANK                                 00080008
      DATA COMM,AST,BLANK/'C','*',' '/                                  00090008
C--------------------------------------                                 00100008
C
C     I.Giese 2007
      SAVE
C
      IRC1=IR                                                           00110010
      W5=W                                                              00120009
      W5=MAX0(5,W5)                                                     00130002
      IF(IRC1.EQ.0)RETURN                                               00140000
      IDUM=1                                                            00150000
      ITEST=0                                                           00160000
      IF(IRC1.LT.0)ITEST=1                                              00170000
      IRC1=IABS(IRC1)                                                   00180000
      I=1                                                               00190000
      IRC2=0                                                            00200000
 10   IRC=0                                                             00210000
C--------------------------------------                                 00220008
      CALL DECODN(ARRAY(I),AW,IRC,W5,72)                                00230003
C--------------------------------------                                 00240008
      IF(IRC.EQ.99) GOTO 99                                             00250003
      IF(IPR.GT.5) WRITE(IPR,15) AW                                     00260000
 15   FORMAT(1X,72A1)                                                   00270000
      IF(AW(1).EQ.AST .OR. (AW(1).EQ.COMM .AND. AW(2).EQ.BLANK))GOTO 10 00280008
      IF(ITEST.EQ.0) GOTO 18                                            00290000
      IF(IRC.GT.0) IDUM=ABS(ARRAY(1))+0.01                              00300000
      IF(IDUM.EQ.0) RETURN                                              00310000
 18   IRC2=IRC2+IRC                                                     00320000
      IF(IRC2.GE.IRC1) GOTO 20                                          00330000
      I=IRC2+1                                                          00340000
      GOTO 10                                                           00350000
 20   RETURN                                                            00360000
 99   IRC1=99                                                           00370003
      RETURN                                                            00380003
      END                                                               00390000
