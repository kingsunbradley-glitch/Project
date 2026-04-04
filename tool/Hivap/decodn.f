      SUBROUTINE DECODN (ARRAY,A,IRC,LOGUN,NN)                          00010005
C
C  Unterprogramm zum Einlesen der Eingabedaten
C  Stand: 3.5.1994
C 
C     'A'      IS READ IN A-FORMAT AND/OR DECODED                       00040000
C     'IRC'    OUTPUT:NUMBER OF NUMBERS FOUND                           00050000
C              OUTPUT:99 END OF FILE FOUND                              00060000
C              OUTPUT:-1 DECIMAL EXPONENT EXCEEDS 99                    00070000
C     'LOGUN'  LOG.UNIT FROM WHICH 'A' IS READ                          00080000
C              0  NO READING                                            00090000
C              1-4  =5                                                  00100000
C     'NN'     'A' WILL BE READ AND/OR DECODED FROM A(1) TO A(N)        00110005
C              DEFAULT FOR ZERO IS 72                                   00120000
C---------------------------------------------------------------------  00130000
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ARRAY(1)                                                00140000
      CHARACTER B,A(1),SY(15),MINUS,ESIGN,PERIOD,BLANK,COMMA,EQUAL,DOPP,00150004
     1          SEMI                                                    00160001
      CHARACTER SLASH                                                   00170004
      EQUIVALENCE (SY(11),MINUS),(SY(14),ESIGN),(SY(12),PERIOD),(SY(15),00180000
     1 COMMA)                                                           00190001
      DATA SY,BLANK/'0','1','2','3','4','5','6','7','8','9','-','.','+',00200000
     1 'E',',',' '/,EQUAL,DOPP  /'=',':'/,SLASH/'/'/,SEMI/';'/          00210003
C---------------------------------------------------------------------  00220000
C
C     I.Giese 2007
      SAVE
C
C                                                         INITIALIZE    00230000
      N=NN                                                              00240005
      IF(NN.LE.0) N=72                                                  00250005
      N1=N+1                                                            00260000
      IFLAG=0                                                           00270000
      ISEP=1                                                            00280000
      DEC=0.                                                            00290000
      NUMB=0                                                            00300000
      IE=0                                                              00310000
      IEST=0                                                            00320000
      S=1.                                                              00330000
      IS=1                                                              00340000
      IW=1                                                              00350000
      IBLANK=1                                                          00360000
      IW1=1                                                             00370000
      L=1                                                               00380000
      F=1.                                                              00390000
      IF(LOGUN.LT.1) GOTO 17                                            00400000
      LOGUN1=MAX0(LOGUN,5)                                              00410000
C---------------------------------------------------------------------- 00420000
C                                           READ UP TO N  SYMBOLS       00430000
      READ(LOGUN1,15,END=99) (A(I),I=1,N)                               00440000
 15   FORMAT(72A1,50A1)                                                 00450000
C---------------------------------------------------------------------- 00460000
C                                                SYMBOL LOOP            00470000
 17   IRC=0                                                             00480000
      DO 80 I=1,N1                                                      00490000
      IF(I.EQ.N1) GOTO 25                                               00500000
      B=A(I)                                                            00510000
C----------------------------- IS IT ONE OF THE SYMBOLS 'SY' (0 TO +)   00520000
      DO 20 K=1,13                                                      00530000
      IF(B.EQ.SY(K) .AND.ISEP.EQ.1) GOTO 30                             00540000
 20   CONTINUE                                                          00550000
C----------------------------- NO IT'S A SEPARATION                     00560000
 25   IF(B.EQ.ESIGN .AND. IBLANK.EQ.0) GOTO 60                          00570000
      ISEP=0                                                            00580000
      IF(B.EQ.COMMA .OR. B.EQ.BLANK .OR. B.EQ.EQUAL .OR. B.EQ.DOPP      00590001
     1   .OR. B.EQ.SLASH .OR. B.EQ.SEMI .OR. I.EQ.N1) ISEP=1            00600001
      IBLANK=IBLANK+1                                                   00610000
      IEST=0                                                            00620000
      IF(IBLANK.NE.1) GOTO 50                                           00630000
      IF(IFLAG.EQ.0) GOTO 28                                            00640000
      IF(ISEP.NE.1) GOTO 55                                             00650001
      IW=IW+1                                                           00660000
      IFLAG=0                                                           00670000
      GO TO 50                                                          00680000
 28   IS=1                                                              00690000
      S=1.                                                              00700000
      L=1                                                               00710000
      GOTO 50                                                           00720000
C----------------------------- IS IT A NUMBER OR IS IT + - E .          00730000
 30   IBLANK=0                                                          00740000
      IF(IEST.EQ.0) GOTO 32                                             00750000
      IEST=0                                                            00760000
      L=3                                                               00770000
 32   IF(K.GT.10) GOTO 70                                               00780000
C----------------------------- IT'S A NUMBER                            00790000
      IFLAG=1                                                           00800000
      GOTO (35,40,45),L                                                 00810000
 35   NUMB=NUMB*10+K-1                                                  00820000
      GOTO 80                                                           00830000
 40   F=F*0.1                                                           00840000
      DEC=DEC+(K-1)*F                                                   00850000
      GOTO 80                                                           00860000
 45   IE=IE*10 + K-1                                                    00870000
      IF(IE.LT.99) GOTO 80                                              00880000
      IRC=-1                                                            00890000
      RETURN                                                            00900000
C---------------------------- IS IT END OF NUMBER-WORD                  00910000
 50   IF(IW.EQ.IW1) GOTO 80                                             00920000
C---------------------------- YES IT IS                                 00930000
      R=NUMB                                                            00940000
      D=10.**FLOAT(IE)                                                  00950000
      IF(IS.LT.0)ARRAY(IW-1)=S*(R+DEC)/D                                00960000
      IF(IS.GT.0)ARRAY(IW-1)=S*(R+DEC)*D                                00970000
 55   IW1=IW                                                            00980001
      DEC=0.                                                            00990000
      F=1.                                                              01000000
      NUMB=0                                                            01010000
      IS=1                                                              01020000
      S=1.                                                              01030000
      IE=0                                                              01040000
      L=1                                                               01050000
      GOTO 80                                                           01060000
C---------------------------- E-SIGN FLAG                               01070000
 60   IEST=1                                                            01080000
      GOTO 80                                                           01090000
C---------------------------- IT'S . OR E OR MINUS                      01100000
 70   IF(B.EQ.PERIOD) L=2                                               01110000
      IF(B.EQ.MINUS.AND.L.EQ.1) S=-1.                                   01120000
      IF(B.EQ.MINUS .AND. L.EQ.3) IS=-1                                 01130000
C----------------------------                                           01140000
 80   CONTINUE                                                          01150000
C------------------------------------------------------------------     01160000
C                                         END OF SYMBOL LOOP            01170000
      IRC=IW-1                                                          01180000
      RETURN                                                            01190000
 99   IRC=99                                                            01200000
      RETURN                                                            01210000
      END                                                               01220000
