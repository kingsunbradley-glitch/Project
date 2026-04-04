      FUNCTION TABIP(TABLX,TABLY,NTAB,IDEG,ILOG,X,IPRINT)               00010006
C
C
C  Stand: 3.5.1994
C
C     INTERPOLATE FROM A TABLE                                          00020000
C     TABLX( ) X-VALUES OF TABLE,MUST BE IN ASCENDING ORDER             00030000
C     TABLY( ) Y-VALUES OF TABLE,LIMITED FOR THE MOMENT TO POSITIVE Y   00040005
C     NTAB     SIZE OF TABLE                                            00050000
C     X        ARGUMENT                                                 00060000
C     IDEG     0,1   LINEAR IP                                          00070000
C              2,>2  QUADRATIC IP                                       00080000
C     ILOG     >0    LOGARITHMIC (IN Y) IP                              00090000
C-------------------------------------------------------------------    00100000
      IMPLICIT REAL*8 (A-H,O-Z)  
      DIMENSION TABLX(NTAB),TABLY(NTAB)                                 00110008
C
C     I.Giese 2007
      SAVE
C
      IF(NTAB.LT.2) GOTO 999                                            00120000
      IDG=IDEG                                                          00130000
      IF(NTAB.LT.3) IDG=1                                               00140000
C-------------------------------------------------                      00150000
      IF(X.LT.TABLX(2)) GOTO 40                                         00160000
      IF(X.GT.TABLX(NTAB-1)) GOTO 50                                    00170000
      GOTO 60                                                           00180000
 40   I1=1                                                              00190000
      GOTO 120                                                          00200000
 50   I1=NTAB-2                                                         00210000
      IF(IDG.LE.1) I1=NTAB-1                                            00220000
      GOTO 120                                                          00230000
 60   DO 80 K=1,NTAB                                                    00240000
      IF(X.LT.TABLX(K)) GOTO 90                                         00250007
 80   CONTINUE                                                          00260000
      GOTO 999                                                          00270000
 90   I1=K-1                                                            00280007
      I1=MIN0(I1,NTAB-1)                                                00290000
      IF(IDG.LE.1) GOTO 120                                             00300000
      IF(X-TABLX(K) .LT. TABLX(K+1)-X) I1=K-1                           00310000
      I1=MAX0(1,I1)                                                     00320000
      I1=MIN0(I1,NTAB-2)                                                00330000
C-------------------------------------------------                      00340000
 120  I2=I1+1                                                           00350000
      I3=I1+2                                                           00360000
      X1=TABLX(I1)                                                      00370000
      Y1=TABLY(I1)                                                      00380000
      X2=TABLX(I2)                                                      00390000
      Y2=TABLY(I2)                                                      00400000
      IF(IDG.GE.2) X3=TABLX(I3)                                         00410000
      IF(IDG.GE.2) Y3=TABLY(I3)                                         00420000
      ILOG1=0                                                           00430004
      IF(ILOG.LE.0) GOTO 125                                            00440000
      IF(Y1.EQ.0.) GOTO 125                                             00450004
      Y1=DLOG(ABS(Y1))                                                  00460000
      Y2=DLOG(ABS(Y2))                                                  00470000
      IF(IDG.GE.2) Y3=DLOG(ABS(Y3))                                     00480000
      ILOG1=1                                                           00490004
 125  IF(IDG.LE.1) CALL IP1(X1,X2,X,Y1,Y2,Y)                            00500000
      IF(IDG.GE.2) CALL IP2(X1,X2,X3,X,Y1,Y2,Y3,Y)                      00510000
      IF(IPRINT.GT.0) WRITE(6,130) X1,X2,X3,X,Y1,Y2,Y3,Y                00520006
 130  FORMAT(' X1',4F8.2,4E11.3)                                        00530006
      TABIP=Y                                                           00540000
      IF(ILOG1.GT.0) TABIP=EXP(Y)                                       00550004
      RETURN                                                            00560000
 999  TABIP=0.                                                          00570000
      RETURN                                                            00580000
      END                                                               00590000
