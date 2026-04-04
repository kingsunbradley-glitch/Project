      SUBROUTINE OUT1(SPEC,N,ZERO,D,LOGUN)                              00020005
C
C  Druck-Routine fuer eindimensionales Array
C  Stand: 3.5.1994
C
C     PRINTS 1-DIMENSIONAL NORMALISED ARRAY 'SPEC' OF DIMENSION 'N'     00030001
C     ALSO OUTPUTS SUM AND FIRST MOMENT                                 00040001
C     'ZERO' IS THE VALUE OF THE ARGUMENT FOR SPEC(1)                   00050001
C     ON OUTPUT 'ZERO' IS THE FIRST MOMENT OF 'SPEC'                    00060007
C     THE ARGUMENT INCREASES IN STEPS 'D'                               00070001
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SPEC(N),IP(10)                                          00080004
      SUM=0.                                                            00090000
      FM=0.                                                             00100000
      DO 100 I=1,N                                                      00110000
      SUM=SUM+SPEC(I)                                                   00120000
      FI=FLOAT(I-1)*D+ZERO                                              00130001
      FM=FM+FI*SPEC(I)                                                  00140000
 100  CONTINUE                                                          00150000
      IF(SUM.LE.0.) GOTO 150                                            00160008
      FM=FM/SUM                                                         00170000
      IF(LOGUN.EQ.0) GOTO 140                                           00180007
      WRITE(LOGUN,105) FM,SUM                                           00190005
 105  FORMAT(             4X,'AVERAGE=',F8.3,4X,'SUM=',E12.3)           00200004
      FNORM=10000./SUM                                                  00210004
      DO 120 I10=1,N,10                                                 00220006
      E=ZERO+(I10-1)*D                                                  00230001
      I2=MIN0(N,I10+9)                                                  00240006
      ISUM=0                                                            00250004
      I3=0                                                              00260006
      DO 110 I=I10,I2                                                   00270006
      I3=I3+1                                                           00280006
      IP(I3)=SPEC(I)*FNORM                                              00290006
 110  ISUM=ISUM+IP(I3)                                                  00300006
      IF(ISUM.LE.0) GOTO 120                                            00310004
      WRITE (LOGUN,130) E,(IP (I),I=1,I3 )                              00320006
 120  CONTINUE                                                          00330000
 130  FORMAT(F8.2,10I7)                                                 00340004
 140  ZERO=FM                                                           00350007
      RETURN                                                            00360000
 150  ZERO=0.                                                           00370008
      RETURN                                                            00380009
      END                                                               00390000
