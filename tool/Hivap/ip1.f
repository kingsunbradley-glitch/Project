      SUBROUTINE IP1(X0,X1,X,Y0,Y1,Y)                                   00010001
C
C
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X0.EQ.X1)GOTO 10                                               00020001
      F10=(Y1-Y0)/(X1-X0)                                               00030001
      Y=Y0+F10*(X-X0)                                                   00040001
      RETURN                                                            00050001
 10   Y=(Y0+Y1)/2.                                                      00060001
      RETURN                                                            00070001
      END                                                               00080001
