      SUBROUTINE IP2(X0,X1,X2,X,F0,F1,F2,F)                             00010002
C
C  Berechnet F durch quadratische Interpolation                         00020002
C  X0,X1,X2 in aufsteigender Folge
C  Stand: 3.5.1994
C                                                                       00040002
      IMPLICIT REAL*8 (A-H,O-Z)
C     SPECIAL CASES                                                     00050002
      IF(X0.EQ.X1 .AND. X1.EQ.X2) GOTO 20                               00060002
      IF(X0.EQ.X1 .AND. X1.NE.X2) GOTO 25                               00070002
      IF(X0.EQ.X2 .AND. X1.NE.X2) GOTO 25                               00080002
      IF(X1.EQ.X2 .AND. X1.NE.X0) GOTO 35                               00090002
      GOTO 50                                                           00100002
 20   F=(F0+F1+F2)/3.                                                   00110002
      RETURN                                                            00120002
 25   CALL IP1(X1,X2,X,F1,F2,F)                                         00130002
      RETURN                                                            00140002
 35   CALL IP1(X0,X1,X,F0,F1,F)                                         00150002
      RETURN                                                            00160002
C                                                                       00170002
C     NORMAL CASE  FROM APPL.NUM.METH.(1969)11                          00180002
C                                                                       00190002
 50   F10=(F1-F0)/(X1-X0)                                               00200002
      F21=(F2-F1)/(X2-X1)                                               00210002
      F210=(F21-F10)/(X2-X0)                                            00220002
      F=F0+(X-X0)*F10+(X-X0)*(X-X1)*F210                                00230002
      RETURN                                                            00240002
      END                                                               00250002
