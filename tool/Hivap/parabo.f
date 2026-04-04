      SUBROUTINE PARABO(X,Y,A,X0,Y0)                                    00020000
C
C  Sucht Ort und Wert des Maximums einer Parabel
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(3),Y(3),A(3)                                          00030000
      X0=0.                                                             00040000
      Y0=0.                                                             00050000
      DUM=X(2)-X(1)                                                     00060000
      IF(DUM.EQ.0.)RETURN                                               00070000
      F21=(Y(2)-Y(1))/DUM                                               00080000
      DUM=X(3)-X(2)                                                     00090000
      IF(DUM.EQ.0.)RETURN                                               00100000
      F32=(Y(3)-Y(2))/DUM                                               00110000
      DUM=X(3)-X(1)                                                     00120000
      IF(DUM.EQ.0.)RETURN                                               00130000
      F321=(F32-F21)/DUM                                                00140001
      A(2)=F21-F321*(X(1)+X(2))                                         00150000
      A(3)=F321                                                         00160000
      A(1)=Y(1)-X(1)*(F21-X(2)*F321)                                    00170001
      IF(A(3).EQ.0.)RETURN                                              00180002
      X0=-0.5*A(2)/A(3)                                                 00190002
      Y0=-0.25*A(2)*A(2)/A(3) + A(1)                                    00200000
      RETURN                                                            00210000
      END                                                               00220000
