      SUBROUTINE TOT(N,M,POP,SPECJ,AMX,AVGJ,TOTAL,SPECE,AVGE,EXMAX)     00020003
C
C  Bestimmt Summe und erstes Moment der E-J-Bevoelkerung
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),SPECJ(1),SPECE(1)                                00030003
      JDIM=100                                                          00040005
      IEDIM=400                                                         00050007
C----------------------------------------------------------             00060006
      SUM=0.                                                            00070000
      AMX=0.                                                            00080000
      INDEX=0                                                           00090000
      AVG=0.                                                            00100002
      AVGE=0.                                                           00110003
      DO 5 I=1,JDIM                                                     00120005
 5    SPECJ(I)=0.                                                       00130005
      DO 6 I=1,IEDIM                                                    00140005
 6    SPECE(I)=0.                                                       00150005
C----------------------------------------------------------------       00160006
      NN=MIN0(JDIM,N)                                                   00170006
      DO 50 JF1=1,NN                                                    00180006
      SUMJ=0.                                                           00190000
      FJ=JF1-1                                                          00200000
      DO 45 KE=1,M                                                      00210000
      JFKE=INDEX+KE                                                     00220000
      DUM=POP(JFKE)                                                     00230000
      SUMJ=SUMJ+DUM                                                     00240000
      AMX=DMAX1(DUM,AMX)                                                00250000
 45   CONTINUE                                                          00260000
      SPECJ(JF1)=SUMJ                                                   00270000
      SUM=SUM+SUMJ                                                      00280000
      AVG=AVG+SUMJ*FJ                                                   00290000
      INDEX=INDEX+M                                                     00300000
 50   CONTINUE                                                          00310000
C                                                                       00320003
      MM=MIN0(IEDIM,M)                                                  00330006
      DO 80 KE=1,MM                                                     00340006
      FKE=1-KE                                                          00350003
      SUME=0.                                                           00360003
      JFKE=KE-M                                                         00370003
      DO 70 JF1=1,N                                                     00380003
      JFKE=JFKE+M                                                       00390003
      SUME=SUME+POP(JFKE)                                               00400003
 70   CONTINUE                                                          00410003
      SPECE(MM+1-KE)=SUME                                               00420006
      AVGE=AVGE+FKE*SUME                                                00430004
 80   CONTINUE                                                          00440003
C                                                                       00450000
      TOTAL=SUM                                                         00460000
      IF(SUM.NE.0.)AVGJ=AVG/SUM                                         00470000
      IF(SUM.NE.0.)AVGE=AVGE/SUM +EXMAX                                 00480003
      RETURN                                                            00490000
      END                                                               00500000
