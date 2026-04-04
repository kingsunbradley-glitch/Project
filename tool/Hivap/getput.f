      SUBROUTINE GETPUT(POP,IPOPS,IUNIT,IREAD,M,N,EXMAX,MZ,MA,ACN,ZCN)  00020003
C  
C  Speichert oder liest die Puffer der E-J - Bevoelkerung
C  (logische Einheiten 10-13)
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION POP(4096)                                               00030003
C
C     I.Giese 2007
      SAVE
C
      L=IPOPS/256                                                       00040002
      I1=1                                                              00050001
      IF(IREAD.EQ.1) GOTO 80                                            00060001
      WRITE(IUNIT,101) M,N,EXMAX,MZ,MA,ACN,ZCN,IUNIT                    00070003
      IF(M.EQ.0) RETURN                                                 00080003
      DO 50 NREC=1,L                                                    00090002
      I2=I1+255                                                         00100001
      DO 79 I=I1,I2
      WRITE(IUNIT,102) POP(I)                                           00110001
 79   CONTINUE
      I1=I1+256                                                         00120001
 50   CONTINUE                                                          00130001
      RETURN                                                            00140001
 80   READ (IUNIT,101) M,N,EXMAX,MZ,MA,ACN,ZCN,IDUM                     00150004
      IF(M.EQ.0) RETURN                                                 00160003
      DO 90 NREC=1,L                                                    00170002
      I2=I1+255                                                         00180001
      DO 81 I=I1,I2
      READ(IUNIT,102) POP(I)                                            00190001
 81   CONTINUE
      I1=I1+256                                                         00200001
 90   CONTINUE                                                          00210001
      RETURN                                                            00220001
 101  FORMAT(2I4,F8.2,2I4,2F7.3,I4)                                     00230001
 102  FORMAT(E10.5)
      END                                                               00240001
