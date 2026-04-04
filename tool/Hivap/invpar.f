      SUBROUTINE INVPAR(EQ,F,RZERO,D,AP,ZP,AT,ZT,SIGML,FLGRAZ)          0000000 
C                                                       
C  Vereinfachte Version von PARAP, wird zusammen mit FUSION benutzt
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIGML(1)                                                0000000 
C                                                                       0000000 
C     MAYBE UNINITIALIZED: TL0,TL1
C
C     I.Giese 2007
      SAVE
C
C     CALCULATION OF NUCLEAR CROSS SECTIONS                             71910837
C     WITH PARABOLIC BARRIER - SEE THOMAS, PHYS REV 116,703 (1959) -    71910838
      EN=EQ                                                             71910841
      DO1 I=1,200                                                       71910842
 1    SIGML(I)=0.                                                       71910845
      LMAX=200                                                          71910846
 2    WRITE(6,1001)                                                     71910848
 1001 FORMAT(1H0,'INVERTED PARABOLA'/)                                  0000000 
      LMAX=LMAX+1                                                       71910850
      COUL=(1.4393)*ZT*ZP                                               71910851
      AR=AT**0.333333+AP**0.333333                                      71910854
      RA=RZERO*AR                                                       0000000 
C                                                                       0000000 
      PRINT 1007,F,D,RZERO,RA                                           0000000 
 1007 FORMAT(1H0,'F,D,RZERO,RA',4F10.3/)                                0000000 
      PRINT 1003                                                        0000000 
 1003 FORMAT(5X,'L',6X,'TL',13X,'RB',10X,'VB', 10X,'WL',10X,'SIGL')     0000000 
      U=AT*AP/(AT+AP)                                                   71910863
      H2=41.814                                                         71910864
      ACN=AP+AT                                                         71910865
      CR2=D*COUL/F                                                      71910867
      CR3P=H2*D/(F*U)                                                   71910868
 6    E=EN*AT/(AT+AP)                                                   71910870
      EL=EN                                                             71910871
      PWAVE=SQRT(2.*931.16*U*E)                                         71910902
      WAVEL=19.732/PWAVE                                                71910903
      ARWAVE=3141.59*(WAVEL**2)                                         71910904
      TOTSIG=0.                                                         0000000 
      TL2=0.                                                            0000000 
      IREP=0                                                            0000000 
      RB=0.                                                             0000000 
C                                                                       0000000 
      DO 29 K=1,LMAX                                                    71910909
      L=K-1                                                             71910910
      R=RB                                                              0000000 
      IF(L.EQ.0)R=RA                                                    0000000 
      RB=0.                                                             0000000 
      AL=FLOAT(L)                                                       71910912
      CR3=CR3P*AL*(AL+1)                                                71910913
 9    DUM=CR2/R**2+CR3/R**3                                             0000000 
      IF(DUM)91,91,90                                                   0000000 
 91   IF(RB)900,900,11                                                  0000000 
 90   RB=RA-D*DLOG(DUM)                                                 0000000 
      DELTR=ABS(RB-R)                                                   71910915
      IF(DELTR-0.01)11,10,10                                            0000000 
 10   R=RB                                                              71910917
      GO TO 9                                                           71910918
 11   CONTINUE                                                          71910919
      R0=RB/AR                                                          71910920
      VLB=COUL/RB+(H2/2./U)*(AL*(AL+1)/RB**2)-F*EXP(-(RB-RA)/D)         71910921
      D2VLB=2.*COUL/RB**3+(H2/U)*3.*AL*(AL+1)/RB**4                     71910922
     1 -(F/D**2)*EXP(-(RB-RA)/D)                                        71910923
      IF(D2VLB) 12,13,14                                                71910924
 12   SD2VLB=-1.0                                                       71910925
      GO TO 15                                                          71910926
 13   SD2VLB=0.0                                                        71910927
      GO TO 15                                                          71910928
 14   SD2VLB=1.0                                                        71910929
 15   CONTINUE                                                          71910930
      WL2=ABS((H2/U)*D2VLB)                                             71910931
      WL=SQRT(WL2)                                                      71910932
      DUM=2.*3.14159*(VLB-E)/WL                                         0000000 
      IF(DUM.LT.-20.) GOTO 152                                          0000000 
      TL=1./(1.+EXP(DUM))                                               0000000 
      GOTO 155                                                          0000000 
 152  TL=1.                                                             0000000 
 155  CONTINUE                                                          0000000 
      RSIGL=TL*(2.*AL+1.)                                               71910934
      SIGL=ARWAVE*RSIGL                                                 71910935
      SIGML(K)=SIGL                                                     71910936
      TOTSIG=TOTSIG+SIGL                                                0000000 
      IF(TL.LT.0.5) IREP=IREP+1                                         0000000 
      IF(IREP.GE.2) GOTO 20                                             0000000 
      FLGRAZ=L-1                                                        0000000 
      TL1=TL2                                                           0000000 
      TL2=TL                                                            0000000 
      IF(L.EQ.0) TL0=TL                                                 0000000 
 20   CONTINUE                                                          0000000 
      IF( (TL.LT.0.99).OR.(L.EQ.1) )                                    0000000 
     1PRINT 1011,L,TL,RB,VLB,WL,SIGL                                    0000000 
 1011 FORMAT(I6,F12.6,4F12.3)                                           0000000 
      IF(TL-.0001)30,29,29                                              71900974
 29   CONTINUE                                                          71900975
 30   CONTINUE                                                          0000000 
      LMAX=K-1                                                          0000000 
      IF(TL0.LT.0.5) GOTO 32                                            0000000 
      FLGRAZ=FLGRAZ+(TL1-0.5)/(TL1-TL2)                                 0000000 
 32   IF(FLGRAZ.LT.0.)FLGRAZ=0.                                         0000000 
      FLGR2=FLGRAZ*FLGRAZ                                               0000000 
      WRITE(6,1012) FLGRAZ,LMAX,TOTSIG,FLGR2,E                          0000000 
 1012 FORMAT(//1H0,'LGRAZE,LMAX,TOTSIG,LGRAZE2,ECM',F8.1,I8,4X,         0000000 
     1 F8.3,F8.0,F10.1//)                                               0000000 
C                                                                       0000000 
 900  PRINT 1020                                                        0000000 
 1020 FORMAT(1H0,'NEG.ARG. FOR LOG'/)                                   0000000 
      RETURN                                                            0000000 
      END                                                               0000000 
