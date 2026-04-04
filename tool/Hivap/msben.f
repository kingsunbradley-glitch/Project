      SUBROUTINE  MSBEN(ZCN,ACN,NA,NZ,CMASS,SHLL)                       00020006
C
C  Berechnet Separationsenergien unter Benutzung der Lysekyl -
C  Parametrisierung (1967) von Myers-Swiatecki
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BIND/BE(25,25,5),EXC(27,27),MC,MP,NUMB,MASSES,DELT         00030010
      DIMENSION EM(10),XK(10),Y(2),F(2),EMP(10),XQ(20)                  00040000
C     NEEDED EXTERNAL ACN,ZCN,BE,MC,MP,NZ,NA                            00050000
C     MP=1 NO PAIRING FLAG                                              00060000
C     MC=1 NO SHELL CORRECTION FLAG                                     00070000
C     DELIVERS SEPARATION ENERGIES OF N,P,A,D,T IN A SQUARE NA*NZ OF THE00080000
C     ISOTOPE TABLE;(ACN,ZCN) IS UPPER RIGHT CORNER;USES MYERS-SWIATECKI00090000
C     1967 MASS FORMULA;BE ARRAY FOR SEPARATION ENERGIES,EXC=ARRAY      00100000
C     CONTAINING MASS EXCESSES                                          00110000
C-----------------------------------------------------------------------00120000
      DATA EM/0.0,2.,8.,14.,28.,50.,82.,126.,184.,258./,                00130008
     1     CAY1,CAY2,CAY3,CAY4,CAY5,CAY6/1.15303,0.0,200.,11.,8.07144,  00140008
     2     7.28899/,GAMMA/1.7826/,A1,A2,A3/15.4941,17.9439,0.7053/,     00150008
     3     D,C,SMALC/0.444,5.8,0.325/                                   00160008
C-----------------------------------------------------------------------00170008
C
C     I.Giese 2007
      SAVE
C
      IF(DELT.LE.0.) DELT=CAY4                                          00180011
C-----------------------------------------------------------------------00190006
      IFLAG=0                                                           00200006
      IF( (NA+NZ).LE.2 ) IFLAG=1                                        00210006
      DO 3 I=1,10                                                       00220000
      EMP(I)=EM(I)**(5.0/3.0)                                           00230000
 3    CONTINUE                                                          00240000
      DO 4 I=1,9                                                        00250000
      XK(I)=0.6*(EMP(I+1)-EMP(I))/(EM(I+1)-EM(I))                       00260000
 4    CONTINUE                                                          00270000
      RZ=.863987/A3                                                     00280000
      L=0                                                               00290000
      Z=1.0                                                             00300000
      KZ=ZCN+0.01                                                       00310012
      KA=ACN+0.01                                                       00320012
      NNZ=1                                                             00330008
      NNA=1                                                             00340008
      IF(IFLAG.EQ.1) GOTO 8                                             00350008
      NNZ=NZ+2                                                          00360008
      NNA=NA+2                                                          00370000
C                                                                       00380000
C********************************** LOOPS NZ,NA *********************** 00390000
 8    DO 32 JZ=1,NNZ                                                    00400000
      DO 32 JA=1,NNA                                                    00410000
 13   IA=KA+2-JA-JZ                                                     00420000
      IZ=KZ+1-JZ                                                        00430000
      N=IA-IZ                                                           00440000
 14   Z=IZ                                                              00450000
      UN=N                                                              00460000
      A=IA                                                              00470000
 15   A3RT=A**(1.0/3.0)                                                 00480000
      A2RT=SQRT(A)                                                      00490000
      A3RT2=A3RT**2.0                                                   00500000
      ZSQ=Z**2.0                                                        00510000
      SYM=((UN-Z)/A)**2                                                 00520000
      ACOR=1.0-GAMMA*SYM                                                00530000
      PARMAS=CAY5*UN+CAY6*Z                                             00540000
      VOLNUC=-1.0*A1*ACOR*A                                             00550000
      SUFNUC=A2*ACOR*A3RT2                                              00560000
      COULMB=A3*ZSQ/A3RT                                                00570000
      FUZSUR=-1.0*CAY1*ZSQ/A                                            00580000
      ODDEV=0.                                                          00590008
      IF(MP.EQ.0)                                                       00600008
     1ODDEV=-1.0*(1.0+2.0*(N/2)-UN+2.*(IZ/2)-Z)/SQRT(A)*DELT            00610010
      WTERM=-1.*CAY2*A3RT2*EXP(-1.*CAY3*SYM)                            00620000
      WOTNUC=PARMAS+COULMB+FUZSUR+ODDEV+WTERM                           00630000
      SMASS=WOTNUC+VOLNUC+SUFNUC                                        00640000
      SHLL=0.                                                           00650008
      XQ(JZ)=SMASS                                                      00660000
C-----------------------------------------------------------------------00670006
C                                         MS   SHELL CORRECTION         00680008
      IF(MC.NE.0) GOTO 31                                               00690008
 16   CONTINUE                                                          00700000
      C2=(SUFNUC+WTERM)/(A**(2.0/3.0))                                  00710000
      X=COULMB/(2.0*(SUFNUC+WTERM))                                     00720000
 17   BARR=0.0                                                          00730000
 18   Y(1)=UN                                                           00740000
      Y(2)=Z                                                            00750000
      DO 22 J=1,2                                                       00760000
      DO 19 I=1,9                                                       00770000
      IF (EM(I+1).GT.(Y(J)-0.1)) GOTO 21                                00780011
 19   CONTINUE                                                          00790000
 20   STOP                                                              00800000
 21   F(J)=XK(I)*(Y(J)-EM(I))-.6*(Y(J)**(5./3.)-EMP(I))                 00810000
 22   CONTINUE                                                          00820000
      S=(2.0/A)**(2.0/3.0)*(F(1)+F(2))-SMALC*A**(1./3.)                 00830000
      EE=2.*C2*D**2*(1.0-X)                                             00840000
      FF=.42591771*C2*D**3*(1.+2.*X)/A3RT                               00850000
      SSHELL=C*S                                                        00860000
      V=SSHELL/EE                                                       00870000
      EPS=1.5*FF/EE                                                     00880000
      IF(EE*(1.-3.*V).LE.0.0) GO TO 23                                  00890000
      QCALC=0.0                                                         00900000
      THETA=0.0                                                         00910000
      SHLL=SSHELL                                                       00920000
      GO TO 31                                                          00930000
 23   TO=1.0                                                            00940000
 24   DO 25 IPQ=1,10                                                    00950000
      EXPDUM=0.                                                         00960005
      DUM=TO**2                                                         00970005
      IF(DUM.LT.25.) EXPDUM=EXP(-DUM)                                   00980011
      T=TO-(1.-EPS*TO-V*(3.-2.*TO**2)*EXPDUM     )/(-EPS+V*(10.*TO-4.   00990005
     1 *TO**3)*EXPDUM     )                                             01000005
      IF (T.LE.0.0) GO TO 27                                            01010000
      IF (ABS(T-TO) .LT.0.0001) GO TO 26                                01020000
      TO=T                                                              01030000
 25   CONTINUE                                                          01040000
      GO TO 29                                                          01050000
 26   EXPDUM=0.                                                         01060005
      DUM=T **2                                                         01070005
      IF(DUM.LT.25.) EXPDUM=EXP(-DUM)                                   01080011
      IF (2.*EE*(1.-2.*EPS*T-V*(3.-12.*T**2+4.*T**4)*EXPDUM    )        01090005
     1 .GT.0.0) GO TO 30                                                01100000
 27   DO 28 I=1,20                                                      01110000
      TO=FLOAT(I)/10.                                                   01120000
      GL=EE*(1.-EPS*TO-V*(3.-2.*TO**2)*EXP(-TO**2))                     01130000
      IF (GL.GE.0.0) GO TO 24                                           01140000
 28   CONTINUE                                                          01150000
 29   CONTINUE                                                          01160000
      GO TO 32                                                          01170000
 30   THETA=T                                                           01180000
      ALPHA0=D*SQRT(5.)/A**(1./3.)                                      01190000
      ALPHA=ALPHA0*THETA                                                01200000
      SIGMA=ALPHA*(1.+ALPHA/14.)                                        01210000
      QCALC=.004*Z*(RZ*A3RT)**2*(EXP(2.*SIGMA)-EXP(-SIGMA))             01220000
      SHLL=EE*T**2-FF*T**3+SSHELL*(1.-2.*T**2)*EXP(-T**2)               01230000
 31   CMASS=SMASS+SHLL                                                  01240000
      IF(IFLAG.EQ.0)      EXC(JA,JZ)=CMASS                              01250006
      XQ(JZ)=CMASS                                                      01260000
 32   CONTINUE                                                          01270000
C                                                                       01280000
C******************************* END LOOPS NA,NZ ********************** 01290000
 36   IF(IFLAG.EQ.1)RETURN                                              01300006
      DO 37 JZ=1,NZ                                                     01310006
      DO 37 JA=1,NA                                                     01320000
      BE(JZ,JA,1)=8.07+EXC(JA+1,JZ)-EXC(JA,JZ)                          01330000
      BE(JZ,JA,2)=7.29+EXC(JA,JZ+1)-EXC(JA,JZ)                          01340000
      BE(JZ,JA,4)=13.14+EXC(JA+1,JZ+1)-EXC(JA,JZ)                       01350000
      BE(JZ,JA,5)=14.95+EXC(JA+2,JZ+1)-EXC(JA,JZ)                       01360000
 37   BE(JZ,JA,3)=2.42+EXC(JA+2,JZ+2)-EXC(JA,JZ)                        01370000
 38   RETURN                                                            01380000
      END                                                               01390000
