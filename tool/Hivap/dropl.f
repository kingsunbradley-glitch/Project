      SUBROUTINE DROPL(Z,A,ICODE,DROP,DEL,EPS,LDM,IPRNT) 
C
C
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BEES/BS,BC,BK,BR,BV,BW                                     00040000
C     MAYBE UNINITIALIZED: A1,A2,A3,ZJ,C1,C2,C3,C4P,C5,T3,T4,T5,T6,
C                          T7,T8,T9,T10
C
C     I.Giese 2007
      SAVE
C
      IF(LDM.EQ.1) GOTO 15                                              00050000
C----------------------------------------                               00060000
             A1=               15.96                                    00070000
             A2=               20.69                                    00080000
             A3=               0.                                       00090000
             ZJ=               36.8                                     00100000
             Q=                17.0                                     00110000
             RZ=               1.18                                     00120000
             ZK=               240.                                     00130000
             ZL=               100.                                     00140000
             ZM=               0.                                       00150000
C-----------------------------------------------------                  00160000
          RNULL = RZ                                                    00170000
         R3 = 1.E0/3.E0                                                 00180000
         C1 = 0.6E0*1.4399784E0/RZ                                      00190000
         C2 = (C1**2.E0/168.E0) * (0.5E0/ZJ + 9.E0/ZK)                  00200000
         C3 = 2.5E0 * C1 * (0.99E0/RZ)**2                               00210000
         C4 = 1.25E0 * (1.5E0 / 3.14159265359E0)**(2.E0/3.E0)*C1        00220000
         C4P = C4/1.25992105E0                                          00230000
         C5 = (C1**2/Q)/64.E0                                           00240000
         T3 = 0.1875E0 * C1 / Q                                         00250000
         T4 = 2.25E0 * ZJ/Q                                             00260000
         T5 = 2.E0 * A2 / ZK                                            00270000
         T6 = ZL / ZK                                                   00280000
         T7 = C1 / ZK                                                   00290000
         T8 = ZJ * T4                                                   00300000
         T9 = 0.5E0 * ZK                                                00310000
         T10 = 0.5E0 * ZM                                               00320000
C------------------------------------                                   00330000
 15      IF(ICODE.EQ.2) GOTO 20                                         00340000
         BS=1.                                                          00350000
         BC=1.                                                          00360000
         BK=1.                                                          00370000
         BR=1.                                                          00380000
         BV=1.                                                          00390000
         BW=1.                                                          00400000
C-------------------------------------                                  00410000
 20      X = A**.333333                                                 00420004
         X2 = X*X                                                       00430000
         AS = 1.E0 - 2.E0 * Z/A                                         00440000
C                                                                       00450000
         IF(LDM.EQ.1) GOTO 40                                           00460000
         DEL = (AS + T3 * Z / X2 * BV)/(1.E0 + T4 / X * BS )            00470000
         D2 = DEL*DEL                                                   00480000
         EPS = (- T5 / X * BS + T6 * D2 + T7 * Z * Z / X2**2 * BC)      00490000
         E2 = EPS*EPS                                                   00500000
C                                                                       00510000
         V1=-A1*A                                                       00520000
         V2=ZJ*D2*A                                                     00530000
         V3=-  E2*A *T9                                                 00540000
         V4=T10*D2*D2*A                                                 00550000
         S1=A2*X2*BS                                                    00560000
         S2=T8*D2*X2*BS                                                 00570000
         CURV=A3*X*BK                                                   00580000
         COUL1=Z*Z*C1*BC/X                                              00590000
         COUL2=- Z*C2*X*BR *Z                                           00600000
         COUL3=-Z*Z*C3/A                                                00610000
         COUL4=-C4P*Z                                                   00620000
         COUL5=- Z*C5*BW*Z                                              00630000
C--------------------------------------------------------               00640000
      IF(LDM.EQ.0) GOTO 50                                              00650000
 40   A1=15.4941                                                        00660000
      RKAPPA=1.7826                                                     00670001
C     RKAPPA=2.815                                                      00680001
      ZJJ=A1*RKAPPA                                                     00690000
      A2=17.9439                                                        00700001
C     A2=18.7081                                                        00710001
      C1=0.7053                                                         00720000
C     C1=0.6765                                                         00730000
      T88=-A2*RKAPPA                                                    00740000
      AS2=AS*AS                                                         00750000
      V1=   -A1*A                                                       00760000
      V2=   ZJJ*AS2*A                                                   00770000
      V3=0.                                                             00780000
      V4=0.                                                             00790000
      S1= A2*X2*BS                                                      00800000
      S2=T88*AS2*X2*BS                                                  00810000
      CURV=0.                                                           00820000
      COUL1= Z*Z*C1*BC/X                                                00830000
      COUL2=0.                                                          00840000
      COUL3=-Z*Z*1.15303/A                                              00850000
      COUL4=0.                                                          00860000
      COUL5=0.                                                          00870000
C----------------------------------------------------------             00880000
  50  DROP=COUL1+COUL2+COUL3+COUL4+COUL5+V1+V2+V3+V4+S1+S2+CURV         00890000
      IF(IPRNT.NE.1)RETURN                                              00900000
      WRITE(6,101) Z,A,LDM,DROP,EPS,DEL                                 00910000
 101  FORMAT(/' Z,A',2F6.0,4X,'LDM',I4,4X,'DROP',F10.2/' EPS,DEL',2E12.400920000
     1)                                                                 00930000
      WRITE(6,102) V1                                                   00940000
 102  FORMAT(' BULK VOLUME', 6X,F10.2)                                  00950000
      WRITE(6,103) V2                                                   00960000
 103  FORMAT(' ASYM VOLUME', 6X,F10.2)                                  00970000
      WRITE(6,104) V3                                                   00980000
 104  FORMAT(' COMP VOLUME', 6X,F10.2)                                  00990000
      WRITE(6,105) V4                                                   01000000
 105  FORMAT(' ASYM VOLUME 2', 4X,F10.2)                                01010000
      VV=V1+V2+V3+V4                                                    01020000
      WRITE(6,115) VV                                                   01030000
 115  FORMAT(' TOTAL VOLUME ',14X,F10.2)                                01040000
      WRITE(6,106) S1                                                   01050000
 106  FORMAT(' BULK SURFACE ', 4X,F10.2)                                01060000
      WRITE(6,107) S2                                                   01070000
 107  FORMAT(' ASYM SURFACE ', 4X,F10.2)                                01080000
      SS=S1+S2                                                          01090000
      WRITE(6,116) SS                                                   01100000
 116  FORMAT(' TOTL SURFACE ',14X,F10.2)                                01110000
      WRITE(6,108) COUL1                                                01120000
 108  FORMAT(' COULOMB      ', 4X,F10.2)                                01130000
      WRITE(6,109) COUL2                                                01140000
 109  FORMAT(' COUL VOL RED ', 4X,F10.2)                                01150000
      WRITE(6,112) COUL5                                                01160000
 112  FORMAT(' COUL SURF RED', 4X,F10.2)                                01170000
      COUL=COUL1+COUL2+COUL5                                            01180000
      WRITE(6,117) COUL                                                 01190000
 117  FORMAT(' TOTAL COUL   ',14X,F10.2)                                01200000
      WRITE(6,113) CURV                                                 01210000
 113  FORMAT(' CURVATURE    ', 4X,F10.2)                                01220000
      WRITE(6,110) COUL3                                                01230000
 110  FORMAT(' COUL DIFFUS  ', 4X,F10.2)                                01240000
      WRITE(6,111) COUL4                                                01250000
 111  FORMAT(' COUL EXCHANGE', 4X,F10.2)                                01260000
         RETURN                                                         01270000
         END                                                            01280000
