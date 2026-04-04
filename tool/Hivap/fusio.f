      SUBROUTINE FUSIO (ELAB,AP,ZP,AT,ZT,Q2,V0,R0,D,XTH,APUSH,FPUSH,    00030005
     2            SIGR0,CRED,NORUTH,NOPROX,NOCURV,                      00040000
     1            ION,SIGML,FLGRAZ,LMAX,SIGF,ITEST,IOPT,MPT,CUTOFF)     00050000
C
C  
C   modifiziert am 8.7.1987 (new fusion) und am 8.10.1987 (limbar)
C   Stand: 3.5.1994
C    
C                      PARAMETERS ARE AS VLRO, AND FUSION               00060000
C     SIGR0 IS THE PERCENT CHANGE OF R0 (STAND.DEV.)                    00070000
      IMPLICIT REAL*8 (A-H,O-Z)
C--------------------------------------------------                     00250000
C
C     I.Giese 2007
C     Zeilen nach oben verschoben
C
      DIMENSION SIGML(1),SIGL(200),SIGL1(200)                           00260000
      COMMON/FUS/VBFUS,RBFUS,HWFUS                                      00270008
C------------------------------------------------                       00280000
C
C     MAYBE UNINITIALIZED: XDIFF
      DIMENSION NPOINT(7),KEY(8),Z(24),WEIGHT(24)                       00080000
      DATA NPOINT/2,3,4,5,6,10,15/,KEY/1,2,4,6,9,12,17,25/,             00090000
     1     Z/0.5773502,0.0,      0.7745966,                             00100000
     2       0.3399810,0.8611363,0.0,      0.5384693,                   00110000
     3       0.9061798,0.2386191,0.6612093,0.9324695,                   00120000
     4       0.1488743,0.4333953,0.6794095,0.8650633,                   00130000
     5       0.9739065,0.0,      0.2011940,0.3941513,                   00140000
     6       0.5709721,0.7244177,0.8482065,0.9372733,                   00150000
     7       0.9879925/,                                                00160000
     8WEIGHT/1.0000000,0.8888888,0.5555555,                             00170000
     9       0.6521451,0.3478548,0.5688888,0.4786286,                   00180000
     A       0.2369268,0.4679139,0.3607615,0.1713244,                   00190000
     B       0.2955242,0.2692667,0.2190863,0.1494513,                   00200000
     C       0.0666713,0.2025782,0.1984314,0.1861610,                   00210000
     D       0.1662692,0.1395706,0.1071592,0.0703660,                   00220000
     E       0.0307532/                                                 00230000
C                                                                       00240000
C------------------------------------------------                       00280000
C
C     I.Giese 2007
      SAVE
C
C     M-POINT GAUSS QUADRATURE,M=2,3,4,5,6,10,15;DEFAULT IS 5           00290000
      M=MPT                                                             00300000
      MTIMES=1                                                          00310000
      DO 7  L=1,7                                                       00320000
      IF(M.EQ.NPOINT(L)) GOTO 9                                         00330000
 7    CONTINUE                                                          00340000
      IF(MPT.LT.16) GOTO 8                                              00350000
      MTIMES=M/5                                                        00360003
      L=4                                                               00370000
      GOTO 9                                                            00380000
 8    M=10                                                              00390000
      L=6                                                               00400000
C     DEFAULT                                                           00410000
 9    JFIRST=KEY(L)                                                     00420000
      JLAST=KEY(L+1)-1                                                  00430000
      IF(CUTOFF.LE.0.) CUTOFF=2.5                                       00440000
C----------------------------------------------                         00450000
      LIMBAR=0                                                          00460009
      IF (LMAX.LT.0) LIMBAR=1                                           00470009
      R0FUS=0.                                                          00480000
      DO 40 L=1,200                                                     00490000
 40   SIGML(L)=0.                                                       00500000
      LMAX2=0                                                           00510000
      NOQ2=0                                                            00520000
      IF(ABS(Q2).LT.0.001) NOQ2=1                                       00530000
 50   ION1=ION                                                          00540000
      IF(NOQ2.EQ.1)ION1=4                                               00550000
      IF(SIGR0.GT.0.) GOTO 80                                           00560000
C----------------------------------------------                         00570000
      LMAX=0                                                            00580000
      IF(LIMBAR.EQ.1) LMAX=-1                                           00590009
      CALL FUSION(ELAB,AP,ZP,AT,ZT,Q2,V0,R0,D,XTH,APUSH,FPUSH,          00600000
     1            CRED,NORUTH,NOPROX,                                   00610008
     2            NOCURV,ION1,SIGML,FLGRAZ,LMAX,SIGF,ITEST,IOPT)        00620000
      RETURN                                                            00630000
C-----------------------------------------------                        00640000
 80   AP3=AP**0.333333                                                  00650000
      AT3=AT**0.333333                                                  00660000
      FLUCT=0.01*SIGR0                                                  00670000
      XLIM=CUTOFF*FLUCT                                                 00680000
      LMAX2=0                                                           00690000
C---------------------------------------------------------------------  00700000
C     RANGE OF INTEGRATION IS XIN TO XFI                                00710000
C--------------------------------------------------                     00720000
                                                                        00730000
 100  XIN=R0*(1.-XLIM)                                                  00740000
      XFI=R0*(1.+XLIM)                                                  00750000
      XXDIFF=(XFI-XIN)/2.                                               00760004
      DXIN=(XFI-XIN)/FLOAT(MTIMES)                                      00770000
      WIDTH=R0*FLUCT                                                    00780000
      XXAV=(XFI+XIN)/2.                                                 00790000
      SUM=0.                                                            00800000
      V0X=V0                                                            00810000
      R0X=R0                                                            00820000
      DX=D                                                              00830000
      QX=Q2                                                             00840000
      ION1=ION                                                          00850000
      IF(NOQ2.EQ.1)ION1=4                                               00860000
      XIN=XIN-DXIN                                                      00870000
C-------------                                                          00880003
      DO 300 MTIM=1,MTIMES                                              00890000
      XIN=XIN+DXIN                                                      00900000
      XFI=XIN+DXIN                                                      00910000
      XDIFF=(XFI-XIN)/2.                                                00920000
      XAV=(XFI+XIN)/2.                                                  00930000
C------------------------------------------                             00940000
      DO 200 J=JFIRST,JLAST                                             00950000
      IF(Z(J).NE.0.) GOTO  150                                          00960000
      R0X=XAV                                                           00970000
      XXX0=XXAV-XAV                                                     00980000
      GAU0=GAU(XXX0,WIDTH,XXDIFF)                                       00990004
      LMAX=0                                                            01000000
      IF(LIMBAR.EQ.1) LMAX=-1                                           01010009
      CALL FUSION(ELAB,AP,ZP,AT,ZT,QX,V0X,R0X,DX,XTH,APUSH,FPUSH,       01020000
     1            CRED,NORUTH,NOPROX,                                   01030008
     2            NOCURV,ION1,SIGL,FLGRAZ,LMAX,SIGF,ITEST,IOPT)         01040000
      SUM=SUM+WEIGHT(J)*SIGF*GAU0                                       01050000
      LMAX2=MAX0(LMAX2,LMAX)                                            01060000
      DO 140 L=1,LMAX2                                                  01070000
      SIGML(L)=SIGML(L)+                                                01080000
     1  WEIGHT(J)*SIGL(L)*GAU0                                          01090000
 140  CONTINUE                                                          01100000
      GOTO 200                                                          01110000
C------------------------                                               01120000
 150  XX1=Z(J)*XDIFF+XAV                                                01130000
      XXX1=XXAV-XX1                                                     01140000
      GAU1=GAU(XXX1,WIDTH,XXDIFF)                                       01150004
      R0X=XX1                                                           01160000
 160  LMAX=0                                                            01170000
      IF(LIMBAR.EQ.1) LMAX=-1                                           01180009
      CALL FUSION(ELAB,AP,ZP,AT,ZT,QX,V0X,R0X,DX,XTH,APUSH,FPUSH,       01190000
     1            CRED,NORUTH,NOPROX,                                   01200008
     2            NOCURV,ION1,SIGL,FLGRAZ,LMAX,SIGF1,ITEST,IOPT)        01210000
      LMAX2=MAX0(LMAX2,LMAX)                                            01220000
      XX2=-Z(J)*XDIFF+XAV                                               01230000
      XXX2=XXAV-XX2                                                     01240000
      GAU2=GAU(XXX2,WIDTH,XXDIFF)                                       01250004
      R0X=XX2                                                           01260000
 170  LMAX=0                                                            01270000
      IF(LIMBAR.EQ.1) LMAX=-1                                           01280009
      CALL FUSION(ELAB,AP,ZP,AT,ZT,QX,V0X,R0X,DX,XTH,APUSH,FPUSH,       01290000
     1            CRED,NORUTH,NOPROX,                                   01300008
     2            NOCURV,ION1,SIGL1,FLGRAZ,LMAX,SIGF2,ITEST,IOPT)       01310000
      SUM=SUM+WEIGHT(J)*(SIGF1*GAU1+SIGF2*GAU2)                         01320000
      LMAX2=MAX0(LMAX2,LMAX)                                            01330000
      DO 180 L=1,LMAX2                                                  01340000
      SIGML(L)=SIGML(L)+                                                01350000
     1  WEIGHT(J)*(SIGL(L)*GAU1+SIGL1(L)*GAU2)                          01360000
 180  CONTINUE                                                          01370000
 200  CONTINUE                                                          01380000
C------------------------------------------                             01390000
 300  CONTINUE                                                          01400003
      SIGF=XDIFF*SUM                                                    01410000
      LMAX2=MIN0(200,LMAX2)                                             01420000
      DO 210 L=1,LMAX2                                                  01430000
 210  SIGML(L)=SIGML(L) *XDIFF                                          01440000
C--------------------------------------------------                     01450000
      LMAX=LMAX2                                                        01460000
      RETURN                                                            01470000
      END                                                               01480000
C---------------------------------------------------------------------  01490000
      FUNCTION GAU(X,SI,XLIM)                                           01500000
C     SPECIAL GAUSS                                                     01510000
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(ABS(X).GT.XLIM .OR. XLIM.LE.0.) GOTO 100                       01520000
      IF(ABS(SI).LT.0.01*XLIM)GOTO 80                                   01530000
      SI2=2.*SI                                                         01540000
      DUM=X*X/(SI2  *SI)                                                01550000
      IF(DUM.GT.100.) GOTO 100                                          01560000
      T=XLIM/(1.41421*SI)                                               01570000
      T=ERF(T)                                                          01580000
      DUM=EXP(-DUM)                                                     01590000
      GAU=DUM/(2.5066*SI*T)                                             01600000
C     SQRT(2*PI)                                                        01610000
      RETURN                                                            01620000
 80   GAU=0.5/XLIM                                                      01630000
 100  GAU=0.                                                            01640000
      RETURN                                                            01650000
      END                                                               01660000
