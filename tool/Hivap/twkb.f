      FUNCTION TWKB(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,V0,R0,D,        00010028
     1           XTH,APUSH,FPUSH,DR,OCOS,ECM,L,RIN,ROUT,RFUS,CRED,IOPT) 00020028
C
C  Berechnet Transmissionskoeffizienten fuer Fusion mittels WKB - 
C  naeherung ( Eingang WKB )
C  letzte Modifizierung 10.8.1987
C  Stand: 3.5.1994
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NPOINT(7),KEY(8),Z(24),WEIGHT(24)                       00040004
      DATA NPOINT/2,3,4,5,6,10,15/,KEY/1,2,4,6,9,12,17,25/,             00050004
     1     Z/0.5773502  ,0.0000000  ,0.7745966  ,                       00060004
     2       0.3399810  ,0.8611363  ,0.0000000  ,0.5384693  ,           00070004
     3       0.9061798  ,0.2386191  ,0.6612093  ,0.9324695  ,           00080004
     4       0.1488743  ,0.4333953  ,0.6794095  ,0.8650633  ,           00090004
     5       0.9739065  ,0.0000000  ,0.2011940  ,0.3941513  ,           00100004
     6       0.5709721  ,0.7244177  ,0.8482065  ,0.9372733  ,           00110004
     7       0.9879925  /,                                              00120004
     8WEIGHT/1.0000000  ,0.8888888  ,0.5555555  ,                       00130004
     9       0.6521451  ,0.3478548  ,0.5688888  ,0.4786286  ,           00140004
     A       0.2369268  ,0.4679139  ,0.3607615  ,0.1713244  ,           00150004
     B       0.2955242  ,0.2692667  ,0.2190863  ,0.1494513  ,           00160004
     C       0.0666713  ,0.2025782  ,0.1984314  ,0.1861610  ,           00170004
     D       0.1662692  ,0.1395706  ,0.1071592  ,0.0703660  ,           00180004
     E       0.0307532  /                                               00190004
      DATA HBARC,AMU/197.329,931.502/                                   00200030
C-----------------------------------------------------------------      00210004
C
C     I.Giese 2007
      SAVE
C
      AMXEXP=100.                                                       00220018
      MTIMES=2                                                          00230027
      MPNT=5                                                            00240004
      DO  5 I=1,7                                                       00250004
      IF(MPNT.EQ.NPOINT(I)) GOTO 6                                      00260004
 5    CONTINUE                                                          00270004
      MPNT=5                                                            00280004
      I=4                                                               00290006
C     DEFAULT                                                           00300004
 6    JFIRST=KEY(I)                                                     00310004
      JLAST=KEY(I+1)-1                                                  00320004
C------------------------------------------------------------------     00330004
      IF(DR.LE.0.)DR=0.2                                                00340004
      IUP=4./DR                                                         00350004
      IUPP=8./DR                                                        00360029
      IF(CRED.LE.0)CRED=1.                                              00370011
      ARED=AP*AT*CRED/(AP+AT)                                           00380030
      F0=2.*ARED*AMU/(HBARC*HBARC)                                      00390030
      F0=SQRT(F0)                                                       00400011
      R=RIN                                                             00410004
      V=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,          00420028
     1       XTH,APUSH,FPUSH,ECM,L,R,IOPT)                              00430028
      GOTO 8                                                            00440006
C-----------------------------------------------------------------      00450006
C
C     I.Giese 2007
C     WIRD ANSCHEINEND NICHT BENUTZT
C
C     ENTRY WKB(OCOS,ECM,L,RIN,ROUT)                                    00460004
      R=RIN                                                             00470004
      V=VO(L,RIN,ECM,OCOS)                                              00480006
 8    DO 10 I=1,IUP                                                     00490006
      R1=R                                                              00500004
      V1=V                                                              00510004
      R=R-DR                                                            00520004
      V=VR(R)                                                           00530004
      IF(V.LT.ECM) GOTO 12                                              00540007
 10   CONTINUE                                                          00550004
      IF(L.LT.3)GOTO 200                                                00560015
      RI=DMAX1(RFUS,RIN)                                                00570017
      GOTO 25                                                           00580015
 12   R3=0.5*(R+R1)                                                     00590007
      V2=VR(R3)                                                         00600004
      IF(V2.LT.ECM)GO TO 15                                             00610004
      V1=V2                                                             00620004
      R1=R3                                                             00630004
      GOTO 20                                                           00640004
 15   V=V2                                                              00650004
      R=R3                                                              00660004
 20   RI=((V1-ECM)*R-(V-ECM)*R1)/(V1-V)                                 00670004
 25   R=ROUT                                                            00680015
      V=VR(ROUT)                                                        00690004
      DO 30 I=1,IUPP                                                    00700004
      R1=R                                                              00710004
      V1=V                                                              00720004
      R=R+DR                                                            00730004
      V=VR(R)                                                           00740004
      IF(V.LT.ECM) GOTO 35                                              00750004
 30   CONTINUE                                                          00760004
      GOTO 201                                                          00770004
 35   R3=0.5*(R+R1)                                                     00780004
      V2=VR(R3)                                                         00790004
      IF(V2.LT.ECM)GO TO 45                                             00800004
      V1=V2                                                             00810004
      R1=R3                                                             00820004
      GOTO 50                                                           00830004
 45   V=V2                                                              00840004
      R=R3                                                              00850004
 50   RO=((V1-ECM)*R-(V-ECM)*R1)/(V1-V)                                 00860004
      RIN  =RI                                                          00870005
      ROUT =RO                                                          00880005
      IF(RFUS.LT.RO)RI=DMAX1(RI,RFUS)                                   00890011
C----------------------------------------------------------             00900001
      DRR=(RO-RI)/FLOAT(MTIMES)                                         00910004
      R1=RI-DRR                                                         00920004
      SY=0.                                                             00930001
      DO 95 I=1,MTIMES                                                  00940004
      R1=R1+DRR                                                         00950004
      R2=R1+DRR                                                         00960004
      C=(R2-R1)/2.                                                      00970004
      D1=(R2+R1)/2.                                                     00980004
      SUM=0.                                                            00990004
      DO 90 J=JFIRST,JLAST                                              01000004
      IF(Z(J).NE.0.) GOTO 80                                            01010004
      DUM=VR(D1)-ECM                                                    01020004
      IF(DUM.GT.0.)DUM=SQRT(DUM)                                        01030004
      IF(DUM.LT.0.)DUM=0.                                               01040004
      SUM=SUM+WEIGHT(J)*DUM                                             01050004
      GOTO 90                                                           01060004
 80   DUM=Z(J)*C+D1                                                     01070004
      DUM1=-Z(J)*C+D1                                                   01080004
      DUM=VR(DUM)-ECM                                                   01090004
      DUM1=VR(DUM1)-ECM                                                 01100004
      IF(DUM.GT.0.)DUM=SQRT(DUM)                                        01110004
      IF(DUM.LE.0.)DUM=0.                                               01120004
      IF(DUM1.GT.0.)DUM1=SQRT(DUM1)                                     01130004
      IF(DUM1.LE.0.)DUM1=0.                                             01140004
                     SUM=SUM+WEIGHT(J)*(DUM+DUM1)                       01150004
 90   CONTINUE                                                          01160004
      GAUSS=C*SUM*F0                                                    01170004
      SY=SY+GAUSS                                                       01180004
 95   CONTINUE                                                          01190004
      Y=SY                                                              01200001
      TWKB=0.                                                           01210001
C----------------------------------------------------------             01220001
      DUM=2.*Y                                                          01230018
      TWKB=0.                                                           01240018
      IF(DUM.LT.AMXEXP)TWKB=EXP(-DUM)                                   01250018
      RETURN                                                            01260001
C----------------------------------------------------------             01270006
 200  IF(L.LE.2) WRITE(6,210)  L,OCOS,RIN,R,V,ECM,IUP,DR                01280029
 210  FORMAT(' TWKB TROUBLE RI  L,OCOS=',I4,F8.3/' RIN,R,V,ECM',4F9.3/  01290009
     1       ' IUP,DR',I6,F8.2)                                         01300009
      TWKB=0.                                                           01310029
      RETURN                                                            01320006
 201  IF(L.LE.2) WRITE(6,220)  L,OCOS,ROUT,R,V,ECM                      01330029
 220  FORMAT(' TWKB TROUBLE RO  L,OCOS=',I4,F8.3/' ROUT,R,V,ECM',4F9.3) 01340008
      TWKB=0.                                                           01350029
      RETURN                                                            01360006
      END                                                               01370001
