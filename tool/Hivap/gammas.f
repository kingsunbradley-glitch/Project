      SUBROUTINE  GAMMAS(ACN,IOUT,DEL,IGAM,NOE2,STRIPE,IOPT)            00020005
C
C  Berechnet Gamma - Staerkefunktionen
C  modifiziert am 1.6.1982 und am 25.4.1983
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/GAM/KEG,ECUBE(26),EFIVE(26),CGIANT,EGIANT,GAMMA,GAMQ,JFACTR00050014
      DIMENSION ARRAY(10)                                               00060012
C
C     I.Giese 2007
      SAVE
C
C     GAMMA NORMALISATION                                               00070000
C--------------------------------------------------------------         00080005
C     DEFAULTS  E2:TEN WEISSKOPF UNITS, 1 TO 4 MEV                      00090019
C               E1:GIANT DIPOLE RESONANCE 1 TO (EGIANT+WGIANT) MEV      00100005
      DEL=1.                                                            00110025
      IOPT=0                                                            00120023
      STRIPE=0.                                                         00130023
      JFACTR=0                                                          00140022
      EG1WU= 0.                                                         00150021
      EG2WU=10.                                                         00160021
      EG2MIN=1.                                                         00170018
      EG2MAX=4.                                                         00180005
      CGIANT=  11.5/3.                                                  00190007
      EGIANT=80./(ACN**0.33333)                                         00200005
      WGIANT=5.                                                         00210005
      DUM=EGIANT+WGIANT                                                 00220005
      EG1MIN= 1.                                                        00230018
      EG1MAX=-1.                                                        00240010
      DO 10 KE=1,KEG                                                    00250014
      ECUBE(KE)=0.                                                      00260014
 10   EFIVE(KE)=0.                                                      00270014
      IF(IGAM.EQ.0)GOTO 60                                              00280005
C---------------------------------------------------------------        00290005
      CALL MYRD(ARRAY,7,23,5)                                           00300020
      EG1WU =ARRAY(1)                                                   00310012
      EG2WU =ARRAY(2)                                                   00320012
      EG1MIN=ARRAY(3)                                                   00330012
      EG1MAX=ARRAY(4)                                                   00340012
      EG2MIN=ARRAY(5)                                                   00350012
      EG2MAX=ARRAY(6)                                                   00360012
      JFACTR=ARRAY(7) + 0.01                                            00370014
      CALL MYRD(ARRAY,5,23,5)                                           00380020
      CGIANT=ARRAY(1)                                                   00390012
      EGIANT=ARRAY(2)                                                   00400012
      WGIANT=ARRAY(3)                                                   00410012
      STRIPE=ARRAY(4)                                                   00420012
      IOPT  =ARRAY(5)+0.01                                              00430012
      IF(EG2WU.LE.0.)NOE2=1                                             00440005
      IF(EG2MAX.LE.0.)EG2MAX=4                                          00450005
      IF(EGIANT.LE.0.)EGIANT=80/(ACN**0.33333)                          00460005
      IF(WGIANT.LE.0.)WGIANT=5.                                         00470005
C---------------------------------------------------------------        00480016
 60   IDEL=DEL+0.01                                                     00490018
      DEL=IDEL                                                          00500018
      IF(DEL.LE.1.) DEL=1.                                              00510018
      IF(KEG.EQ.0)KEG=26                                                00520002
C---------------------------------------------------                    00530005
      IF(CGIANT.LE.0.) GOTO 180                                         00540005
C                                                                       00550000
C**************************** LORENTZIAN ****************************** 00560000
 160  CONTINUE                                                          00570000
      EG2=EGIANT*EGIANT                                                 00580000
      W2=WGIANT*WGIANT                                                  00590000
      CONST=(1.63E-6)*CGIANT*ACN*WGIANT                                 00600009
      FK=DEL                                                            00610004
      DO 163 K=1,KEG                                                    00620002
      FK2=FK*FK                                                         00630000
      EFIVE(K)=0.                                                       00640000
      DUM=EG2-FK2                                                       00650000
      ECUBE(K)=CONST*FK2*FK2/(DUM*DUM+FK2*W2)                           00660000
      FK=FK+1.                                                          00670004
 163  CONTINUE                                                          00680004
C-------------------------------------------------------------------    00690006
 180  GNORM=(6.75E-8)*(ACN**0.6667)*6.283*EG1WU                         00700007
      QNOR= (4.77E-14)*(ACN**1.3333)*6.283*EG2WU                        00710007
 95   CONTINUE                                                          00720006
      FK=DEL                                                            00730006
      DO 80 K=1,KEG                                                     00740006
      IF(FK.LT.EG1MIN .OR. FK.GT.EG1MAX)GOTO 70                         00750006
      ECUBE(K)=GNORM*FK**3                                              00760006
 70   IF(FK.LT.EG2MIN .OR. FK.GT.EG2MAX)GOTO 80                         00770006
      EFIVE(K)=QNOR*FK**5                                               00780006
 80   FK=FK+1.                                                          00790006
      FKSI1=GNORM/6.283                                                 00800006
      FKSI2=QNOR/6.283                                                  00810006
C--------------------------------------------------------------         00820005
      IF(IOUT.EQ.0) THEN                                                00830026
         WRITE(6,230) (ECUBE(I),I=1,20)                                 00840026
 230     FORMAT(/' GAMMAS:   ECUBE'/4(5E12.3/))                         00850026
      ENDIF                                                             00860026
      RETURN                                                            00870000
      END                                                               00880000
