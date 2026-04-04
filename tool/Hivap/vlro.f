      FUNCTION VLRO(P_AP,P_AT,ZP,ZT,P_Q2,L_NORUTH,L_NOPROX,L_NOCURV,    00020030
     1              OCOS,P_V0,R0,P_D,P_XTH,P_APUSH,FPUSH,ECM,           00030030
     2              L,R,L_IOP)                                          00040030
C
C     I.Giese 2007
C     Wegen Entries veraenderbare Parameter umsetzen
C
C
C  Berechnet Kern-Kern-Potential, verschiedene Eingaenge V0,VL,VR
C  letzte Modifikation 29.11.1983
C  Stand: 7.6.1994
C
C  Berechnung 'Extra-Push' aus x_eq_mean = (x_cn*x_eq_eff)**0.5
C  'Gaeggler - Feldmeier - Methode, Formeln cf. z.B. GSI-85-11 
C 
C     2 NOV 81  CORRECTED ERROR IN CURVATURE                            00050000
C     29 APR 83  INSERTED DEFAULTS COMPATIBLE WITH BLOCKI 77            00060014
C        AND CORRECTED ERROR RTS SCALING (AFFECTS OPT 1,2,4 IF Q2=0)    00070014
C     DEFORMED HI PROXIMITY POT                                         00080000
C     AP PROJECTILE MASS (SPHERIC)                                      00090000
C     AT TARGET MASS (DEFORMED)                                         00100000
C     ZP PROJECTILE Z  (SPHERIC)                                        00110000
C     ZT TARGET(DEFORMED) Z                                             00120000
C     Q2 QUADRUPOLE MOMENT OF TARGET IN (FM**2)                         00130000
C     NORUTH #0:DOES NOT USE CHANGE OF ORIENTATION ANGLE ALONG TRAJ.    00140025
C            =2:DOES NOT USE QUADRATIC COULOMB INSIDE                   00150025
C            =3:ONLY DOES COULOMB CORRECTION FOR FINITE Q2              00160025
C     NOPROX =1:USES CENTER LINE PRESCRIPTION(INSTEAD OF PROX.PRESCR.)  00170000
C     NOCURV =1:IGNORES CURVATURE DEPENDANCE OF POT                     00180000
C     OCOS COS OF   ANGLE          SYM.AXIS WITH CENTER LINE AT INFINITY00190000
C     V0,RO,D POT PARAMETRES                                            00200000
C     XTH,APUSH,FPUSH:PARAMETERS,XTH,A,F OF BJORNHOM-SWIATECKI 82 EXTRAP00210020
C     ECM CENTER-OF-MASS ENERGY                                         00220000
C     L ORBITAL ANG.MOM.                                                00230000
C     R CENTER DISTANCE(FM)                                             00240000
C     IOPT =-1:NUCLEAR IS WS,NO RADII CORRECTION ,USE CENTRAL RADII     00250025
C     IOPT =-2:SAME AS -1 BUT FOR AP=1 SPIN ORBIT L+1/2 IS ADDED        00260003
C     IOPT =0:VNUCLEAR IS EXPONENTIAL,NO RADII CORRECTION,USE CENTRAL RA00270025
C     IOPT =1:VNUCLEAR IS EXPONENTIAL,WITH RADII CORRECTION,USE SHARP RA00280025
C     IOPT =2:FULL PROXIMITY POT,WITH RADII CORRECTION,USE SHARP RADII  00290010
C     IOPT =3:FULL PROXIMITY POT,NO RADII CORRECTION,USE CENTRAL RADII  00300010
C     IOPT =4: EXP PROXIMITY POT,WITH RADII CORRECTION,USE SHARP RADII  00310010
C     IOPT =5: BASS POTENTIAL,BERLIN 79, SET D=1 OR 0,USE SHARP RADII   00320025
C------------------------------------------------------------------     00330000
C
C     I.Giese 2007
C     Statement vorgezogen
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     I.Giese 2007
C     Wegen Entry-Feststellung ueber BLOCK DATA
      COMMON/TESTS/L0_ENTRY,L1_ENTRY,L2_ENTRY
C     Wegen Entries folgende Variable abspeichern
      COMMON/VLROCOM/AP,AT,Q2,NORUTH,NOPROX,NOCURV,V0,D,XTH,APUSH,IOP,
     1               AL,AL2,ARED,C01,CENT,CENTL,COUL,DUM1,EPS2,EPS3,
     2               EXPUSH,EZP,FLCH2,IBASS,IOPT,IREP,OCOS1,PHI,PHIFU,
     3               PHIR,PUSHK,P2,Q4A,RA,RA2,RBAR,RC,RC2,RELL,RELL2,
     4               RPS,RTS,SINR,VQ2,VQ20,VQ4,X,X2,ZETA
      COMMON/PUSHPA/SPUSH,SCREIPA
C
C     I.Giese 2007
C     Wert 1...E-44 nur fuer RISC6000-Architektur
C     PC/Linux: 1.1754994E-38 kleinste sinnvolle Zahl
C
      DATA ESQUAR,HBAR,AMU/1.4399, 0.65819E-21, 1.03641E-44/            00340012
C
C     I.Giese 2007
C     Hier kein SAVE notwendig, da alle Variable in VLROCOM aufgelistet
C
C     I.Giese 2007
C     Lokale Kopien der benutzten Parameter innerhalb ENTRY-Funktionen
C     Festhalten, ueber welchen Entry man in VLRO kam
C
      L0_ENTRY = 1
      L1_ENTRY = 0
      L2_ENTRY = 0
      AP = P_AP
      AT = P_AT
      Q2 = P_Q2
      NORUTH = L_NORUTH
      NOPROX = L_NOPROX
      NOCURV = L_NOCURV
      V0 = P_V0
      D = P_D
      XTH = P_XTH
      APUSH = P_APUSH
      IOP = L_IOP
C     Ende Umspeichern der veraenderbaren Groessen der Parameterliste
C
C------------------------------------------------------------------     00350011
      IBASS=1                                                           00360015
      IF(IOP.NE.5) IBASS=0                                              00370015
      IOPT=IOP                                                          00380015
      IF(IBASS.NE.1) GOTO 5                                             00390026
      IOPT=4                                                            00400026
      V0=1.                                                             00410026
      D=1.                                                              00420026
      NOCURV=0                                                          00430026
      NORUTH=1                                                          00440026
 5    ARED=AP*AT/(AP+AT)                                                00450026
      EZP=ESQUAR*ZP                                                     00460000
      COUL=EZP*ZT                                                       00470000
      C01=ARED*AMU*COUL/HBAR                                            00480000
      CENT=0.5*HBAR*HBAR/(ARED*AMU)                                     00490000
      AP3=AP**0.333333                                                  00500020
      AT3=AT**0.333333                                                  00510020
      DUM=(AP+AT-2.*(ZP+ZT))/(AP+AT)                                    00520020
      DEPI= 1.-1.7826*DUM*DUM                                           00530020
C                                                                       00540020
      EXPUSH=0.                                                         00550020
      IF(XTH.LE.0.) GOTO 7                                              00560020
      ZETACR=50.883*DEPI                                                00570020
      DEPA=AP3*AT3*(AP3+AT3)                                            00580020
      IF(SPUSH.GT.0) GOTO 6
      ZETA=4.*ZP*ZT/DEPA-ZETACR*XTH                                     00590020
      PUSHK=7.601E-4*DEPA*(AP3+AT3)*APUSH*APUSH/(AP+AT)                 00600020
      IF(ZETA.GT.0.) EXPUSH=PUSHK*ZETA*ZETA                             00610020
      GOTO 8
 6    AMAPA=DEPA*(AP+AT)*(AP+AT)    
      ZETAF=4*AP*AT*(ZP+ZT)*(ZP+ZT)/AMAPA 
      ZSPALT=(ZP+ZT)*(ZP+ZT)/(AP+AT)
      ZETFM=SQRT(ZETAF*ZSPALT)
      ZETA=ZETFM-ZETACR*XTH
      PUSHK=7.601E-4*DEPA*(AP3+AT3)*APUSH*APUSH/(AP+AT)                 00600020
      IF(ZETA.GT.0.) EXPUSH=PUSHK*ZETA*ZETA                
      IF(SCREIPA.GT.0.) GOTO 8
      SCREIPA=1
      XEQEFF=ZETAF/ZETACR
      XCN=ZSPALT/ZETACR
      XEQM=ZETFM/ZETACR  
      ECHPU=7.601E-4*DEPA*(AP3+AT3)*ZETACR*ZETACR/(AP+AT)  
      WRITE(6,1001)EXPUSH,XEQM,XTH,APUSH,ECHPU,XCN,XEQEFF
 1001 FORMAT(/'Extra-Push Parameter: E_push =',F8.2,' MeV',3X,
     a'x_eq_m=',F6.3,3x,'xth=',F6.3,3x,'Apush=',F6.3,3X,'E_ch=',
     bF8.3,3X,'xcn=',F6.3,3X,'x_eq_eff=',F6.3) 
 8    ALCH=0.10270*DEPA*AP3*AT3/SQRT(AP+AT)                             00620020
      FLCH=FPUSH/ALCH                                                   00630020
      FLCH2=FLCH*FLCH                                                   00640020
 7    IF(XTH.LT.0.) EXPUSH=APUSH                                        00650024
C                                                                       00660011
      IF(V0.GT.0.) GOTO 9                                               00670024
      DUM=0.9517*DEPI                                                   00680020
      IF(IOPT.GT.1)V0=12.5664   *DUM                                    00690013
      IF(IOPT.EQ.0 .OR. IOPT.EQ.1 ) V0=43.1906*DUM                      00700013
 9    IF(IOPT.GT.1 .AND. D.LE.0.) D=1.                                  00710011
      IF(IOPT.LT.2 .AND. D.LE.0.) D=0.75                                00720011
      RTS=R0*AT3                                                        00730000
      RPS=R0*AP3                                                        00740000
      IF(AP.LT.4.01 .AND. AP.GT.3.99) RPS=1.6                           00750000
      IF(RTS.GT.0.) GOTO 10                                             00760018
      RPS=1.28*AP3-0.76+0.8/AP3                                         00770011
      RTS=1.28*AT3-0.76+0.8/AT3                                         00780011
      IF(IOPT.EQ. 1 .OR.IOPT.EQ.2 .OR. IOPT.EQ.4) GOTO 9000             00790023
      RPS=RPS- 1./RPS                                                   00800013
      RTS=RTS- 1./RTS                                                   00810013
 9000 R0=(RPS+RTS)/(AP3+AT3)                                            00820023
 10   RA=RTS                                                            00830023
      RC=RTS                                                            00840000
      RBAR=RTS*RPS/(RTS+RPS)                                            00850000
      Q2A=0.                                                            00860000
      Q4A=0.                                                            00870000
      IF(ABS(Q2).LT.0.01) GOTO 25                                       00880025
      RTQ=1.15*AT3                                                      00890000
      Q2A=0.5*Q2                                                        00900000
      VQ20=EZP*Q2A                                                      00910000
      IF(NORUTH.EQ.3) GOTO 25                                           00920025
C-------------------------------                                        00930000
C     AXES OF ELLIPSOID                                                 00940000
      DUM=2.5*Q2/ZT                                                     00950000
      RQ3=RTQ*RTQ*RTQ                                                   00960000
      RC0=RTQ                                                           00970000
      IREP=0                                                            00980000
 15   RC1=DUM+RQ3/RC0                                                   00990000
      RC1=SQRT(RC1)                                                     01000000
      DIFF=RC1-RC0                                                      01010000
      IF(ABS(RC1-RC0).LT.0.01) GOTO 18                                  01020000
      IREP=IREP+1                                                       01030000
      IF(IREP.GT.20) GOTO 18                                            01040000
      RC0=RC1                                                           01050000
      GOTO 15                                                           01060000
 18   RC=0.5*(RC1+RC0)                                                  01070000
      RA=SQRT(RQ3/RC)                                                   01080000
      RC=RTS*RC/RTQ                                                     01090000
      RA=RTS*RA/RTQ                                                     01100000
 25   IF(IOPT.LE.0) GOTO 26                                             01110000
      IF(IOPT.EQ.4) GOTO 28                                             01120000
      IF(IOPT.GT.2) GOTO 26                                             01130000
C     RADII CORRECTION FOR IOPT=1,2,4                                   01140000
 28   DSQ=D*D                                                           01150000
      RPS=RPS-DSQ/RPS                                                   01160000
      RTS=RTS-DSQ/RTS                                                   01170014
      RA=RA-DSQ/RC                                                      01180000
      RC=RC-DSQ*RC/(RA*RA)                                              01190000
 26   RA2=RA*RA                                                         01200000
      RELL2=RA2                                                         01210000
      RELL=RA                                                           01220000
      RC2=RC*RC                                                         01230000
      EPS2=(RC2-RA2)/RC2                                                01240000
      EPS3=EPS2*0.5/RA2                                                 01250000
 30   CONTINUE                                                          01260000
C-------------------------------------------------                      01270000
      ENTRY VO(L,R,ECM,OCOS)                                            01280000
C
C     I.Giese 2007
      L1_ENTRY = 1
C
      IF(ABS(Q2).LT.0.01) GOTO 40                                       01290025
      IF(OCOS.NE.OCOS1)PHI=ACOS(OCOS)                                   01300000
      OCOS1=OCOS                                                        01310000
      IF(NORUTH.EQ.0) GOTO 40                                           01320000
      PHIR=PHI                                                          01330000
      X=OCOS                                                            01340000
      SINR=SIN(PHIR)                                                    01350000
      X2=X*X                                                            01360000
      P2=1.5*X2-0.5                                                     01370000
      VQ2=VQ20*P2                                                       01380000
      IF(Q4A.EQ.0.) GOTO 40                                             01390000
      P4=4.3750*X2*X2-3.75*X2+0.375                                     01400000
      VQ4=EZP*Q4A*P4                                                    01410000
 40   CONTINUE                                                          01420000
C---------------------------------------------------                    01430000
C
C     I.Giese 2007
C
C     WIRD NICHT BENUTZT...
C     ENTRY VL(L,R,ECM,OCOS)                                            01440000
      AL=L                                                              01450000
      AL2=AL*(AL+1.)                                                    01460000
      CENTL=CENT*AL2                                                    01470000
C                                                                       01480020
      EXPUSH=0.                                                         01490020
      IF(XTH.LE.0.) GOTO 50                                             01500020
      DUM=ZETA+FLCH2*AL2                                                01510020
      IF(DUM.GT.0.) EXPUSH=PUSHK*DUM*DUM                                01520020
 50   IF(XTH.LT.0.) EXPUSH=APUSH                                        01530024
C                                                                       01540020
      IF(ABS(Q2).LT.0.01) GOTO 150                                      01550025
      IF(OCOS.NE.OCOS1)PHI=ACOS(OCOS)                                   01560000
      OCOS1=OCOS                                                        01570000
      IF(NORUTH.NE.0) GOTO 150                                          01580000
      CHIAS=0.                                                          01590000
      IF(L.EQ.0) GOTO  80                                               01600000
      C02=2.*ARED*AMU*ECM                                               01610000
      C01L=C01/AL                                                       01620000
      C03L=SQRT(C02+C01L*C01L)                                          01630000
      DUM2=C01L/C03L                                                    01640000
      CHIAS=ACOS(DUM2)                                                  01650000
 80   PHIR=PHI-CHIAS                                                    01660000
      X=COS(PHIR)                                                       01670000
      SINR=SIN(PHIR)                                                    01680000
      X2=X*X                                                            01690000
      P2=1.5*X2-0.5                                                     01700000
      VQ2=VQ20*P2                                                       01710000
      IF(Q4A.EQ.0.) GOTO 150                                            01720000
      P4=4.3750*X2*X2-3.75*X2+0.375                                     01730000
      VQ4=EZP*Q4A*P4                                                    01740000
C----------------------------------------------------                   01750000
      ENTRY VR(R)                                                       01760000
C
C     I.Giese 2007
      L2_ENTRY = 1
C                                                                       01770000
 150  R2=R*R                                                            01780000
      VCENT=CENTL/R2                                                    01790000
      VCOUL=COUL/R                                                      01800000
      IF(ABS(Q2) .LT.0.01) GOTO 280                                     01810025
      R3=R2*R                                                           01820000
      VCOUL=VCOUL+VQ2/R3                                                01830000
      IF(Q4A.EQ.0.) GOTO 200                                            01840000
      R5=R3*R2                                                          01850000
      VCOUL=VCOUL+VQ4/R5                                                01860000
C                                                                       01870000
 200  IF(ABS(Q2).LT.0.01 .OR. NORUTH.EQ.3) GOTO 280                     01880025
      PHIMIN=PHIR                                                       01890000
      COSMIN=X                                                          01900000
      IREP=0                                                            01910000
      IF(NOPROX.EQ.0) GOTO 220                                          01920000
C--------------------------------------------------------               01930000
C     GET MIN.DIST. 'S0' AND ANGLE 'PHIMIN'                             01940000
      RELL2=RA2/(1.-EPS2*X2)                                            01950000
      RELL=SQRT(RELL2)                                                  01960000
      GOTO 230                                                          01970000
 225  COSMIN=COS(PHIMIN)                                                01980000
 220  DUM2=COSMIN*COSMIN                                                01990000
      RELL2=RA2/(1.-EPS2*DUM2)                                          02000000
      RELL=SQRT(RELL2)                                                  02010000
      IREP=IREP+1                                                       02020000
      IF(IREP.EQ.2) GOTO 230                                            02030000
      DUM=2.*PHIMIN                                                     02040000
      DELTA=(R-RELL)*RELL2*SIN(DUM)*EPS3/R                              02050000
      IF(ABS(DELTA).GT.ABS(PHIR))DELTA=0.                               02060000
      PHIMIN=PHIR-DELTA                                                 02070000
      GOTO 225                                                          02080000
 230  XELL=RELL*COSMIN                                                  02090000
      YELL=RELL*SIN(PHIMIN)                                             02100000
      DUM=R*X        -XELL                                              02110000
      DU1=R*SINR     -YELL                                              02120000
      S0=SQRT(DUM*DUM+DU1*DU1)-RPS                                      02130000
      GOTO 290                                                          02140000
 280  S0=R-RPS-RTS                                                      02150000
 290  IF(NORUTH.EQ.2)GOTO 291                                           02160000
      RTOUCH=RELL+RPS                                                   02170000
      IF(S0.LT.0.)VCOUL=VCOUL-COUL/R +1.5*COUL*(1.-0.333333*(R/RTOUCH)**02180000
     1 2.)/RTOUCH                                                       02190000
 291  IF(NOCURV.EQ.1) GOTO 310                                          02200000
C---------------------------------------------------                    02210000
C     GET REDUCED CURVATURE  ELLIPSOID(RC,RA) + SPHERE(RPS)             02220000
      IF(ABS(Q2).LE.0.01 .OR. NORUTH.EQ.3) GOTO 305                     02230025
      COS2V=(RELL2-RA2)/(RC2-RA2)                                       02240002
      SIN2V=1.-COS2V                                                    02250002
      DUM=RA2*COS2V+RC2*SIN2V                                           02260002
      DUM=ABS(DUM)                                                      02270002
      R1R2=DUM*DUM/RC2                                                  02280002
      R1PR2=(RA2+DUM)*SQRT(DUM)/(RA*RC)                                 02290002
      GOTO 306                                                          02300007
 305  R1R2=RTS*RTS                                                      02310007
      R1PR2=RTS+RTS                                                     02320007
 306  RBAR2=R1R2*RPS*RPS/(R1R2+RPS*R1PR2+RPS*RPS)                       02330007
      RBAR=SQRT(RBAR2)                                                  02340000
      GOTO 320                                                          02350000
 310  RBAR=1.                                                           02360000
C------------------------------------------------------                 02370000
 320  IF(IOPT.LT.2) GOTO 325                                            02380000
      ZETA=S0/D                                                         02390000
      IF(IOPT.EQ.4) GOTO 321                                            02400000
      IF(ZETA.LT.1.2511)GOTO 322                                        02410000
 321  IF(IBASS.EQ.0) PHIFU=-3.437*EXP(-ZETA/0.75)                       02420015
      IF(IBASS.EQ.1) PHIFU=-1./(.033*EXP(ZETA/3.5)+0.007*EXP(ZETA/.65)) 02430016
      GOTO 323                                                          02440000
 322  DUM=ZETA-2.54                                                     02450000
      DUM2=DUM*DUM                                                      02460000
      PHIFU=-0.5*DUM2-0.0852*DUM*DUM2                                   02470000
 323  VN=V0*D*RBAR*PHIFU                                                02480000
      GOTO 328                                                          02490000
 325  IF(IOPT.EQ.(-1))GOTO 335                                          02500000
      IF(IOPT.EQ.(-2))GOTO 335                                          02510002
      DUM=S0/D                                                          02520000
      VN=0.                                                             02530000
      IF(DUM.LT.150.) VN=-V0*RBAR*EXP(-DUM)                             02540004
 328  VLRO=VCOUL+VCENT+VN+EXPUSH                                        02550020
C
C     I.Giese 2007
C
C     Parameter-Rueckgabe, falls veraendert (4) und vom Haupt-Entry aus
      IF(L0_ENTRY.NE.1) GOTO 4000
C     P_AP = AP
C     P_AT = AT
C     P_Q2 = Q2
      L_NORUTH = NORUTH
C     L_NOPROX = NOPROX
      L_NOCURV = NOCURV
      P_V0 = V0
      P_D = D
C     P_XTH = XTH
C     P_APUSH = APUSH
C     L_IOP = IOP
 4000 L0_ENTRY = 0
      L1_ENTRY = 0
      L2_ENTRY = 0
C
      RETURN                                                            02560000
 335  DUM=S0/D                                                          02570000
      VN=0.                                                             02580000
      IF(DUM.LT.150.) VN=-V0*RBAR/(1.+EXP(DUM))                         02590004
      IF(IOPT.NE.(-2)) GOTO 328                                         02600002
C     FOR TEST PURPOSE ADD SPIN ORBIT FOR PROTONS J=L+1/2               02610002
C     BECCHETTI 1969                                                    02620002
      IF(AP.GT.1.01) GOTO 328                                           02630002
      DUM=(R-1.01*AT**(1./3))/0.75                                      02640005
      VSO=0.                                                            02650002
      IF(DUM.LT.80.) DUM1=1.+EXP(DUM)                                   02660006
      IF(DUM.LT.80.) VSO=-6.2*AL*2.0*EXP(DUM)/(0.75*R*DUM1*DUM1)        02670006
      VN=VN+VSO                                                         02680002
      GOTO 328                                                          02690000
      END                                                               02700000
