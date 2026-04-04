      SUBROUTINE SEPENF(MTGT,IZT,MPROJ,IZP,QVALUE,NA,NZ,POP,KMAX1,      00020045
     1                 ISHELL,SHELLS,SHELL0,DELTAS,IPAIR,KHS)           00030037
C
C Berechnet Q-Werte, Separationsenergien fuer n,p,alphas, Spalt-
C barrieren, Schalenkorrektur-Energien
C ruft MSBEN, FISROT, liest Massentafel ( von logischer Einheit 9 )
C Einfuehrung experimenteller Paarungskorrekturen IPAIR=4 (6.10.1981)
C Einfuehrung gg-Kerne als Referenz fuer Schalenkorrekturen
C IPAIR=4 (18.12.1981)
C Stand: 5.5.1994
C
C     NEEDED EXTERNAL BE,MP,MC,NA,NZ,EXC,MTGT,IZT,MPROJ,IZP,QVALUE,     00060000
C     POP(DIM ABOUT 100),IND(DIM=15),ZCN,ACN,NUMB,MASSES,ELEMNT         00070000
C     'MASSES'  LOG UNIT MASS TABLE                                     00080000
C                                                                       00090000
C******* GET SEPARATION ENERGIES N,P,ALPHA,D,T *************************00100000
C                                                                       00110000
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IND(27),POP(100),SHELLS(1),DELTAS(1),SMAS(27,27)        00120041
      DIMENSION MAGIC(6),MAGICN(6)                                      00130037
      COMMON/YR/JDIM,EJAY(500),JUPYR,IYR,IRAST,RATIOS(6),               00140000
     1          RATIOF(6),INERF,FINERT,FINERF,NOROTF                    00150000
      COMMON/BIND/BE(25,25,5),EXC(27,27),MC,MP,NUMB,MASSES,DELT         00160014
      COMMON/ELEM/ELEMNT(128)                                           00170000
      COMMON/FISRT/ROT0,A2MS,RKAPPA,X0,Y0,RKAPA4                        00180000
      COMMON/BARR/BFLDM(25,25),BF(25,25),BAR0,BARFAC,IBF,IFISRT         00190000
C     MAYBE UNINITIALIZED: INIT,R0MS,C3MS,EXCCNT,ODDEV
      CHARACTER*4 ELEMNT                                                00200040
      DATA MAGIC/114,82,64,50,28,20/                                    00210037
      DATA MAGICN/184,126,82,50,28,20/                                  00220037
C
C     I.Giese 2007
      SAVE
C
      NZMX=25                                                           00230000
      NAMX=25                                                           00240000
      KAKZMX=25*25                                                      00250000
      NZMX2=NZMX+2                                                      00260010
      NAMX2=NAMX+2                                                      00270010
      NA=MIN0(NAMX-1,NA)                                                00280010
      NZ=MIN0(NZMX-1,NZ)                                                00290010
C*********************************************************************  00300000
C                                                                       00310000
      IF(INIT.EQ.99) GOTO 15                                            00320027
      A2MS=17.9439                                                      00330027
      R0MS=1.2249                                                       00340027
      C3MS=0.863987/R0MS                                                00350027
      RKAPPA=1.7826                                                     00360027
      NOROTF=0                                                          00370046
      FINERF=1.                                                         00380046
      RKAPA4=0.                                                         00390046
 15   WRITE(6,106)                                                      00400027
 106  FORMAT(/' Unterprogramm SEPENF')                                  00410000
      IF(IPAIR.NE.4) GOTO 84                                            00420001
      DO 80 KAKZ=1,KAKZMX                                               00430000
 80   DELTAS(KAKZ)=0.                                                   00440000
C-------------------------------------                                  00450000
 84   IF(IFISRT.NE.1) GOTO 100                                          00460001
      CALL MYRD(POP,7,23,5)                                             00470041
      A2MS=POP(1)                                                       00480000
      R0MS=POP(2)                                                       00490000
      C3MS=POP(3)                                                       00500000
      RKAPPA=POP(4)                                                     00510000
      NOROTF=POP(5)                                                     00520000
      FINERF=POP(6)                                                     00530000
      RKAPA4=POP(7)                                                     00540000
C-------------------------------------                                  00550000
 100  IF(A2MS.LE.0.)A2MS=17.9439                                        00560000
      IF(R0MS.LE.0.)R0MS=1.2249                                         00570000
      IF(C3MS.LE.0)C3MS=0.863987/R0MS                                   00580000
      IF(RKAPPA.LE.0.) RKAPPA=1.7826                                    00590000
      ROT0=51.827/(R0MS*R0MS)                                           00600000
      Y0=ROT0/A2MS                                                      00610000
      X0=C3MS/(2.*A2MS)                                                 00620000
      IF(FINERF.LE.0.01) FINERF=1.                                      00630000
      IF(IFISRT.EQ.2) GOTO 110                                          00640036
C      WRITE(6,105) A2MS,R0MS,C3MS,RKAPPA,ROT0,X0,Y0,FINERF,NOROTF,      00650000
C     1             RKAPA4                                               00660000
C 105  FORMAT(' FISROT PARAMETERS'/' A2MS=',F9.4,4X,'R0MS=',F8.4,4X,     00670028
C     1 'C3MS=',F8.4,4X,'KAPPA=',F8.4,4X,'ROT0,X0,Y0=',F 9.3,2F9.5/      00680028
C     2 'FINERF,NOROTF=',F8.3,I4,4X,'KAPPA4=',F9.4/)                     00690000
 107  FORMAT(18I4)                                                      00700000
      GOTO 115                                                          00710036
 110  WRITE(6,111)                                                      00720036
 111  FORMAT(/' benutzt Sierks Programm BARFIT ')                       00730036
C-------------------------------------- BE'S WITH MS67                  00740029
 115  IZCN=IZT+IZP                                                      00750036
      IACN=MTGT+MPROJ                                                   00760000
      INCN=IACN-IZCN                                                    00770012
      ZCN=IZCN                                                          00780000
      ACN=IACN                                                          00790000
      ZP=IZP                                                            00800000
      AP=MPROJ                                                          00810000
      ZT=IZT                                                            00820000
      AT=MTGT                                                           00830000
      KMAX=KMAX1                                                        00840045
      IF(KMAX.EQ.0)KMAX=3                                               00850000
      DO 120 KK=1,NZMX2                                                 00860000
      DO 120 JJ=1,NAMX2                                                 00870000
 120  EXC(KK,JJ)=0.                                                     00880000
      CALL  MSBEN(ZCN,ACN,NA,NZ,DUM,DUM)                                00890000
      DO 122 KK=1,NZMX2                                                 00900000
      DO 122 JJ=1,NAMX2                                                 00910000
      SMAS(KK,JJ)=0.                                                    00920010
 122  EXC(KK,JJ)=0.                                                     00930000
C-------------------------------------------- USER BE'S --              00940029
      IF(NUMB.LE.0)  GOTO 130                                           00950000
      CALL MYRD(POP,1,23,5)                                             00960041
      QVALUE=POP(1)                                                     00970000
      NUMB1=4*NUMB                                                      00980000
      CALL MYRD(POP,NUMB1,23,5)                                         00990041
      DO 123 L=1,NUMB1,4                                                01000000
      IZ=IZCN+1-POP(L)                                                  01010000
      IA=IACN+2-IZ-POP(L+1)                                             01020000
      K=POP(L+2)+0.01                                                   01030000
      IF(IZ.LE.NZMX .AND. IA.LE.NAMX .AND. K.LE.KMAX)                   01040000
     1             BE(IZ,IA,K)=POP(L+3)                                 01050000
 123  CONTINUE                                                          01060000
 125  FORMAT(3(3I4,F8.3))                                               01070000
      GOTO 25                                                           01080000
C----------------------------------- QVAL WITH MS67 +SHELLS AND PAIRING 01090029
 130  MCC=MC                                                            01100000
      MPP=MP                                                            01110000
      MP=0                                                              01120000
      MC=0                                                              01130000
      CALL MSBEN(ZP,AP,0,0,EXCP ,SHLL)                                  01140000
      CALL MSBEN(ZT,AT,0,0,EXCT ,SHLL)                                  01150000
      CALL MSBEN(ZCN,ACN,0,0,EXCCN,SHLLCN)                              01160029
      MP=MPP                                                            01170000
      MC=MCC                                                            01180000
      QVALUE=EXCP+EXCT-EXCCN                                            01190000
      IF(MASSES.EQ.0)                                                   01200029
     1 WRITE(6,127) MP,MC,ZP,AP,EXCP,ZT,AT,EXCT,ZCN,ACN,EXCCN,QVALUE    01210029
 127  FORMAT(' MSBEN MP,MC=',2I2,3(2F6.0,F9.3,3X),'QVALUE=',F9.3)       01220025
C-------------------------------------------------------------------    01230007
C                 EXP MASSES FROM DISC(EXC),LD MASSES MS67(SMAS)        01240029
      IF(MASSES.EQ.0) GOTO 25                                           01250029
 132  MCC=MC                                                            01260007
      MPP=MP                                                            01270007
      MP=1                                                              01280007
      MC=1                                                              01290007
      DO 135 KK=1,NZMX2                                                 01300007
      DO 135 JJ=1,NAMX2                                                 01310000
 135  EXC(KK,JJ)=0.                                                     01320000
      NA4=NA+4                                                          01330022
      NZ4=NZ+4                                                          01340022
C--------------------------                                             01350007
 140  READ(MASSES,107 ) IZZ,IADN,IAUP                                   01360000
      IF(IZZ.EQ.0) GOTO 152                                             01370016
      NPOW=IAUP-IADN+1                                                  01380000
      READ(MASSES,145 )(POP(I),I=1,NPOW)                                01390000
 145  FORMAT(10F8.3)                                                    01400000
      IF(IZZ.EQ.IZP) EXCP=POP(MPROJ-IADN+1)                             01410000
      IF(IZZ.EQ.IZT) EXCT=POP(MTGT-IADN+1)                              01420000
      IF(IZZ.EQ.IZCN) EXCCNT=POP(IACN-IADN+1)                           01430031
      IZ=IZCN-IZZ+2                                                     01440007
      RZZ=IZZ                                                           01450000
      IF(IZ.GT.NZ4) GOTO 140                                            01460022
      IF(IZ.LT.1) GOTO 152                                              01470017
C--------------------------------                                       01480007
      DO 150  IA=1,NA4                                                  01490022
      IAA=IZZ+INCN+1 + 1-IA                                             01500012
      RAA=IAA                                                           01510012
      CALL MSBEN(RZZ,RAA,0,0,SMAS(IZ,IA),DUM)                           01520012
C                       AMONIX HERE                                     01530039
      I=IAA-IADN+1                                                      01540012
      IF(I.LT.1  .OR.  I.GT.NPOW) GOTO 150                              01550012
      EXC(IZ,IA)=POP(I)                                                 01560000
 150  CONTINUE                                                          01570000
C--------------------------------                                       01580007
      GOTO 140                                                          01590000
C---------------------------------- TEST OUTPUT FOR TESTSE              01600029
 152  IF(JDIM.NE.1) GOTO 160                                            01610021
      WRITE (6,153 )                                                    01620021
 153  FORMAT(10X,'MASS EXCESSES TABLE AND LDM')                         01630016
      NA0=NA+4                                                          01640016
      NA2=MIN0(12,NA4)                                                  01650022
      NA1=MIN0(NA0,13)                                                  01660016
 154  DO 159 IZ = 1,NZ4                                                 01670022
      IIZ=IZCN+1   +1-IZ                                                01680016
      IDUM=IIZ+INCN+1   +1                                              01690016
      DO 155   IA=1,NA4                                                 01700022
 155  IND (IA)=IDUM-IA                                                  01710016
      WRITE(6,290 ) (ELEMNT(IIZ),IND (NA1-IA),IA=1,NA2)                 01720016
      WRITE(6,300 ) (EXC(IZ,NA1-IA),IA=1,NA2)                           01730016
      WRITE(6,300 ) (SMAS(IZ,NA1-IA),IA=1,NA2)                          01740016
 159  CONTINUE                                                          01750016
      IF(NA1.EQ.NA0) GOTO 160                                           01760016
      NA2=MIN0(12,NA0-NA1)                                              01770016
      NA1=NA1+12                                                        01780016
      NA1=MIN0(NA0,NA1)                                                 01790016
      GOTO 154                                                          01800016
C-----------------------------------   GET EXP Q-VALUE                  01810029
 160  MP=MPP                                                            01820007
      MC=MCC                                                            01830007
      DUM=ABS(EXCCNT)                                                   01840031
      ACNLD=EXCCN-SHLLCN                                                01850029
      IF(DUM.GT.0.0001)EXCCN=EXCCNT                                     01860031
      SHLLCN=EXCCN-ACNLD                                                01870032
      IF(MASSES.NE.5) REWIND MASSES                                     01880000
      QVALUE=EXCP+EXCT-EXCCN                                            01890000
      WRITE (6,164) QVALUE                                              01900000
 164  FORMAT(' Q-Wert',F10.3,1X,'MeV')                                   01910000
      IF(ISHELL.EQ.0) GOTO 25                                           01920031
C-----------------------------------   GET EXP BE'S                     01930029
      DO 220  IZ=1,NZ                                                   01940000
      DO 210  IA=1,NA                                                   01950000
      DUM=ABS( EXC(IZ+1,IA+1) )                                         01960007
      IF(DUM.LT.0.0001) GOTO 210                                        01970000
      DUM=ABS(EXC(IZ+1,IA+2))                                           01980007
      IF(DUM.LT.0.0001) GOTO 170                                        01990000
      BE(IZ,IA,1)=8.072 + EXC(IZ+1,IA+2) - EXC(IZ+1,IA+1)               02000007
 170  DUM=ABS(EXC(IZ+2,IA+1))                                           02010007
      IF(DUM.LT.0.0001) GOTO 180                                        02020000
      BE(IZ,IA,2)=7.289 + EXC(IZ+2,IA+1) - EXC(IZ+1,IA+1)               02030007
 180  DUM=ABS(EXC(IZ+3,IA+3))                                           02040007
      IF(DUM.LT.0.0001) GOTO 190                                        02050000
      BE(IZ,IA,3)=2.425 + EXC(IZ+3,IA+3) - EXC(IZ+1,IA+1)               02060007
 190  DUM=ABS(EXC(IZ+2,IA+2))                                           02070007
      IF(DUM.LT.0.0001) GOTO 200                                        02080000
      BE(IZ,IA,4)=13.136 + EXC(IZ+2,IA+2) - EXC(IZ+1,IA+1)              02090007
 200  DUM=ABS(EXC(IZ+2,IA+3))                                           02100007
      IF(DUM.LT.0.0001) GOTO 210                                        02110000
      BE(IZ,IA,5)=14.950 + EXC(IZ+2,IA+3) - EXC(IZ+1,IA+1)              02120007
 210  CONTINUE                                                          02130000
 220  CONTINUE                                                          02140000
 230  WRITE(6,234)                                                      02150000
 234  FORMAT(' Massen aus der Datei "mlz.dat" ')                        02160000
      DO 240  I=1,100                                                   02170000
 240  POP(I)=0.                                                         02180000
C                                                                       02190000
C---------------------------------------------------                    02200000
C                                    GET SHELLS AND DELTAS              02210007
 25   IFLAG=0                                                           02220033
      IF(ISHELL.EQ.0 .AND. IPAIR.EQ.0) IFLAG=1                          02230033
      IF(IFLAG.EQ.0) GOTO 2505                                          02240034
      ISHELL=2                                                          02250034
      IPAIR =4                                                          02260034
      IF( MASSES.EQ.0) ISHELL=1                                         02270034
      IF( MASSES.EQ.0) IPAIR=2                                          02280034
      GOTO 2510                                                         02290033
 2505 IF(MASSES.EQ.0 .AND. ISHELL.EQ.2) ISHELL=1                        02300034
      IF(IPAIR.EQ.4  .AND. ISHELL.NE.2) IPAIR=2                         02310010
 2510 INZ=0                                                             02320033
      INA=0                                                             02330000
      MCC=MC                                                            02340000
      MPP=MP                                                            02350000
      MC=0                                                              02360000
      MP=0                                                              02370000
C-------------------------------                                        02380007
      NZ2=NZ+2                                                          02390022
      DO 40 KZ=1,NZ2                                                    02400022
      KZ1=(KZ-1)*NAMX2                                                  02410000
      IZ=IZCN+1-KZ                                                      02420000
      IZOD=MOD(IZ,2)                                                    02430002
      Z=IZ                                                              02440000
C---------------------                                                  02450007
      NA2=NA+2                                                          02460022
      DO 30 KA=1,NA2                                                    02470022
      KAKZ=KZ1+KA                                                       02480000
      IA=IACN+2-KZ-KA                                                   02490000
      IN=IA-IZ                                                          02500000
      INOD=MOD(IN,2)                                                    02510000
      A=IA                                                              02520000
      CALL MSBEN(Z,A,INA,INZ,CMASS,SHELL)                               02530000
      SHELLS(KAKZ)=SHELL                                                02540000
      SMASS=CMASS-SHELL                                                 02550000
      IF(IPAIR.NE.4) GOTO 26                                            02560014
      ODDEV=-DELT*(1-IZOD-INOD)/SQRT(A)                                 02570014
      SMASS=SMASS-ODDEV                                                 02580038
 26   IF(ISHELL.EQ.2 .AND. EXC(KZ+1,KA+1).NE.0.)SHELLS(KAKZ)=           02590007
     1                     EXC(KZ+1,KA+1)-SMASS                         02600007
      IF(IPAIR.NE.4) GOTO 27                                            02610004
      IDUM=IZOD+INOD                                                    02620004
      IF(IDUM.EQ.KHS+KHS)DELTAS(KAKZ)=0.                                02630037
      DUM1=EXC(KZ,KA+1)                                                 02640007
      DUM2=EXC(KZ+2,KA+1)                                               02650007
      DUM3=EXC(KZ+1,KA)                                                 02660007
      DUM4=EXC(KZ+1,KA+2)                                               02670007
      IF(ABS(DUM1).LT.0.0001) GOTO 2620                                 02680019
      IF(ABS(DUM2).LT.0.0001) GOTO 2620                                 02690019
      IF(ABS(DUM3).LT.0.0001) GOTO 2620                                 02700019
      IF(ABS(DUM4).LT.0.0001) GOTO 2620                                 02710019
      DUM1=DUM1-SMAS(KZ,KA+1)                                           02720007
      DUM2=DUM2-SMAS(KZ+2,KA+1)                                         02730007
      DUM3=DUM3-SMAS(KZ+1,KA)                                           02740007
      DUM4=DUM4-SMAS(KZ+1,KA+2)                                         02750007
      DUM5=0.25*(DUM1+DUM2+DUM3+DUM4)-SHELLS(KAKZ)                      02760021
      IF(IZOD.EQ.KHS) GOTO 2610                                         02770037
      DELTAS(KAKZ)=0.5*(DUM1+DUM2)-SHELLS(KAKZ)                         02780008
      IF(KHS.EQ.0) GOTO 2610                                            02790037
      DO 2605 MAG=1,6                                                   02800037
      IF(IZ.EQ.MAGIC(MAG)) DELTAS(KAKZ)=0.                              02810037
 2605 CONTINUE                                                          02820037
 2610 IF(INOD.EQ.KHS) GOTO 2630                                         02830037
      DELTAN=0.5*(DUM3+DUM4)-SHELLS(KAKZ)                               02840037
      DO 2615 MAG=1,6                                                   02850037
      IF(IN.EQ.MAGICN(MAG)) DELTAN=0.                                   02860037
 2615 CONTINUE                                                          02870037
      DELTAS(KAKZ)=DELTAS(KAKZ)+DELTAN                                  02880037
      GOTO 2630                                                         02890007
 2620 DELTAS(KAKZ)=-ODDEV-(1-2*KHS)*DELT/SQRT(A)                        02900037
C------------------------------------------------------------C          02910037
 2630 SHELLS(KAKZ)=SHELLS(KAKZ)+DELTAS(KAKZ)-KHS*DELT/SQRT(A)           02920038
C-----FOLLOWING FOR OVERALL SHIFT OF REFERENCE SURFACE ------C          02930023
 27   SHELLS(KAKZ)=SHELLS(KAKZ)+SHELL0                                  02940005
 30   CONTINUE                                                          02950000
C---------------------                                                  02960007
 40   CONTINUE                                                          02970000
      IF(IFLAG.EQ.0) GOTO 50                                            02980033
      IF(ABS(SHELLS(1)).GT.0.001) SHLLCN=SHELLS(1)                      02990034
      DO 44 KAKZ=1,KAKZMX                                               03000033
      DELTAS(KAKZ)=0.                                                   03010033
 44   SHELLS(KAKZ)=0.                                                   03020033
      ISHELL=0                                                          03030033
      IPAIR=0                                                           03040033
      IFLAG=0                                                           03050033
C--------------------------------------                                 03060007
 50   MP=MPP                                                            03070033
      MC=MCC                                                            03080000
C------------------------------------------------------------------     03090000
      IF(ISHELL.EQ.1)WRITE(6,60) SHELL0                                 03100033
      IF(ISHELL.EQ.2)WRITE(6,61) SHELL0                                 03110000
 60   FORMAT(/9X,'MS SHELLS  SHELL0=',F8.2)                             03120000
 61   FORMAT(/9X,'Exp. Schaleneffekte,  SHELL0=',F8.2)                  03130000
      NA2=MIN0(12,NA)                                                   03140000
      NA0=NA+1                                                          03150000
      NA1=MIN0(NA0,13)                                                  03160000
 65   DO 90 IZ = 1,NZ                                                   03170010
      IIZ=IZCN+1-IZ                                                     03180000
      IDUM=IACN-IZ+2                                                    03190000
      IZ1=(IZ-1)*NAMX2+NA1                                              03200000
      DO 70    IA=1,NA                                                  03210000
 70   IND (IA)=IDUM-IA                                                  03220000
      WRITE(6,75  ) (ELEMNT(IIZ),IND (NA1-IA),IA=1,NA2)                 03230000
 75   FORMAT(8X,12(2X,A4,I3,1X))                                        03240000
      WRITE(6,85  ) (SHELLS(IZ1-IA) ,IA=1,NA2)                          03250000
 85   FORMAT(8X,12F10.2/)                                               03260000
 90   CONTINUE                                                          03270000
      IF(NA1.EQ.NA0) GOTO 410                                           03280000
      NA2=MIN0(12,NA0-NA1)                                              03290000
      NA1=NA1+12                                                        03300000
      NA1=MIN0(NA0,NA1)                                                 03310000
      GOTO 65                                                           03320000
C----------------------------------- GET F BARRIERS BF,BFLDM----        03330029
 410  AL=0.                                                             03340000
      DUM=0.                                                            03350000
      DO 420 IZ=1,NZ                                                    03360000
      Z=ZCN+1-IZ                                                        03370000
      IZEE=Z+0.01                                                       03380036
      KZ1=(IZ-1)*NAMX2                                                  03390000
      DO 415 IA=1,NA                                                    03400000
      A=ACN+2-IZ-IA                                                     03410000
      MASS=A+0.01                                                       03420036
      KAKZ=KZ1+IA                                                       03430000
      AN=A-Z                                                            03440000
      IF(IFISRT.NE.2) CALL FISROT(A,Z,AN,AL,DELR,DELSP,ERO,DUM)         03450036
      IF(IFISRT.NE.2)  BFLDM(IA,IZ)=DELSP-DELR                          03460036
      IF(IFISRT.EQ.2) CALL BARFIT(IZEE,MASS,0,BFIS,DELR,ELMAX)          03470036
      IF(IFISRT.EQ.2)  BFLDM(IA,IZ)=BFIS                                03480036
      IF(IBF.EQ.0)BF(IA,IZ)=BARFAC*BFLDM(IA,IZ)+BAR0                    03490000
      IF(IBF.EQ.1)BF(IA,IZ)=BARFAC*BFLDM(IA,IZ)+BAR0 - SHELLS(KAKZ)     03500000
 415  CONTINUE                                                          03510000
 420  CONTINUE                                                          03520000
      IF(IBF.LT.2) GOTO 249                                             03530000
C     IF(IBF.EQ.2) READ(5,1070) BF(1,1)                                 03540026
      IF(IBF.EQ.2) CALL MYRD(BF,1,23,5)                                 03550041
 1070 FORMAT(9F8.3)                                                     03560000
      DO 425 IZ=1,NZ                                                    03570000
      DO 424 IA=1,NA                                                    03580000
 424  BF(IA,IZ)=BF(1,1)                                                 03590000
 425  CONTINUE                                                          03600000
C     IF(IBF.EQ.3) READ(5,1070) ((BF(IA,IZ),IA=1,NA),IZ=1,NZ)           03610026
      IRC=NA*NZ                                                         03620026
      IF(IBF.EQ.3) CALL MYRD(BF,IRC,23,5)                               03630041
C---------------------------------------------------------------        03640000
 249  IF(MASSES.NE.0) GOTO 280                                          03650000
 250  WRITE(6,251)                                                      03660000
 251  FORMAT(' MS LYSEKIL 1967')                                        03670000
      IF(MC)260,260,255                                                 03680000
 255  WRITE (6,256 )                                                    03690000
 256  FORMAT(' LIQUID DROP MASSES')                                     03700000
      GO TO 265                                                         03710000
 260  WRITE (6,261 )                                                    03720000
 261  FORMAT(' SHELL CORRECTED MASSES')                                 03730000
 265  IF (MP) 275,275,270                                               03740000
 270  WRITE (6,271 )                                                    03750000
 271  FORMAT(' ZERO PAIRING')                                           03760000
      GO TO 280                                                         03770000
 275  WRITE (6,276 )                                                    03780000
 276  FORMAT(' WITH PAIRING')                                           03790000
C-------------------------------------                                  03800000
 280  WRITE(6,282) BARFAC,BAR0                                          03810000
 282  FORMAT(/10X,' BARFAC=',F8.3,6X,'BAR0=',F8.3)                      03820000
      WRITE (6,281 )                                                    03830000
 281  FORMAT(/10X,'Separationsenergien (BN,BP,BA) und Spaltbarrieren')  03840002
      NA0=NA+1                                                          03850000
      NA2=MIN0(12,NA)                                                   03860000
      NA1=MIN0(NA0,13)                                                  03870000
 283  DO 310 IZ = 1,NZ                                                  03880000
      IIZ=IZCN+1-IZ                                                     03890000
      IDUM=IACN-IZ+2                                                    03900000
      KZ1=(IZ-1)*NAMX2                                                  03910000
      DO 285   IA=1,NA                                                  03920000
 285  IND (IA)=IDUM-IA                                                  03930000
      WRITE(6,290 ) (ELEMNT(IIZ),IND (NA1-IA),IA=1,NA2)                 03940000
 290  FORMAT(8X,12(2X,A4,I3,1X))                                        03950000
      DO 295  K=1,KMAX                                                  03960000
 295  WRITE(6,300 ) (BE(IZ,NA1-IA,K),IA=1,NA2)                          03970000
 300  FORMAT(8X,12F10.2)                                                03980000
      WRITE(6,300)(BF(NA1-IA,IZ),IA=1,NA2)                              03990000
      IF(IPAIR.EQ.4) WRITE(6,300)(DELTAS(KZ1+NA1-IA),IA=1,NA2)          04000000
 310  CONTINUE                                                          04010000
      IF(NA1.EQ.NA0) GOTO 313                                           04020000
      NA2=MIN0(12,NA0-NA1)                                              04030000
      NA1=NA1+12                                                        04040000
      NA1=MIN0(NA0,NA1)                                                 04050000
      GOTO 283                                                          04060000
C----------------------------------------------------                   04070000
 313  QVALUE=EXCP+EXCT-EXCCN                                            04080000
      IF(ISHELL.EQ.0) SHELLS(1)=SHLLCN                                  04090030
      INIT=99                                                           04100027
      RETURN                                                            04110000
      END                                                               04120000
