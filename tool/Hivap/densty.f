      SUBROUTINE  DENSTY (NROW, NSTEP,LDBM,LEP1,LEVPRN)                 00020015
C  
C  Berechnet die Zustandsdichten, wird von HIVAP und EVA aufgerufen
C  Stand: 3.5.1994
C
      INCLUDE 'common.f'                                                00040014
      DIMENSION ARRAY(100),YJK(100),XLDM(9),THMIN(9)                    00050004
C     DIMENSION OF ARRAY MUST BE >= JDIM                                00060000
C     MAYBE UNINITIALIZED: SOR,KI,SPT0,SPT30,HHF,HHG,IAX,KHS,I7
C
C     I.Giese 2007
C     WICHTIG AUFZUHEBEN: RJLIM,RJDIF,ARRAY wegen NAN sonst
C
      DATA XLDM,THMIN/.3 ,.4  , .5, .6 ,.67, .7 , .8, .9 , 1.,          00070004
     1                .57,.555,.55,.495,.47,.485,.62,.795,1.0/          00080005
C--------------------------------------------------------------------   00090014
C
C     I.Giese 2007
      SAVE
C
C     MINIMUM MOMENT OF INERTIA FROM CPS74                              00100009
      AMPAR1=0.0
      AMPAR2=30.
      AMPAR3=2.0
      AMPAR4=0.3
      AMPAR5=10.  
      KJC1=4*JDIM +1                                                    00110000
      PINERF=(YRSMO(KJC1+30)-YRSMO(KJC1))/(31.*30.)                     00120002
      PINERT=(YRSMO(31)-YRSMO(1))/(31.*30.)                             00130002
      IF (PINERF.LT.0. .OR. PINERT.LT.0.)WRITE(6,210)YRSMO(1),YRSMO(2), 00140001
     1YRSMO(KJC1),YRSMO(KJC1+1),NROW,NSTEP                              00150000
 210  FORMAT(' DENSTY YRSMO ERROR',4E12.4,2I6)                          00160000
C     ABOVE IS HBAR-SQUARE OVER INERTIA FOR SADDLE AND GS               00170000
      BETAF=(PINERT/PINERF-1.)/0.28                                     00180001
      LEP=LEP1                                                          00190015
      LEP=MAX0(1,LEP)                                                   00200000
      IDUM=0                                                            00210000
      DO 10 K=1,5                                                       00220000
 10   IDUM=IDUM+LADJ(K)                                                 00230000
      IF(IDUM.NE.0)  GO TO 316                                          00240000
C      LADJ IS USED IN ENTRANCE FROM EVAP, TO CONTROL THE TABLE THAT    00250000
C           IS BEING ADJUSTED                                           00260000
C---------------------------------------------------------------------  00270000
C                                        RESET                          00280000
      SOR=0.                                                            00290000
      IF(EXMAX1.LT.200.) GOTO 15                                        00300000
      SOR=SQRT(ACN*EXMAX1*0.1)                                          00310000
  15  DO 20  K=1,IDEN                                                   00320000
      OMEGF(K)=0.                                                       00330000
      OMEGG(K)=0.                                                       00340000
      OMEGN(K)=0.                                                       00350000
  20  OMEGP(K)=0.                                                       00360000
      DO 30  K=1,IDNA                                                   00370000
  30  OMEGA(K)=0.                                                       00380000
      DO 40 K=1,IEDIM                                                   00390000
      ENHANC(K)=0.                                                      00400000
  40  ENHANS(K)=0.                                                      00410000
C---------------------------------------------------------------------- 00420000
C                  SET DIMENSIONS         EXMAX(N,P,A,G,F)              00430000
      EXMAXG=EXMAX1-DEL(1)                                              00440000
      EXMAXN=EXMAX1+Q(2)-DEL(2)                                         00450000
      EXMAXP=EXMAX1+Q(3)-DEL(3)                                         00460000
      EXMAXA=EXMAX1+Q(4)-DEL(4)                                         00470000
      EXMAXF=EXMAX1-DEL(5)                                              00480000
C---------------------------------------------------------------------- 00490000
C                                           NW(G,N,P,A,F)               00500000
      IT=1                                                              00510000
      EXMAX=EXMAXG                                                      00520000
      KLAST=MIN0(N1+2,JDIM)                                             00530000
      KMAX=KLAST+1                                                      00540000
    4 DO 302 K=1,KLAST                                                  00550000
      KI=KMAX-K                                                         00560000
      IF(EJAY(KI).LE.EXMAX) GO TO 303                                   00570000
  302 CONTINUE                                                          00580000
  303 GO TO (5,6,7,8,9),IT                                              00590000
    5 IT=2                                                              00600000
      EXMAX=EXMAXN                                                      00610000
      KLAST=MIN0(N1+LN-1,JDIM)                                          00620000
      KMAX=KLAST+JDIM+1                                                 00630000
      NWG=KI                                                            00640000
      GO TO 4                                                           00650000
    6 IT=3                                                              00660000
      EXMAX=EXMAXP                                                      00670000
      KLAST=MIN0(N1+LP-1,JDIM)                                          00680000
      KMAX=KLAST+2*JDIM+1                                               00690000
      NWN=KI-JDIM                                                       00700000
      GO TO 4                                                           00710000
    7 IT=4                                                              00720000
      EXMAX=EXMAXA                                                      00730000
      KLAST=MIN0(N1+LA-1,JDIM)                                          00740000
      KMAX=KLAST+3*JDIM+1                                               00750000
      NWP=KI-2*JDIM                                                     00760000
      GO TO 4                                                           00770000
    8 IT=5                                                              00780000
      EXMAX=EXMAXF                                                      00790000
      KLAST=MIN0(N1,JDIM)                                               00800000
      KMAX=KLAST+4*JDIM+1                                               00810000
      NWA=KI-3*JDIM                                                     00820000
      GOTO 4                                                            00830000
 9    NWF=KI-4*JDIM                                                     00840000
C-----------------------------------------------------------------------00850000
C                                              EXMIN(G,N,P,A,F)         00860000
C                                         DECIMAL PARTS OF EXMAX        00870000
      EXMING=EXMAXG-AINT(EXMAXG)                                        00880000
      EXMINN=EXMAXN-AINT(EXMAXN)                                        00890000
      EXMINP=EXMAXP-AINT(EXMAXP)                                        00900000
      EXMINA=EXMAXA-AINT(EXMAXA)                                        00910000
      EXMINF=EXMAXF-AINT(EXMAXF)                                        00920000
C-----------------------------------------------------------------------00930000
C                                              MW(G,N,P,A,F)            00940000
      MWG=EXMAXG-EXMING+1.0001                                          00950000
      MWG=MIN0(MWG,IDEN/NWG)                                            00960000
      MWG=MAX0(MWG,1)                                                   00970000
      MWN=EXMAXN-EXMINN+1.0001                                          00980000
      MWN=MIN0(MWN,IDEN/NWN)                                            00990000
      MWN=MAX0(MWN,1)                                                   01000000
      MWP=EXMAXP-EXMINP+1.0001                                          01010000
      MWP=MIN0(MWP,IDEN/NWP)                                            01020000
      MWP=MAX0(MWP,1)                                                   01030000
      MWA=EXMAXA-EXMINA+1.0001                                          01040000
      MWA=MIN0(MWA,IDNA/NWA)                                            01050000
      MWA=MAX0(MWA,1)                                                   01060000
      MWF=EXMAXF-EXMINF+1.0001                                          01070000
      MWF=MIN0(MWF,IDEN/NWF)                                            01080000
      MWF=MAX0(MWF,1)                                                   01090000
C-----------------------------------------------------------------------01100000
C                                   AM DEPENDENCE AF/AN AND SHELLS      01110000
      AMASS=MTGT+MPROJ+2-NSTEP-NROW                                     01120000
      ZEE=IZT+IZP+1-NROW                                                01130000
      AN=AMASS-ZEE                                                      01140000
      JSHAP=97                                                          01150000
      JBAR=99                                                           01160000
      DLIM=120/AMASS                                                    01170000
      KJC=4*JDIM                                                        01180000
      DO 3040 KJ=1,JDIM                                                 01190000
      J=KJC+KJ                                                          01200000
      AA=KJ-1                                                           01210000
C
C     I.Giese 2007
C     0.0D0
C
      CALL FISROT(AMASS,ZEE,AN,AA,DELR,SPT,EROT,0.0D0)                  01220000
      IF(KJ.EQ.1)SPT0=SPT                                               01230009
      IF(KJ.EQ.31)SPT30=SPT                                             01240009
      IF(SPT0.LT.5.) GOTO 3030                                          01250000
      IF(KJ.EQ.2) HHF=SPT-SPT0                                          01260000
      IF(KJ.EQ.2) HHG=DELR                                              01270000
      SPT1=AA*(AA+1)*HHF+SPT0 -SPT                                      01280000
      DELR1=AA*(AA+1)*HHG -DELR                                         01290000
      IF((SPT1.GT.DLIM .OR. DELR1.GT.DLIM) .AND. JSHAP.EQ.97)JSHAP=KJ-1 01300000
      IF(SPT-DELR .LT. DLIM .AND.JBAR.EQ.99) JBAR=KJ-1                  01310000
 3030 DUM=(SPT-DELR)/SPT0                                               01320000
      DUM=DMAX1(AMPAR1,DUM)                                             01330000
 3040 ARRAY(KJ)=DUM                                                     01340000
      PINERF=(SPT30-SPT0)/(30.*30.)                                     01350009
      RJLIM=0.5*(JBAR+JSHAP)                                            01360000
      RJLIM=DMAX1(AMPAR2,RJLIM)                                         01370000
      RJDIF=0.2*(JBAR-JSHAP)                                            01380000
      RJDIF=DMAX1(AMPAR3,RJDIF)                                         01390000
      RJDIF=DMIN1(AMPAR5,RJDIF)                                         01400000
      IF(LPRINT.LT.3) WRITE(6,3041) JSHAP,JBAR,RJLIM,RJDIF              01410000
 3041 FORMAT(' DENSTY: JSHAP=',I3,4X,'JBAR=',I3,4X,2F8.1)               01420000
C-----------------------------------------------------------------------01430000
C                                   GAMMA PARAMETERS                    01440000
      LE=1                                                              01450000
      IT=1                                                              01460000
      IT1=IT                                                            01470000
  304 CR=R(1)                                                           01480000
      CAL=AL(1)                                                         01490000
      MA=KZ(NSTEP)+IN(NSTEP)                                            01500000
      IZ=KZ(NSTEP)                                                      01510000
      IMAX=MWG                                                          01520000
      KE0=0                                                             01530000
      NJ=NWG                                                            01540000
      DELTAC=DELTA(1)                                                   01550000
      SHELL=SHELK(1)                                                    01560000
      EXMAX=EXMAXG                                                      01570000
      KNDEX=0                                                           01580000
      KJC=0                                                             01590000
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)                     01600002
      PINERT=PINERT**1.5 /24.                                           01610000
      IAX=0                                                             01620006
C     KNDEX IS OMEG OFFSET,KJC IS EJAY OFFSET                           01630000
C***********************************************************************01640000
C                              CALCULATE LEVEL DENSITY                  01650000
  305 RAL=CR*CAL                                                        01660000
      IF(EXMAX.LE.0.)GOTO 3111                                          01670000
      IF(CR.LE.0.)      RAL=1.                                          01680000
      SC=CAL                                                            01690000
C-----------------------------------------------------------------------01700000
C                                  ENHANCEMENT FUNCTION                 01710000
C                 USE ONLY FOR EXCIT < 100 MEV AND LOW ANG.MOMENTA      01720000
C                 USES SAME FOR GAMMA AND N CHANNELS                    01730000
      IF(IENH.EQ.0) GOTO 360                                            01740000
      IF(LE.NE.1) GOTO 360                                              01750000
      KEMAX=EXMAX+5.                                                    01760000
      KEMAX=MIN0(100,KEMAX)                                             01770000
      DO 350 KE=1,KEMAX                                                 01780000
      E=KE                                                              01790000
      KKE=KE+KE0                                                        01800000
      IF(IT1.NE.5) GOTO 345                                             01810000
      BETAMX=BETAF                                                      01820000
      UEFF=E                                                            01830000
      IF(ESHELL.NE.0.) UEFF=E+SHELL*(1.-EXP(-E/ESHELL))                 01840000
      GOTO 348                                                          01850000
 345  CALL BETA(IZ,MA,SHELL,DUM,E,ESHELL,BETAMX,FROT,UEFF)              01860000
 348  IF(UEFF.LE.0.) GOTO 350                                           01870000
      TEMP=SQRT(UEFF/CAL)                                               01880000
      SIG2=TEMP/PINERF                                                  01890001
      IF(IT1.NE.5) SIG2=(1.+0.28*BETAMX)*TEMP/PINERT                    01900001
      TCROT=37.8*BETAMX/FLOAT(MA)**0.33333                              01910000
      EDROT=CAL*TCROT*TCROT                                             01920000
      IF(EDROT.LE.0.)GOTO 350                                           01930000
      DUM=UEFF/EDROT                                                    01940000
      IF(EDCOLL.GT.0.) DUM=DUM+UEFF/EDCOLL                              01950000
      ENHANC(KKE)=0.                                                    01960000
      IF(DUM.GT.150.) GOTO 350                                          01970000
      ENHANC(KKE)=SIG2*EXP(-DUM)                                        01980000
      ENHANS(KKE)=UEFF                                                  01990000
 350  CONTINUE                                                          02000000
C     IF(IENHP.EQ.0 ) GOTO 360                                          02010000
C     WRITE(6,352) IT1,SHELL,BETA0,BETAF                                02020000
 352  FORMAT(/' DENSTY  ENHANC-FACTOR  IT,SHELL,BETA0=',I3,F8.2,F8.3,   02030000
     1       2X,'BETAF=',F8.3/)                                         02040000
C     WRITE(6,354) (ENHANC(KE+KE0),KE=1,KEMAX)                          02050000
 354  FORMAT(10F8.2)                                                    02060000
C-----------------------------------------------------------------------02070000
C                                        LEVDENS AS IN BM69             02080000
C--------------------------------------------------                     02090007
 360  DO 490 KJ=1,NJ                                                    02100007
      RJ=KJ-1                                                           02110000
      IROW=(KJ-1)*IMAX                                                  02120000
      J=KJC+KJ                                                          02130000
      YJ=EJAY(J)                                                        02140000
      IF(LDBM.EQ.3 .OR. LDBM.EQ.4) YJ=YRSMO(J)                          02150000
      PREEX=1.                                                          02160007
      IF(IAX.EQ.0) PREEX=KJ+KJ-1                                        02170000
      IF(IPREEX.EQ.1) PREEX= PREEX*PINERT                               02180000
      SHELJ=YJ-YRSMO(J)+SHELL                                           02190000
 405  IF(IT.NE.5) GOTO 410                                              02200000
C      IT=5 IS FISSION CHANNEL                                          02210000
      SSS=SHELK(1)-SHELK(5)                                             02220000
      IF(ISHELL.EQ.0) SSS=0.                                            02230000
      IF(LDBM.EQ.3 .OR. LDBM.EQ.4) YJ=YRSMO(J)-SSS                      02240000
      ALJF=(AL(5)-AL(1))*(1./(1.+EXP((RJ-RJLIM)/RJDIF) ))               02250000
C
C     I.Giese 2007
C     PC/Linux: RJLIM,RJDIF haben ihre Werte nicht behalten
C
C     WIRKUNG: EXP(0.0/0.0) => NAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF(LDBM.EQ.2 .OR.LDBM.EQ.4) ALJF=AL(5)-AL(1)                      02260000
      CAL=AL(1)+ALJF                                                    02270000
      SHELJ=SHELL*ARRAY(KJ)                                             02280000
C
C     I.Giese 2007
C     PC/Linux: ARRAY-Elemente haben ihre Werte nicht behalten
C
C     WIRKUNG: ARRAY(KJ) = NAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
 410  DUM=1.05*ABS(SHELJ)                                               02290000
      ESHELJ=DMAX1(ESHELL,DUM)                                          02300000
      KUP=1                                                             02310000
      IF(IAX.EQ.1) KUP=KJ                                               02320013
      IU=1                                                              02330000
C--------------------------------------------                           02340000
      DO 470 KE=LE,IMAX                                                 02350000
      INDEX=KNDEX+IROW+KE                                               02360000
      E=EXMAX+1.-FLOAT(KE)                                              02370000
      OM=0.                                                             02380000
      PREEX1=PREEX                                                      02390013
      DO 455 KVAL=1,KUP                                                 02400000
      U=E-YJ                                                            02410009
      IF(IAX.EQ.1) U=U-YJK(KVAL)                                        02420009
      IF(U.LE.0.) GOTO 455                                              02430010
      HE=1.                                                             02440007
      IF(U.LT.UCRIT .AND. KHS.EQ.1) HE=(1.-U*U/(UCRIT*UCRIT))           02450007
      U=U -DELTAC*HE                                                    02460010
      IF(U.LE.0.) GOTO 455                                              02470010
      IF(IENH.EQ.0) GOTO 450                                            02480000
      IU=U+0.5                                                          02490000
      IU=MAX0(1,IU)                                                     02500000
      IU=MIN0(100,IU)                                                   02510000
      RU=0.                                                             02520000
      IF(U.LT.100.) RU=U-IU                                             02530000
      U=ENHANS(IU+KE0)+RU                                               02540000
      IF(U.LE.0.) GOTO 455                                              02550010
 450  CALS=CAL                                                          02560000
      U1=DMAX1(AMPAR4,U)                                                02570000
      FF=1.                                                             02580000
      IF(IENH.EQ.0) FF=1.+SHELJ*(1.-EXP(-U1/ESHELJ))/U1                 02590000
      IF(ISHELL.NE.0)CALS=CAL*FF                                        02600000
      SQ=SQRT(CALS*U)                                                   02610000
      SQ=SQ+SQ-SOR                                                      02620000
      OM= OM+       SQRT(CALS)*PREEX1*EXP(SQ)/(U1*U1)                   02630013
      IF(KVAL.EQ.2) PREEX1=2.*PREEX                                     02640013
 455  CONTINUE                                                          02650010
C
C     I.Giese 2007
C     Nicht lesbar: INDEX > 2500, bis zu 15000 !!!!!!!!!!!!!!!!!!!!!!!!
C common.f: OMEGG(2500),OMEGN(2500),OMEGP(2500),OMEGA(5000),OMEGF(2500)
C
      OMEGG(INDEX) =OM * (1.+ENHANC(IU+KE0))                            02660000
C     GOTO 470                                                          02670011
C460  OMEGG(INDEX)=0.                                                   02680011
 470  CONTINUE                                                          02690000
C---------------------------------------------                          02700000
 490  CONTINUE                                                          02710000
C                                                                       02720000
      IF(LEVPRN.NE.1) GOTO 3111                                         02730000
      WRITE(6,493) MA,KZ(NSTEP),EXMAX,CAL,UCRIT,ESHELL,SHELL,           02740000
     1             DELTAC                                               02750000
 493  FORMAT(' DENSTY MASS,Z',2I4/' EXMAX AL UCRIT ESHELL',4F10.3/      02760000
     1         ' SHELL,DELTA',2F12.5,2F10.3/)                           02770000
      DO 495 KJ=1,NJ                                                    02780000
      IROW=(KJ-1)*IMAX                                                  02790000
      WRITE(6,491) KJ,IT,LE,IMAX,NJ                                     02800000
 491  FORMAT(' DENSTY KJ,IT,LE,IMAX,NJ',5I6)                            02810000
      DO 492 KE=LE,IMAX                                                 02820000
      POPA(KE)=0.                                                       02830000
      INDEX=KNDEX+IROW+KE                                               02840000
      IF(OMEGG(INDEX).LE.0.) GOTO 492                                   02850000
      POPA(KE)=DLOG (OMEGG(INDEX))                                      02860000
 492  CONTINUE                                                          02870000
      WRITE(6,496) (POPA(KE),KE=LE,IMAX)                                02880000
 495  CONTINUE                                                          02890000
 496  FORMAT(10F10.3)                                                   02900000
      DO 497 KE=LE,IMAX                                                 02910000
 497  POPA(KE)=0.                                                       02920000
C                                                                       02930000
C-----------------------------------------------------------------------02940000
C                                                  SWITCH               02950000
 3111 GO TO (312,313,314,515,315,320,322,324,330),IT                    02960000
C-----------------------------------------------------------------------02970000
C                                     NEUTRON PARAMETERS                02980000
  312 IT=2                                                              02990000
      IT1=IT                                                            03000000
 3121 CR=R(2)                                                           03010000
      CAL=AL(2)                                                         03020000
      MA=KZ(NSTEP)+IN(NSTEP)-1                                          03030000
      IZ=KZ(NSTEP)                                                      03040000
      DELTAC=DELTA(2)                                                   03050000
      SHELL=SHELK(2)                                                    03060000
      IMAX=MWN                                                          03070000
      KE0=0                                                             03080000
      NJ=NWN                                                            03090000
      EXMAX=EXMAXN                                                      03100000
      KNDEX=IDEN                                                        03110000
      KJC=JDIM                                                          03120000
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)                     03130002
      PINERT=PINERT**1.5 /24.                                           03140000
      IAX=0                                                             03150006
      GO TO 305                                                         03160000
C-----------------------------------------------------------------------03170000
C                                       PROTON PARAMETERS               03180000
  313 IT=3                                                              03190000
      IT1=IT                                                            03200000
 3131 CR=R(3)                                                           03210000
      IU0=0                                                             03220000
      CAL=AL(3)                                                         03230000
      MA=KZ(NSTEP)+IN(NSTEP)-1                                          03240000
      IZ=KZ(NSTEP)-1                                                    03250000
      DELTAC=DELTA(3)                                                   03260000
      SHELL=SHELK(3)                                                    03270000
      IMAX=MWP                                                          03280000
      KE0=100                                                           03290000
      NJ=NWP                                                            03300000
      EXMAX=EXMAXP                                                      03310000
      KNDEX=2*IDEN                                                      03320000
      KJC=2*JDIM                                                        03330000
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)                     03340002
      PINERT=PINERT**1.5 /24.                                           03350000
      GO TO 305                                                         03360000
C-----------------------------------------------------------------------03370000
C                                      ALPHA  PARAMETERS                03380000
  314 IT=4                                                              03390000
      IT1=IT                                                            03400000
 3141 CR=R(4)                                                           03410000
      CAL=AL(4)                                                         03420000
      MA=KZ(NSTEP)+IN(NSTEP) -4                                         03430000
      IZ=KZ(NSTEP)-2                                                    03440000
      DELTAC=DELTA(4)                                                   03450000
      SHELL=SHELK(4)                                                    03460000
      IMAX=MWA                                                          03470000
      KE0=200                                                           03480000
      NJ=NWA                                                            03490000
      EXMAX=EXMAXA                                                      03500000
      KNDEX=3*IDEN                                                      03510000
      KJC=3*JDIM                                                        03520000
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)                     03530002
      PINERT=PINERT**1.5 /24.                                           03540000
      IAX=0                                                             03550006
      GO TO 305                                                         03560000
C-----------------------------------------------------------------------03570000
C                                    FISSION  PARAMETERS                03580000
  515 IT=5                                                              03590000
      IT1=IT                                                            03600000
 5151 CR=R(5)                                                           03610000
      CAL=AL(5)                                                         03620000
      MA=KZ(NSTEP)+IN(NSTEP)                                            03630000
      IZ=KZ(NSTEP)                                                      03640000
      DELTAC=DELTA(5)                                                   03650000
      SHELL=SHELK(5)                                                    03660000
      IMAX=MWF                                                          03670000
      KE0=300                                                           03680000
      NJ=NWF                                                            03690000
      EXMAX=EXMAXF                                                      03700000
      KNDEX=3*IDEN+IDNA                                                 03710000
      KJC=4*JDIM                                                        03720000
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)                     03730012
      IF(AX.LE.0.) GOTO 520                                             03740008
C     APPROXIMATE DETERMINATION OF INERTIA PARAMETER PARALLEL TO K-AXIS 03750000
      IAX=1                                                             03760006
      XLD=FLOAT(IZ*IZ)/FLOAT(MA)                                        03770006
      DUM=(MA-2*IZ)/FLOAT(MA)                                           03780006
      DUM=50.883*(1.-1.7826*DUM*DUM)                                    03790006
      XLD=XLD/DUM                                                       03800006
      PINERK=TABIP(XLDM,THMIN,9,2,0,XLD,0)                              03810010
      DUM=(YRSMO(31)-YRSMO(1))/(31.*30.)                                03820002
      THMN=PINERK                                                       03830009
      PINERK=DUM/PINERK                                                 03840006
      PINER0=PINERK-PINERT                                              03850000
      IF(I7.NE.99)                                                      03860010
     1WRITE(6,517) PINERT,PINERK,DUM,XLD,THMN                           03870010
 517  FORMAT(' DENSTY  INERF, INERK, INERGS',3F9.4,4X,'XLD THMIN',2F9.4)03880009
      I7=99                                                             03890010
      DO 518 KVAL=1,NJ                                                  03900000
      DUM=KVAL-1                                                        03910000
 518  YJK(KVAL)= PINER0*(DUM+1.)*DUM                                    03920009
C                                                                       03930000
      PINERT=PINERT*SQRT(PINERK)/24.                                    03940013
      GOTO 305                                                          03950000
 520  PINERT=PINERT**1.5/24.                                            03960000
      GO TO 305                                                         03970000
C***********************************************************************03980000
C                                            EXIT                       03990000
  315 CONTINUE                                                          04000000
      LEVPRN=0                                                          04010000
      RETURN                                                            04020000
C***********************************************************************04030000
C                                     ENTRY POINT FROM EVA              04040000
C-----------------------------------------------------------------------04050000
C                                         READJUST GAMMAS               04060000
  316 IF(LADJ(1).NE.1)     GO TO 320                                    04070000
      IT=6                                                              04080000
      IT1=1                                                             04090000
      LE=MWG-LEP                                                        04100000
      NJ=NWG                                                            04110000
      IMAX=MWG                                                          04120000
      KNDEX=0                                                           04130000
C                                                                       04140000
  317 DO 318 KJ=1,NJ                                                    04150000
      INDEX=KNDEX+(KJ-1)*IMAX                                           04160000
      DO 318 K=1,LE                                                     04170000
      INDEX=INDEX+1                                                     04180000
      ITM=INDEX+LEP                                                     04190000
C     I.Giese 2007: INDEX > 2500, bis zu 15000 !!!!!!!!!!!!!!!!!!!!!!!!
      OMEGG(INDEX)=OMEGG(ITM)                                           04200000
  318 OMEGG(ITM)=0.                                                     04210000
      GO TO (319,321,323,325, 335),IT1                                  04220000
C--------------------------------------------                           04230000
  319 LE=MWG+1-LEP                                                      04240000
      LADJ(1)=0                                                         04250000
      EXMAXG=EXMAXG-FLOAT(LEP)                                          04260000
      GO TO 304                                                         04270000
C-----------------------------------------------------------------------04280000
C                                          READJUST NEUTRONS            04290000
  320 IF(LADJ(2).NE.1)      GO TO 322                                   04300000
      IT=7                                                              04310000
      IT1=2                                                             04320000
      LE=MWN-LEP                                                        04330000
      NJ=NWN                                                            04340000
      IMAX=MWN                                                          04350000
      KNDEX=IDEN                                                        04360000
      GO TO 317                                                         04370000
C-------------------------------------------                            04380000
  321 LE=MWN+1-LEP                                                      04390000
      LADJ(2)=0                                                         04400000
      EXMAXN=EXMAXN-FLOAT(LEP)                                          04410000
      GO TO 3121                                                        04420000
C-----------------------------------------------------------------------04430000
C                                         READJUST PROTONS              04440000
  322 IF(LADJ(3).NE.1)     GO TO 324                                    04450000
      IT=8                                                              04460000
      IT1=3                                                             04470000
      LE=MWP-LEP                                                        04480000
      NJ=NWP                                                            04490000
      IMAX=MWP                                                          04500000
      KNDEX=2*IDEN                                                      04510000
      GO TO 317                                                         04520000
C--------------------------------------------                           04530000
  323 LE=MWP+1-LEP                                                      04540000
      LADJ(3)=0                                                         04550000
      EXMAXP=EXMAXP-FLOAT(LEP)                                          04560000
      GO TO 3131                                                        04570000
C-----------------------------------------------------------------------04580000
C                                         READJUST ALPHAS               04590000
  324 IF(LADJ(4).NE.1) GOTO 330                                         04600000
      IT=9                                                              04610000
      IT1=4                                                             04620000
      LE=MWA-LEP                                                        04630000
      NJ=NWA                                                            04640000
      IMAX=MWA                                                          04650000
      KNDEX=3*IDEN                                                      04660000
      GO TO 317                                                         04670000
C------------------------------------------                             04680000
  325 LE=MWA+1-LEP                                                      04690000
      LADJ(4)=0                                                         04700000
      EXMAXA=EXMAXA-FLOAT(LEP)                                          04710000
      GO TO 3141                                                        04720000
C---------------------------------------------------------------------- 04730000
C                                         READJUST FISSION              04740000
  330 IF(LADJ(5).NE.1) GOTO 315                                         04750000
      IT1=5                                                             04760000
      LE=MWF-LEP                                                        04770000
      NJ=NWF                                                            04780000
      IMAX=MWF                                                          04790000
      KNDEX=3*IDEN+IDNA                                                 04800000
      GO TO 317                                                         04810000
C------------------------------------------                             04820000
  335 LE=MWF+1-LEP                                                      04830000
      LADJ(5)=0                                                         04840000
      EXMAXF=EXMAXF-FLOAT(LEP)                                          04850000
      GO TO 515                                                         04860000
      END                                                               04870000
