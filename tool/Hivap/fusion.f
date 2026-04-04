      SUBROUTINE FUSION(ELAB,AP,ZP,AT,ZT,QQ,V0,R00,D,XTH,APUSH,FPUSH,   00020077
     1                   CRED,NORUTH,                                   00030071
     2             NOPROX,NOCURV,ION,SIGML,FLGRAZ,LMAX,SIGF,ITEST,IOPT) 00040065
C
C  Berechnet Fusionsquerschnitte nach der WKB-Methode fuer verschiedene
C  Kern-Kern-Potentiale, erlaubt deformierte Target-Kerne und 
C  erhoehte Sub-Barrieren - Transparenz
C  ruft VLRO(Eingaenge V0, VL, VR), TWKB(Eingang WKB), PARABO, FISROT
C  auf; maechtigere Version FUSIO wurde entwickelt
C  Modifikation 8.10.1987
C  Modifikation 7.6.1994 Drehimpulsbegrenzung von BAR>0.5 auf 
C  BAR>0.1 gesenkt, da vorher LCRIT = 1 fuer Z_CN > 105 
C  Stand: 7.6.1994
C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      COMMON/FUS/VBFUS,RBFUS,HWFUS                                      00050072
C--------------------------------------------------------------------   00060057
C     VERSION OF FUSION ALLOWING TO COMMUNICATE BARRIER PARAMETERS      00070073
C--------------------------------------------------------------------   00080073
C     24 PARAMETRS                                                      00090065
C     AP PROJECTILE MASS (SPHERIC)                                      00100066
C     AT TARGET MASS (DEFORMED)                                         00110066
C     ZP PROJECTILE Z  (SPHERIC)                                        00120066
C     ZT TARGET(DEFORMED) Z                                             00130066
C     QQ  ALIAS Q2 :QUADRUPOLE MOMENT OF TARGET IN (FM**2)              00140067
C     V0,R00,D : POT PARAMETRES                                         00150079
C     FOR MASSES 1,2,3: RPS+RTS=R00*AT3                                 00160079
C     FOR MASS   4    : RPS+RTS=R00*AT3+1.6                             00170079
C     FOR MASSES >4   : RPS+RTS=R00*(AT3 + AP3)                         00180079
C     XTH,APUSH,FPUSH:PARAMETERS,XTH,A,F OF BJORNHOM-SWIATECKI 82 EXTRAP00190067
C     CRED : ALLOWS TO MANIPULATE BARRIER TRANSPARENCY,USUALLY=1        00200067
C            IN PARABOLIC APPROX CRED IS PROPORT TO 1./SQRT(OMEGA)      00210067
C     NORUTH =1:DOES NOT USE CHANGE OF ORIENTATION ANGLE ALONG TRAJ.    00220066
C     NORUTH =2:SAME AS 1,BUT ALSO DOES NOT USE QUADRATIC COULOMB INSIDE00230066
C     NOPROX =1:USES CENTER LINE PRESCRIPTION(INSTEAD OF PROX.PRESCR.)  00240066
C     NOCURV =1:IGNORES CURVATURE DEPENDANCE OF POT                     00250066
C     LMAX   INPUT:  -1   MAXIMUM L WILL BE DETERMINED FROM             00260057
C                         BFLDM(L)<0.5 MEV (NEEDS FISROT)               00270060
C                    <-1  MAXIMUM L CALCULATED IS |LMAX|                00280057
C                    >=0  MAXIMUM L IS DETERMINED INTERNALLY            00290057
C             OUTPUT      MAXIMUM L CALCULATED                          00300057
C     SIGF    OUTPUT: FUSION CROSS SECTION IN MB                        00310067
C     ITEST :IS SET 99 IF POT HAS NO POCKET FOR L=0                     00320067
C            2,3,4,5    BIG PRINT OUT                                   00330067
C            1          SMALL PRINT OUT                                 00340067
C            <1  OR >5  NOPRINT OUT                                     00350067
C            3 OR -1    PARABOLIC APPROX.(RATHER THAN WKB)              00360067
C            -2         CLASSICAL APPROX.(NO TUNNELLING)                00370067
C     IOPT =-2:SAME AS -1 BUT FOR AP=1 SPIN ORBIT L+1/2 IS ADDED        00380067
C     IOPT =-1:NUCLEAR IS WS,NO RADII CORRECTION,USE CENTRAL RADII      00390067
C          =0:VNUCLEAR IS EXPONENTIAL,NO RADII CORRECTION,CENTRAL RADII 00400068
C          =1:VNUCLEAR IS EXPONENTIAL,WITH RADII CORRECTION,SHARP RADII 00410068
C          =2:FULL PROXIMITY POT,WITH RADII CORRECTION,USE SHARP RADII  00420068
C          =3:FULL PROXIMITY POT,NO RADII CORRECTION,USE CENTRAL RADII  00430068
C          =4: EXP PROXIMITY POT,WITH RADII CORRECTION,USE SHARP RADII  00440068
C          =5: BASS POTENTIAL,BERLIN 79, SET V0,R0,D=0.                 00450068
C     ION  =1: 5 POINT ORIENTATIONAL INTEGRATION                        00460068
C          =2:10 POINT ORIENTATIONAL INTEGRATION                        00470068
C          =3:OCOS=0  OR THETA=90 DEGREE                                00480068
C          =4:OCOS=1  OR THETA= 0 DEGREE                                00490068
C------------------------------------------------------------------     00500057
      DIMENSION Z(17),WEIGHT(17),X(3),Y(3),A1(3),SIGML(1),RMX(17),      00510000
     1 T0(17),VMX(17),TMX(17),HWMX(17),IWKB(17),RINN(17),ROUTT(17)      00520011
     2 ,TPARAP(17),SIGOR(17)                                            00530045
      COMMON/FISRT/ROT0,A2MS,RKAPPA,X0,Y0                               00540030
C     MAYBE UNINITIALIZED: TCL,JWKB
      DATA WEIGHT/.23693,.47863,.56889,.47863,.23693,.06667,.14945,     00550000
     1 .21909,.26927,.29552,.29552,.26927,.21909,.14945,.06667,2.0,2.0/,00560000
     2 Z/.04691,.23077,.50000,.76923,.95309,.01305,.06747,.16030,.28330,00570000
     3 .42556,.57444,.71670,.83397,.93253,.98693,0.,1./,                00580000
     4 HBARC,AMU/197.329,931.502/                                       00590079
C
C     I.Giese 2007
      SAVE
C
      AMPAR1=0.2
      AMPAR2=0.02
      LUP=200                                                           00600005
      AMXEXP=100.                                                       00610055
      TMIN=EXP(-AMXEXP)                                                 00620055
C-----------------------------------------S1--------------------------  00630000
      R0FUS=0.                                                          00640069
      Q2=QQ                                                             00650066
      NOWKB=0                                                           00660041
      ICLASS=0                                                          00670057
      IF(ITEST.EQ.3)NOWKB=1                                             00680041
      IF(ITEST.EQ.-1 .OR. ITEST.EQ. -2) NOWKB=1                         00690057
      IF(ITEST.EQ.-2) ICLASS=1                                          00700057
      LIMBAR=0                                                          00710061
      IF(LMAX.EQ.(-1))LIMBAR=1                                          00720053
      IF(LMAX.LT.-1)LUP=IABS(LMAX)                                      00730053
      DO 5 L=1,LUP                                                      00740038
 5    SIGML(L)=0.                                                       00750038
      DO 6 IO=1,17                                                      00760039
      SIGOR(IO)=0.                                                      00770045
 6    ROUTT(IO)=0.                                                      00780039
      ARED=AP*AT/(AP+AT)                                                00790079
      ECM=ELAB*AT/(AP+AT)                                               00800005
      IF(ITEST.GT.0 .AND. ITEST.LT.6)WRITE(6,10) ELAB,ECM               00810056
 10   FORMAT(' FUSION ELAB=',F8.1,4X,'ECM=',F8.1)                       00820018
      HW0=HBARC*HBARC/(AMU*ARED)                                        00830079
      SIG0=15.708*HW0/ECM                                               00840000
      ACN=AP+AT                                                         00850005
      ZCN=ZP+ZT                                                         00860005
      AN=ACN-ZCN                                                        00870005
      MPNT=5                                                            00880010
      MTIMES=1                                                          00890010
      F0=2./HW0                                                         00900010
C     LIMBAR=0                                                          00910060
      ONLYGS=0.                                                         00920030
C--------------------------------------   S2                            00930000
      K=MAX0(ION,1)                                                     00940000
      K=MIN0(ION,4)                                                     00950000
      IF(ABS(Q2).LT.0.01) K=4                                           00960068
      GOTO(30,32,34,36),K                                               00970002
 30   IO1=1                                                             00980000
      IO2=5                                                             00990000
      GOTO 38                                                           01000000
 32   IO1=6                                                             01010000
      IO2=15                                                            01020000
      GOTO 38                                                           01030000
 34   IO1=16                                                            01040000
      IO2=16                                                            01050000
      GOTO 38                                                           01060000
 36   IO1=17                                                            01070000
      IO2=17                                                            01080000
C---------------------------------------  S3                            01090000
 38   AP3=AP**0.333333                                                  01100034
      AT3=AT**0.333333                                                  01110034
      IO=IO1                                                            01120034
      SUM=0.                                                            01130000
      SUMCL=0.                                                          01140057
      L=0                                                               01150000
      IF(ZP.LE.0.)L=1                                                   01160039
      R0=R00                                                            01170034
      IF(R0.GT.0.) GOTO 50                                              01180064
      RPS=1.28*AP3-0.76+0.8/AP3                                         01190064
      RTS=1.28*AT3-0.76+0.8/AT3                                         01200064
      IF(IOPT.EQ.1 .OR. IOPT.EQ.2  .OR. IOPT.EQ.4 .OR. IOPT.EQ.5)       01210075
     1               GOTO 40                                            01220074
C     VLRO NEEDS SHARP RADII FOR OPTIONS 1,2,4,5(BASS)                  01230076
      RPS=RPS-1./RPS                                                    01240064
      RTS=RTS-1./RTS                                                    01250064
C     VLRO NEEDS CENTRAL RADII FOR OPTIONS -1(WS),0(EXP), AND 3(PROX)   01260076
 40   R0=(RPS+RTS)/(AP3+AT3)                                            01270074
 50   IF(AP.LT.3.01)R0=R00*AT3/(AT3+AP3)                                01280078
      RPS=R0*AP3                                                        01290064
      RTS=R0*AT3                                                        01300064
      IF(AP.LT.4.01 .AND. AP.GT.3.99)RPS=1.6                            01310064
      R=RPS+RTS                                                         01320064
      RPTS=R                                                            01330064
      RFUS=R0FUS*R/R0                                                   01340064
      OCOS=Z(IO)                                                        01350064
      V=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,          01360065
     1       XTH,APUSH,FPUSH,ECM,L,R,IOPT)                              01370065
      DR=0.2                                                            01380000
      DR2=0.1                                                           01390003
      RUP=50.*DR+R                                                      01400044
      R=R+DR                                                            01410000
      V1=VR(R)                                                          01420002
      IUP=50                                                            01430044
      DO 60 I=1,IUP                                                     01440039
      R=R+DR                                                            01450000
      V2=VR(R)                                                          01460002
      IF(V2.LT.V1 .AND.V.LT.V1) GOTO 70                                 01470000
      V=V1                                                              01480000
      V1=V2                                                             01490000
 60   CONTINUE                                                          01500000
 62   SIGF=0.                                                           01510003
      LMAX=0                                                            01520000
      WRITE(6,64) OCOS                                                  01530033
 64   FORMAT(' FUSION: POTENTIAL FOR L=0 HAS NO POCKET COS=',F8.3)      01540033
      WRITE(6,66) AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,                  01550045
     1      OCOS, V0,R0,D,ECM,L,R,IOPT                                  01560045
 66   FORMAT(4F6.0,F8.1,3I4/4F8.3,F8.1,I4,F8.2,I4)                      01570045
      ITEST=99                                                          01580056
      RETURN                                                            01590013
 70   X(3)=R                                                            01600000
      Y(3)=V2                                                           01610000
      X(2)=R-DR                                                         01620000
      Y(2)=V1                                                           01630000
      X(1)=X(2)-DR                                                      01640000
      Y(1)=V                                                            01650000
      IREP=0                                                            01660010
 71   CALL PARABO(X,Y,A1,RMAX,VLB)                                      01670010
      VBFUS=VLB                                                         01680069
      RBFUS=RMAX                                                        01690071
      IF(IREP.EQ.1) GOTO 72                                             01700010
      X(3)=RMAX+0.15                                                    01710010
      Y(3)=VR(X(3))                                                     01720010
      X(2)=RMAX                                                         01730010
      Y(2)=VR(RMAX)                                                     01740010
      X(1)=RMAX-0.15                                                    01750010
      Y(1)=VR(X(1))                                                     01760010
      IREP=1                                                            01770010
      GOTO 71                                                           01780010
 72   D2VLB=2. *A1(3)                                                   01790040
      HW=HW0*ABS(D2VLB)                                                 01800005
      HW=SQRT(HW)                                                       01810040
      HWFUS=HW                                                          01820071
      HWMX(IO)=HW                                                       01830002
      VMX(IO)=VLB                                                       01840002
      RMX(IO)=RMAX                                                      01850000
      DUM=       6.28319*(VLB-ECM)/HW                                   01860005
      IF(ECM.LT.VLB)DUM=DUM*SQRT(CRED)                                  01870043
      IF(ABS(DUM).LT.AMXEXP)T0(IO)=1./(1.+EXP(DUM))                     01880055
      IF(DUM.GT.AMXEXP)T0(IO)=0.                                        01890055
      IF(DUM.LT.(-AMXEXP))T0(IO)=1.                                     01900055
      TPARAP(IO)=T0(IO)                                                 01910025
      IF(NOWKB.EQ.1) GOTO 74                                            01920040
      IWKB(IO)=0                                                        01930010
      IF(T0(IO).GT.0.49) GOTO 74                                        01940031
      IWKB(IO)=1                                                        01950020
      IF(T0(IO).EQ.0.) GOTO 74                                          01960010
      RIN=RMAX                                                          01970011
      ROUT=1.44*ZP*ZT/ECM-2.                                            01980033
      IF(ZP.LE.0.)ROUT=SQRT(HW0/ECM)-3.                                 01990039
      ROUT=DMAX1(RMAX,ROUT)                                             02000033
      DR=0.1                                                            02010027
      DUM =TWKB(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,V0,R0,D,XTH,APUSH,  02020065
     1          FPUSH,DR,OCOS,ECM,L,RIN,ROUT,RFUS,CRED,IOPT)            02030065
      T0(IO)=0.                                                         02040013
      IF(DUM.LE.0.)GOTO 73                                              02050013
      T0(IO)=1./(1.+(1./DUM))                                           02060010
 73   RINN(IO)=RIN                                                      02070013
      ROUTT(IO)=ROUT                                                    02080011
      IF(ITEST.GT.1 .AND.ITEST.LT.6)                                    02090056
     1WRITE(6,7300) ECM,OCOS,RMAX,RPTS,RIN,ROUT,RFUS                    02100056
 7300 FORMAT(' WKB L=0:ECM,OCOS,RMAX,RPTS,RIN,ROUT,RFUS',7F9.3)         02110016
 74   IF(ICLASS.EQ.0) GOTO 75                                           02120057
      IF(T0(IO).GE.0.5) TCL=1.                                          02130057
      IF(T0(IO).LT.0.5) TCL=0.                                          02140057
      SUMCL= SUMCL+WEIGHT(IO)*TCL                                       02150057
 75   SUM=SUM+WEIGHT(IO)*T0(IO)                                         02160057
 7500 IF(IO.GE.IO2) GOTO 90                                             02170057
      IO=IO+1                                                           02180000
      OCOS=Z(IO)                                                        02190000
C--------------------------------------                                 02200003
      IF(Q2.LT.0) GOTO 80                                               02210003
      R=RMX(IO-1)-DR2                                                   02220008
      V=VO(L,R,ECM,OCOS)                                                02230008
      R=R+DR                                                            02240003
      V1=VR(R)                                                          02250003
      R=R+DR                                                            02260003
      V2=VR(R)                                                          02270003
 76   IF(V2.LT.V1 .AND. V.LT.V1)GOTO 70                                 02280003
      V=V1                                                              02290003
      V1=V2                                                             02300003
      R=R+DR                                                            02310003
      IF(R.GT.RUP) GOTO 62                                              02320003
      V2=VR(R)                                                          02330003
      GOTO 76                                                           02340003
 80   R=RMX(IO-1)+DR2                                                   02350008
      V2=VO(L,R,ECM,OCOS)                                               02360008
      R=R-DR                                                            02370003
      V1=VR(R)                                                          02380003
      R=R-DR                                                            02390003
      V=VR(R)                                                           02400003
 86   IF(V2.LT.V1 .AND. V.LT.V1) GOTO 88                                02410005
      V2=V1                                                             02420003
      V1=V                                                              02430003
      R=R-DR                                                            02440003
      IF(R.LT.RPS) GOTO 62                                              02450003
      V=VR(R)                                                           02460003
      GOTO 86                                                           02470003
 88   R=R+DR+DR                                                         02480005
      GOTO 70                                                           02490005
C-------------------------------------                                  02500005
 90   T=SUM*0.5                                                         02510003
      TCL=SUMCL*0.5                                                     02520057
      T00=0.01*T                                                        02530010
      IF(ITEST.LT.2) GOTO 105                                           02540005
      IF(IO2.GT.15) GOTO 95                                             02550002
      WRITE(6,101) T,L                                                  02560005
      WRITE(6,100) (VMX(IO),IO=IO1,IO2),(RMX(IO),IO=IO1,IO2),           02570002
     1             (HWMX(IO),IO=IO1,IO2)                                02580005
      WRITE(6,103) (T0(IO),IO=IO1,IO2)                                  02590005
      WRITE(6,103) (TPARAP(IO),IO=IO1,IO2)                              02600025
      DO 92 IO=IO1,IO2                                                  02610045
 92   SIGOR(IO)=SIG0*T0(IO)                                             02620045
      GOTO 105                                                          02630002
 95   WRITE(6,102) VMX(IO),RMX(IO),HW,T,L                               02640002
 100  FORMAT(5F12.2)                                                    02650002
 103  FORMAT(5E12.3)                                                    02660005
 101  FORMAT( E12.3,I4)                                                 02670025
 102  FORMAT(3F12.3,E12.3,I4)                                           02680002
C----------------------------------------------                         02690028
 105  IFLAG=0                                                           02700002
      SIGML(1)=SIG0*T                                                   02710039
      IF(ICLASS.EQ.1) SIGML(1)=SIG0*TCL                                 02720057
      IF(ZP.GT.0.)GOTO 106                                              02730047
      TN=1.-EXP(-1.085*SQRT(ECM))                                       02740046
      IF(ITEST.GT.1 .AND. ITEST.LT.6) WRITE(6,107) ECM,TN               02750056
 107  FORMAT( ' ECM,T0',F8.2,E12.3)                                     02760047
      SIGML(1)=SIG0*TN                                                  02770047
 106  SIGF=SIGML(1)                                                     02780047
C     IF(SIGF.LT..1E-6)GOTO 400                                         02790055
      T=1.                                                              02800000
      DO 110 IO=IO1,IO2                                                 02810028
      TMX(IO)=T0(IO)                                                    02820028
 110  T0(IO)=0.01*T0(IO)                                                02830028
C----------------------------------------PARTIAL WAVE LOOP L            02840028
      DO 250 L=1,LUP                                                    02850005
      AL=L                                                              02860005
      IF(LIMBAR.EQ.0) GOTO 120                                          02870005
      CALL FISROT(ACN,ZCN,AN,AL,DELR,DELSP,ERO,ONLYGS)                  02880030
      BAR=DELSP-DELR                                                    02890005
      IF(BAR.GT.0.1) GOTO 120                                           02900059
      FLGRAZ=AL                                                         02910005
      IFLAG=2                                                           02920005
      GOTO 260                                                          02930005
 120  T1=T                                                              02940005
      SUM=0.                                                            02950000
      SUMCL=0.                                                          02960058
C
C     I.Giese 2007
C     Variable R12 bekommt nie einen Wert, also im Linux auf 0.0 setzen
      R12 = 0.0D0
C
C------------------------------------ORIENTATION LOOP IO                02970028
      DO 150 IO=IO1,IO2                                                 02980000
      IF(TMX(IO).LT.TMIN) GOTO 125                                      02990055
      IF(LUP.LT.200) GOTO 130                                           03000037
      IF(TMX(IO).GT.T0(IO)) GOTO 130                                    03010028
 125  TMX(IO)=0.                                                        03020037
      GOTO 150                                                          03030028
 130  OCOS=Z(IO)                                                        03040028
      R=RMX(IO)                                                         03050000
      X(2)=R                                                            03060003
      Y(2)=VO(L,R,ECM,OCOS)                                             03070008
      R=R+0.15                                                          03080047
      X(3)=R                                                            03090000
      Y(3)=VR(R)                                                        03100002
      R=R-0.30                                                          03110047
 132  X(1)=R                                                            03120047
      Y(1)=VR(R)                                                        03130047
      CALL PARABO(X,Y,A1,RMAX,VLB)                                      03140000
      IF(RMAX.GT.X(3))GOTO 135                                          03150051
      IF(RMAX.GT.X(1))GOTO 134                                          03160051
      X(3)=X(2)                                                         03170047
      Y(3)=Y(2)                                                         03180047
      X(2)=X(1)                                                         03190047
      Y(2)=Y(1)                                                         03200047
      R=R-0.15                                                          03210047
      IF(R.LT.R12) GOTO 135                                             03220047
      GOTO 132                                                          03230047
 134  D2VLB= 2.*A1(3)                                                   03240047
      HW=HW0*ABS(D2VLB)                                                 03250005
      HW=SQRT(HW)                                                       03260040
      IF(A1(3).GT.0.) GOTO 135                                          03270033
      DUM=       6.28319*(VLB-ECM)/HW                                   03280005
      IF(ECM.LT.VLB)DUM=DUM*SQRT(CRED)                                  03290043
      IF(ABS(DUM).LT.AMXEXP)T =1./(1.+EXP(DUM))                         03300055
      IF(DUM.GT.AMXEXP)T =0.                                            03310055
      IF(DUM.LT.(-AMXEXP))T =1.                                         03320055
      JWKB=0                                                            03330048
      GOTO 136                                                          03340033
 135  T=0.                                                              03350033
      IF(ITEST.GT.1 .AND. ITEST.LT.6) WRITE(6,1350) OCOS,L              03360056
 1350 FORMAT(' OCOS,L=',F8.3,I4,4X,'PROBABLY NO POCKET')                03370051
 136  RMX(IO)=RMAX                                                      03380033
      VMX(IO)=VLB                                                       03390002
      TPARAP(IO)=0.                                                     03400045
      IF(T.LE.0.)GOTO 145                                               03410033
      TPARAP(IO)=T                                                      03420025
      IF(NOWKB.EQ.1) GOTO 145                                           03430040
      IF(T.GT.0.49) GOTO 145                                            03440031
      IF(T.EQ.0.) GOTO 145                                              03450020
C---------------------------                                            03460039
      IF(ROUTT(IO).LE.0.) GOTO 138                                      03470039
      IF(RINN (IO).LE.0.) GOTO 138                                      03480050
      DUM=(ROUTT(IO)-RINN(IO))*0.1                                      03490033
      DR=DMIN1(AMPAR1,DUM)                                              03500033
      DR=DMAX1(AMPAR2,DR)                                               03510033
      RIN=RINN(IO)+DR                                                   03520033
      RIN=DMIN1(RMX(IO),RIN)                                            03530033
      ROUT=ROUTT(IO)-DR                                                 03540033
      ROUT=DMAX1(ROUT,RMX(IO))                                          03550033
      GOTO 139                                                          03560047
 138  RIN=          RMX(IO)                                             03570039
      ROUT=RIN                                                          03580039
      DR=0.1                                                            03590039
 139  IF(ZP.GT.0.) GOTO 140                                             03600047
      ROUT=SQRT(HW0*L*(L+1)/(2.*ECM))-0.5                               03610047
 140  IF(R0FUS.LE.0.)RFUS=-RINN(IO)                                     03620039
         T=TWKB(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,V0,R0,D,XTH,APUSH,  03630065
     1          FPUSH,DR,OCOS,ECM,L,RIN,ROUT,RFUS,CRED,IOPT)            03640065
      IF(T.GT.0.)T=1./(1.+(1./T))                                       03650025
      IF(T.LE.0.0 .AND. TPARAP(IO).GT.0.05) T=TPARAP(IO)                03660026
      RINN(IO)=RIN                                                      03670011
      ROUTT(IO)=ROUT                                                    03680011
      JWKB=1                                                            03690048
      IWKB(IO)=1.                                                       03700022
 145  IF(ICLASS.EQ.0) GOTO 146                                          03710057
      IF(T.GE.0.5) TCL=1.                                               03720057
      IF(T.LT.0.5) TCL=0.                                               03730057
      SUMCL= SUMCL+WEIGHT(IO)*TCL                                       03740057
 146  TMX(IO)=T                                                         03750057
      HWMX(IO)=HW                                                       03760002
      SUM=SUM+WEIGHT(IO)*T                                              03770000
 150  CONTINUE                                                          03780000
C--------------------------------------- END ORIENTATION LOOP IO        03790028
      L1=L+1                                                            03800000
      T=SUM*0.5                                                         03810000
      TCL=SUMCL*0.5                                                     03820057
      IF(ITEST.LT.2) GOTO 205                                           03830005
      IF(IO2.GT.15) GOTO195                                             03840002
      WRITE(6,101) T,L                                                  03850025
      WRITE(6,100) (VMX(IO),IO=IO1,IO2),(RMX(IO),IO=IO1,IO2),           03860002
     1             (RINN(IO),IO=IO1,IO2),                               03870027
     2            (ROUTT(IO),IO=IO1,IO2),                               03880027
     3             (HWMX(IO),IO=IO1,IO2)                                03890027
      WRITE(6,103) (TMX(IO),IO=IO1,IO2)                                 03900009
      WRITE(6,103) (TPARAP(IO),IO=IO1,IO2)                              03910025
      DUM=SIG0*(2*L+1)                                                  03920045
      DO 170 IO=IO1,IO2                                                 03930045
 170  SIGOR(IO)=SIGOR(IO)+DUM*TMX(IO)                                   03940045
      GOTO 205                                                          03950002
 195  IF(JWKB.EQ.0)                                                     03960048
     1WRITE(6,102) VMX(IO1),RMX(IO1),HW,T,L                             03970048
      IF(JWKB.EQ.1)                                                     03980048
     1WRITE(6,202) VMX(IO1),RMX(IO1),HW,T,L,RIN,ROUT                    03990048
 202  FORMAT(3F12.3,E12.3,I4,2X,'RIN,ROUT=',2F8.2)                      04000048
 205  IF(T.LT..5 .AND.T1.GT.0.5)FLGRAZ=L-1                              04010002
      IF(T.LT.TMIN) GOTO 260                                            04020055
      IF(LUP.LT.200)GOTO 210                                            04030036
      IF(T.LT.T00)        GOTO 260                                      04040007
 210  SIGML(L1)=SIG0*(2*L+1)   *T                                       04050036
      IF(ICLASS.EQ.1) SIGML(L1)=SIG0*(2*L+1)   *TCL                     04060057
      IF(ZP.LE.0.)SIGML(L1)=SIGML(L1)*TN                                04070047
      SIGF=SIGF+SIGML(L1)                                               04080000
 250  CONTINUE                                                          04090000
C--------------------------------------------END L LOOP                 04100005
 255  IFLAG=1                                                           04110000
 260  LMAX=MIN0(L,LUP)                                                  04120049
      IF(ITEST.LT.1)GOTO 999                                            04130005
      WRITE(6,301) ELAB,AP,ZP,AT,ZT,   V0,R0,D,Q2,CRED,NOPROX,NOCURV,   04140031
     1             ION                                                  04150002
 301  FORMAT(/' FUSION'/' ELAB,AP,ZP,AT,ZT',F10.3,F8.0,3F6.0,4X,        04160031
     1        'V0,R0,D,Q2,CRED',5F10.3/  ' NOPROX,NOCURV,ION',3I2)      04170032
      WRITE(6,302) SIGF,LMAX,FLGRAZ,IFLAG                               04180002
 302  FORMAT(' SIGFUS,LMAX,LGRAZ,FLAG',E12.3,I10,F10.0,I10/             04190031
     1      ' PARTIAL SIG')                                             04200032
      WRITE(6,303) (SIGML(L),L=1,LMAX)                                  04210002
 303  FORMAT(10E12.3)                                                   04220002
      IF(ITEST.LT.2) GOTO 999                                           04230045
      IF(IO1.GT.15) GOTO 999                                            04240045
      WRITE(6,310)(SIGOR(IO),IO=IO1,IO2)                                04250045
 310  FORMAT(/' SIGFUS FOR ORIENT.'/5E12.3/5E12.3)                      04260045
      GOTO 999                                                          04270010
C400  LMAX=0                                                            04280070
C     FLGRAZ=0.                                                         04290070
C     DO 410 L=1,LUP                                                    04300070
C410  SIGML(L)=0.                                                       04310070
C     SIGF=0.                                                           04320070
C     WRITE(6,420) ELAB                                                 04330054
C420  FORMAT(//F8.2,' MEV IS WELL BELOW COULOMB BARRIER,SIGFUS=0.'//)   04340070
 999  RETURN                                                            04350000
      END                                                               04360000
