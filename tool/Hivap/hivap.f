      PROGRAM HIVAP                                                     00020053
C
C  Programm 'HIVAP' zur Berechnung von Fusions- und Verdampfungs-
C  restkernquerschnitten
C  Autor: W. Reisdorf, GSI
C  modifiziert zur Verwendung auf der 'rzri6f' durch F.P.Hessberger
C  Information ueber Modifikationen auf 
C     rzri6f:/u/hess/fortpro/hivapn/aainfo.dat
C  Option zur Behinderung der Spaltung eingebaut; laeuft bis ca. 130 MeV
C  Anregungsenergie, dann gibt es offenbar Schwierigkeiten mit den
C  Zustandsdichten ( >3x10**39)
C  Stand: 7.6.1994
C                                                                       00040000
      INCLUDE 'common.f'                                                00050054
      COMMON/DELFIS/EDELFIS,ADELFIS,INOF 
      COMMON/PUSHPA/SPUSH,SCREIPA
      DIMENSION ARRAY(40),SHEL(5)                                       00060059
      CHARACTER*1 AW(72)                                                00070048
      CHARACTER*4 BLNK                                                  00080052
      CHARACTER*10 BLNK2
      CHARACTER*1 SZ
      CHARACTER*6 AUT1  
      CHARACTER*6 AUT2
      CHARACTER*35 AUTOR
      CHARACTER*12 BEARB
      DATA BLNK/'    '/ BLNK2/'          '/                             00090000
      DATA AUT1/'F.P.He'/,AUT2/'berger'/,BEARB/'Bearbeitung '/
C
C     UNUSED: AW
C     MAYBE UNINITIALIZED: JJ,ICHECK,CUTPOP,SIGFU
      JJ=0
      ICHECK=0
      CUTPOP=0.0D0
      SIGFU=0.0D0
      SZ=CHAR(223)
      AUTOR=BLNK2//BEARB//AUT1//SZ//AUT2
C-----------------------------------------------------------------------00100000
C                                                                       00110000
      SCREIPA=0. 
      AMPAR1=0.
      INOF=0
      LEVPRN=0                                                          00120047
      CALL ALLOCH                                                       00150000
      EPS=0.01                                                          00180000
      IPOPS=4096                                                        00190000
C     DIMENSION OF POP,POPN,POPP,POPA                                   00200000
      MASSES=9                                                          00210000
C     LOGICAL UNIT FOR MASS EXCESS TABLE                                00220000
      IYR=14                                                            00230000
C     LOGICAL UNIT FOR TRANSMISSION COEFF IS IUNIT=21                   00240000
C     LOGICAL UNIT FOR INPUT SUMMARY IS 23                              00250000
      IDEN=2500                                                         00260000
      IDNA=5000                                                         00270000
      NAUP=25                                                           00280000
      NZUP=25                                                           00290000
      NAUP2=NAUP+2                                                      00300000
      NZUP2=NZUP+2                                                      00310000
      NUMIS=5                                                           00320000
      JDIM=100                                                          00330000
C     ALSO DEFINED IN 'TOT'                                             00340000
      IEDIM=400                                                         00350000
C     ALSO DEFINED IN 'TOT'  (SPECE)                                    00360000
      KEN=25                                                            00370000
      KEP=25                                                            00380000
      KEP1=20+KEP                                                       00390000
      KEG=20                                                            00400000
      KEA=35                                                            00410000
      KEA1=30+KEA                                                       00420000
      KEF=20                                                            00430000
      LN=12                                                             00440000
      LP=12                                                             00450000
      LA=20                                                             00460000
      CUT=0.                                                            00470062
      FRACT2=0.                                                         00480062
      ABSMIN=0.                                                         00490062
      PRCN=0.                                                           00500062
      SIGLOW=0.                                                         00510062
      VRN=0.                                                            00520000
      R0RN=0.                                                           00530000
      ADIFRN=0.                                                         00540000
      VIN=0.                                                            00550000
      R0IN=0.                                                           00560000
      ADIFIN=0.                                                         00570000
      VRP=0.                                                            00580000
      R0RP=0.                                                           00590000
      ADIFRP=0.                                                         00600000
      VIP=0.                                                            00610000
      R0IP=0.                                                           00620000
      ADIFIP=0.                                                         00630000
      RCLMBP=0.                                                         00640000
      CBFACP=1.                                                         00650000
      VRA=0.                                                            00660000
      R0RA=0.                                                           00670000
      ADIFRA=0.                                                         00680000
      VIA=0.                                                            00690000
      R0IA=0.                                                           00700000
      ADIFIA=0.                                                         00710000
      RCLMBA=0.                                                         00720000
      CBFACA=1.                                                         00730000
      SIGR0=0.                                                          00740009
      NOLEP=1                                                           00750023
      LEP=1                                                             00760073
      NOLJI=1                                                           00770023
      NOLJF=1                                                           00780023
      IPTR=0                                                            00790062
      WRITE(6,10)                                                       00800000
 10   FORMAT(10X,'Programm HIVAP,  Autor  W.Reisdorf  ',
     1' Version fuer "rzri6f" vom  7. 6. 1994'/)   
      WRITE(6,*) AUTOR
      WRITE(6,*) BLNK2
C                                                                       00820000
C******************************************************************     00830000
C                                                 NEW CASE              00840000
 99   CONTINUE                                                          00850000
 100  FORMAT(18A4)                                                      00870000
 101  FORMAT(18I4)                                                      00880000
C 102  FORMAT(9F8.3)                                                    00890000
      DO 105 K=1,IPOPS                                                  00900000
      POP (K)=0.                                                        00910000
      POPN(K)=0.                                                        00920000
      POPP(K)=0.                                                        00930000
      POPA(K)=0.                                                        00940000
      PENTRY(K)=0.                                                      00950000
 105  CONTINUE                                                          00960000
      DO 107 IE=1,NENMX                                                 00970005
      DO 107 NSTEP=1,NAUP                                               00980005
      DO 107 NROW=1,NZUP                                                00990005
 107  XSECTE(NROW,NSTEP,IE)=0.                                          01000005
      READ(5, 100,END=9999)(TITLE(I),I=1,18)                            01010021
      WRITE(6,100)(TITLE(I),I=1,18)                                     01020000
      WRITE( 23  ,100)(TITLE(I),I=1,18)                                 01030000
      WRITE(6,*) BLNK2
      NOMORE=0                                                          01040037
      DO 110 I=1,18                                                     01050000
      IF(TITLE(I).NE.BLNK) GOTO 120                                     01060000
 110  CONTINUE                                                          01070000
      GOTO 9999                                                         01080000
C---------------------------------------------------------------        01090000
 120  CALL MYRD(ARRAY,18,23,5)                                          01100055
      MPROJ = ARRAY(1 )+EPS                                             01110000
      IZP   = ARRAY(2 )+EPS                                             01120000
      MTGT  = ARRAY(3 )+EPS                                             01130000
      IZT   = ARRAY(4 )+EPS                                             01140000
      ISHELL= ARRAY(5 )+EPS                                             01150000
      ISHELF= ARRAY(6 )+EPS                                             01160000
      IPAIR = ARRAY(7 )+EPS                                             01170000
      IF(ISHELL.EQ.0) IPAIR=0                                           01180030
      MC    = ARRAY(8 )+EPS                                             01190000
      IF(ISHELL.EQ.0) MC=1                                              01200030
      MP    = ARRAY(9 )+EPS                                             01210000
      IF(ISHELL.EQ.0) MP=1                                              01220030
      IBF   = ARRAY(10)+EPS                                             01230000
      IFISRT= ARRAY(11)+EPS                                             01240000
      NOF   = ARRAY(12)+EPS                                             01250000
      NON   = ARRAY(13)+EPS                                             01260000
      NOP   = ARRAY(14)+EPS                                             01270000
      NOA   = ARRAY(15)+EPS                                             01280000
      NOG   = ARRAY(16)+EPS                                             01290000
      IDISC = ARRAY(17)+EPS                                             01300000
      IGAM  = ARRAY(18)+EPS                                             01310000
C-----------------------------                                          01320006
C  Eingefuegt  fuer Verzoegerung der Spaltung  NOF=2, 2.5.94 FPH
C-----------------------------
      IF(NOF.NE.2) GOTO 1201
      INOF=1
      NOF=0
      CALL MYRD(ARRAY,2,23,5)
      EDELFIS=ARRAY(1)+EPS
      ADELFIS=ARRAY(2)+EPS
C-----------------------------  
 1201 IF(IDISC.NE.0) WRITE(25, 9993)(TITLE(I),I=1,18)                   01330033
      NZ1=0                                                             01340000
      ACN=MTGT+MPROJ                                                    01350000
      AP=MPROJ                                                          01360000
      ZP=IZP                                                            01370000
      AT=MTGT                                                           01380000
      ZT=IZT                                                            01390000
      IACN=ACN                                                          01400000
      ZCN=IZT+IZP                                                       01410000
      IZCN=ZCN                                                          01420000
      NOPUT=0                                                           01430000
      IF(NZ1.LT.0) NOPUT=1                                              01440000
      NZ1=IABS(NZ1)                                                     01450000
      NZ1=MAX0(1,NZ1)                                                   01460000
C-----------------------------                                          01470075
      CALL MYRD(ARRAY,18,23,5)                                          01480055
C-----------------------------                                          01490075
      NA    = ARRAY(1 ) + EPS                                           01500000
      NZ    = ARRAY(2 ) + EPS                                           01510000
      MASSES= ARRAY(3 ) + EPS                                           01520000
      NUMB  = ARRAY(4 ) + EPS                                           01530000
      IOVER = ABS(ARRAY(5))+EPS                                         01540000
      IF(ARRAY(5).LT.0.) IOVER=-IOVER                                   01550000
      INERT = ARRAY(6 ) + EPS                                           01560000
      INERF = ARRAY(7 ) + EPS                                           01570000
      FINERT= ARRAY(8 )                                                 01580000
      IF(FINERT.LE.0.)FINERT=1.                                         01590000
      ILIM  = ARRAY(9 ) + EPS                                           01600000
      LPRINT= ARRAY(10) + EPS                                           01610000
      LOGUN = ARRAY(11) + EPS                                           01620000
      ICOR  = ARRAY(12) + EPS                                           01630000
      NUMISO= ARRAY(13) + EPS                                           01640000
      LDBM  = ABS(ARRAY(14)) + EPS                                      01650000
      IF(ARRAY(14).LT.0.) LDBM=-LDBM                                    01660000
      LDBM5=0                                                           01670086
      IF(LDBM.EQ.5) LDBM5=1                                             01680086
      IRAST = ARRAY(15) + EPS                                           01690000
      IOWKB = ARRAY(16) + EPS                                           01700000
      ITRANS= ARRAY(17) + EPS                                           01710000
      JFJI  = ARRAY(18) + EPS                                           01720000
      IF(NA.EQ.0) NA=14                                                 01730000
      NZ=MAX0(1,NZ)                                                     01740000
      NZ=MIN0(14,NZ)                                                    01750000
      LOGUN=0                                                           01770000
      IF(MASSES.NE.0)MASSES=9                                           01780000
      DEL(5)=0.                                                         01790000
      DO 130 K=1,5                                                      01800000
 130  DJ(K)=0.                                                          01810000
      IF(JFJI.NE.1) GOTO 132                                            01820000
C------------------------------------                                   01830075
      CALL MYRD(ARRAY,3,23,5)                                           01840055
C------------------------------------                                   01850075
      DJ(2)=ARRAY(1)                                                    01860000
      DJ(3)=ARRAY(2)                                                    01870000
      DJ(4)=ARRAY(3)                                                    01880000
C      IF(JFJI.EQ.1) WRITE(6,102) DJ(2),DJ(3),DJ(4)                     01890000
C-----------------------------------------------------------------------01900000
C                                                  GAMMA NORM           01910000
 132  IOUT=LPRINT                                                       01920077
      CALL GAMMAS(ACN,IOUT,DELG,IGAM,NOE2,STRIPE,IOPT)                  01930077
      DEL(1)=DELG                                                       01940071
C 1320 WRITE(6,1321) DEL(1)                                              01950072
C 1321 FORMAT(' DEL(1) AFTER GAMMAS',F8.1)                               01960071
C-----------------------------------------------------------------------01970000
C                                       LEVEL DENSITY PARAMETERS        01980000
      CALL MYRD(ARRAY,8,23,5)                                           01990055
C------------------------------                                         02000075
      CST   =ARRAY(1 )                                                  02010000
      CLD   =ARRAY(2 )                                                  02020000
      BARFAC=ARRAY(3 )                                                  02030000
      ESHELL=ARRAY(4 )                                                  02040000
      BAR0  =ARRAY(5 )                                                  02050000
      SHELL0=ARRAY(6 )                                                  02060000
      DELT  =ARRAY(7 )                                                  02070000
      QVALUE=ARRAY(8 )                                                  02080000
C------------------------------                                         02090016
      IPREEX=1                                                          02100086
      AX=0.                                                             02110086
      BETA0=0.                                                          02120086
      UCRIT=0.                                                          02130086
      EDCOLL=0.                                                         02140086
      IENH=0                                                            02150086
      IF(LDBM5.EQ.1) THEN                                               02160086
        CALL MYRD(ARRAY,7,23,5)                                         02170086
        LDBM=ABS(ARRAY(1))+EPS                                          02180086
        IF(ARRAY(1).LT.0.) LDBM=-LDBM                                   02190086
        IPREEX=ARRAY(2)+EPS                                             02200086
        AX=ARRAY(3)                                                     02210086
        IENH=ARRAY(4)+EPS                                               02220086
        BETA0=ARRAY(5)                                                  02230086
        EDCOLL=ARRAY(6)                                                 02240086
        UCRIT=ARRAY(7)                                                  02250086
      ENDIF                                                             02260086
 133  IF(BARFAC.LE.0.) BARFAC=1.                                        02270084
      IF(ESHELL.LE.0.) ESHELL=18.5                                      02280000
      IF(DELT.LE.0.)DELT=11.                                            02290000
 135  FORMAT(4F8.3,6I4)                                                 02300000
      IF(UCRIT.LE.0.)UCRIT=5.                                           02310038
      KHS=0                                                             02320044
C     KHS=1 IS PAIRING ACCORDING TO KHS, NOT YET TESTED                 02330045
      IF(IENH.NE.0) AX=1.                                               02340045
      IF(TZERO.EQ.0.)TZERO=0.2                                          02350000
      IF(LDBM.NE.(-1)) GOTO 220                                         02360000
      IRC=NA*NZ                                                         02370008
      CALL MYRD(ALS,IRC,23,5)                                           02380055
      IRC=NA*NZ                                                         02390008
      CALL MYRD(DELTAS,IRC,23,5)                                        02400055
C-------------------------------------------------------------          02410000
C                                         CUT-OFF PARAMETERS            02420000
 220  IF(ILIM.NE.1)GOTO 226                                             02430000
      CALL MYRD(ARRAY,9,23,5)                                           02440055
      CUT   =ARRAY(1 )                                                  02450000
      FRACT2=ARRAY(2 )                                                  02460000
      ABSMIN=ARRAY(3 )                                                  02470000
      PRCN  =ARRAY(4 )                                                  02480000
      SIGLOW=ARRAY(5 )                                                  02490000
C     DEL(1)=ARRAY(6 )                                                  02500064
      DUM   =ARRAY(6 )                                                  02510064
      NOLEP =ARRAY(7 )+EPS                                              02520000
      NOLJI =ARRAY(8 )+EPS                                              02530000
      NOLJF =ARRAY(9 )+EPS                                              02540023
      IF(JFJI.EQ.1)NOLEP=1                                              02550023
      IF(JFJI.EQ.1)NOLJI=1                                              02560023
      IF(JFJI.EQ.1)NOLJF=1                                              02570023
C      IF(ILIM.EQ.1)WRITE(6,225) CUT,FRACT2,ABSMIN,PRCN ,SIGLOW,DEL(1), 02580000
C     1                          NOLEP,NOLJI,NOLJF                      02590000
C 225  FORMAT(' CUT,FRACT2,ABSMIN=',3E12.4/' PRCN,SIGLOW,DEL1',3F8.3,4X,02600010
C     1       'NOLEP,NOLJI,NOLJF=',3I4)                                 02610006
 226  IF(CUT.LE.0. .OR. CUT.GT.0.15)CUT=0.1E-2                          02620000
      PRCNT=1.-PRCN*0.01                                                02630000
      IF(PRCN.LE.0. .OR. PRCN.GT.20.) PRCNT=1.-CUT                      02640000
      IF(FRACT2.LE.0.)FRACT2=0.001                                      02650000
      IF(DEL(1).LE.0.)DEL(1)=1.                                         02660000
C-------------------------------------------------------------          02670000
C                                       TRANSMISSION                    02680000
      TLOW=0.                                                           02690064
      DEL(2)=0.1                                                        02700064
C      WRITE(6,2261) DEL(2)                                             02710071
C 2261 FORMAT(' DEL(2)=',F8.2)                                          02720071
      QT=0.                                                             02730064
      QQ2=0.                                                            02740064
      IOMN=0                                                            02750064
      IOMP=0                                                            02760064
      IOMA=0                                                            02770064
      IPTR=0                                                            02780064
      IF(IOVER.GE.0) GOTO 229                                           02790000
       WRITE(6,2260)                                                    02800071
 2260 FORMAT(' IOVER nicht Standart OM')                                02810057
      CALL MYRD(ARRAY,11,23,5)                                          02820055
      IOVER =ARRAY(1 ) + EPS                                            02830000
      IOWKB =ARRAY(2 ) + EPS                                            02840000
      ITRANS=ARRAY(3 ) + EPS                                            02850000
      SIGLOW=ARRAY(4 )                                                  02860000
      TLOW  =ARRAY(5 )                                                  02870000
      DEL(2)=ARRAY(6)                                                   02880000
      QT    =ARRAY(7 )                                                  02890000
      QQ2=QT                                                            02900000
      IOMN=ARRAY(8)+EPS                                                 02910000
      IOMP=ARRAY(9)+EPS                                                 02920000
      IOMA=ARRAY(10)+EPS                                                02930000
      IPTR  =ARRAY(11) + EPS                                            02940000
C------------------------------  OPT MOD N                              02950000
      IF(IOMN.EQ.0) GOTO 2280                                           02960000
      CALL MYRD(ARRAY,6,23,5)                                           02970055
      VRN   = ARRAY(1)                                                  02980000
      R0RN  = ARRAY(2)                                                  02990000
      ADIFRN= ARRAY(3)                                                  03000000
      VIN   = ARRAY(4)                                                  03010000
      R0IN  = ARRAY(5)                                                  03020000
      ADIFIN= ARRAY(6)                                                  03030000
 2280 IF(IOMP.EQ.0) GOTO 2281                                           03040000
      CALL MYRD(ARRAY,8,23,5)                                           03050055
      VRP   = ARRAY(1)                                                  03060000
      R0RP  = ARRAY(2)                                                  03070000
      ADIFRP= ARRAY(3)                                                  03080000
      VIP   = ARRAY(4)                                                  03090000
      R0IP  = ARRAY(5)                                                  03100000
      ADIFIP= ARRAY(6)                                                  03110000
      RCLMBP= ARRAY(7)                                                  03120000
      CBFACP= ARRAY(8)                                                  03130000
 2281 IF(IOMA.EQ.0) GOTO 229                                            03140000
      CALL MYRD(ARRAY,8,23,5)                                           03150055
      VRA   = ARRAY(1)                                                  03160000
      R0RA  = ARRAY(2)                                                  03170000
      ADIFRA= ARRAY(3)                                                  03180000
      VIA   = ARRAY(4)                                                  03190000
      R0IA  = ARRAY(5)                                                  03200000
      ADIFIA= ARRAY(6)                                                  03210000
      RCLMBA= ARRAY(7)                                                  03220000
      CBFACA= ARRAY(8)                                                  03230000
 229  IF(SIGLOW.EQ.0.)SIGLOW=0.001                                      03240000
C-------------------------------------------------------------          03250000
C                                       ISOMERS                         03260000
      IF(NUMISO.EQ.0) GOTO  230                                         03270000
      NUMISO=MIN0(NUMIS,NUMISO)                                         03280000
      IRC=NUMISO*3                                                      03290000
      ISO=0                                                             03300000
      CALL MYRD(ARRAY,IRC,23,5)                                         03310055
      DO 2290 I=1,IRC,3                                                 03320000
      ISO=ISO+1                                                         03330000
      ISOZ(ISO)=ARRAY(I  )+EPS                                          03340000
      ISOA(ISO)=ARRAY(I+1)+EPS                                          03350000
 2290 ISOJ(ISO)=ARRAY(I+2)+EPS                                          03360000
      IF(LPRINT.EQ.0) WRITE(6,2291) NUMISO,(ISOZ(I),ISOA(I),ISOJ(I),    03370057
     1                I=1,NUMISO)                                       03380057
 2291 FORMAT(' NUMISO',I3,4X,12I4)                                      03390057
C-------------------------------------------------------------          03400000
C                                        INERTIAS                       03410000
 230  IF(INERT.EQ.0) GOTO 255                                           03420000
      CALL MYRD(RATIOS,6,23,5)                                          03430055
      IF(LPRINT.EQ.0) WRITE(6,2301) RATIOS                              03440057
 2301 FORMAT(' RATIOS',6F8.3)                                           03450057
C                                                                       03480000
C-----------------------------------------------------------------------03490000
C                              BARRIERS AND    SEPARATION ENERGIES      03500000
  255 CALL SEPENF(MTGT,IZT,MPROJ,IZP,QC,NA,NZ,POP,3,ISHELL,SHELLS,      03510000
     1            SHELL0,DELTAS,IPAIR,KHS)                              03520043
      IF(ABS(QVALUE).LT.0.1) GOTO 270                                   03530000
      QC=QVALUE                                                         03540000
 270  CONTINUE                                                          03550000
      QEFF=QC                                                           03560028
      IF(ISHELL.EQ.0)QEFF=QC+SHELLS(1,1)                                03570028
      IF(ISHELL.EQ.0)SHELLS(1,1)=0.                                     03580028
      WRITE(6,269) QC,QEFF                                              03590028
 269  FORMAT(/' benutzter Q-Wert',F10.3,1X,'MeV',4X,
     1'effektiver Q-Wert',F10.3,1X,'MeV')                               03600028
C                                                                       03650000
C***********************************************************************03660000
C                                                  ELAB LOOP            03670000
 400  CONTINUE                                                          03680000
      IRC=-10                                                           03690021
      IEOF=0                                                            03700021
      CALL MYRD(ARRAY,IRC,23,5)                                         03710055
      IF(IRC.EQ.99) IEOF=1                                              03720021
      IF(IEOF.EQ.1) GOTO 3000                                           03730021
      IF(ARRAY(1).LE.0.) GOTO 3000                                      03740000
      EQ    =ARRAY(1)                                                   03750000
      IEXC  =ARRAY(2 ) +EPS                                             03760000
      IFUS  =ARRAY(3 ) +EPS                                             03770000
      LIMBAR=ARRAY(4 ) +EPS                                             03780000
      JLOWER=ARRAY(5 ) +EPS                                             03790000
      JUPPER=ARRAY(6 ) +EPS                                             03800000
      NEWFIS=ARRAY(7 ) +EPS                                             03810000
      ITSTRT=ARRAY(8 ) +EPS                                             03820000
      JFIS  =ARRAY(9 ) +EPS                                             03830000
      EQ1   =ARRAY(10)                                                  03840000
      IFUS1=0                                                           03850009
      IF(IFUS.EQ.10) IFUS1=1                                            03860009
      IF(IFUS.EQ.10) IFUS=8                                             03870009
      IF(IFUS.EQ.11) IFUS1=2                                            03880026
      IF(IFUS.EQ.11) IFUS=8                                             03890026
      IF(LPRINT.LT.3)WRITE(6,403) EQ,EQ1                                03900000
 403  FORMAT(1H1,4X,' ENERGY',2F10.3)                                   03910007
      IF(EQ1.GT.0.) IFUS=8                                              03920000
      IF(EQ1.LT.EQ) GOTO 404                                            03930000
      H=EQ1                                                             03940000
      EQ1=EQ                                                            03950000
      EQ =H                                                             03960000
 404  JLOWER=MAX0(JLOWER,1)                                             03970007
      JUPPER=MIN0(JUPPER,JDIM)                                          03980000
      IF(JFIS.GT.0 .AND. NOMORE.EQ.0) WRITE(19,100)(TITLE(I),I=1,18)    03990037
      NOMORE=1                                                          04000037
      IF(NEWFIS.NE.1) GOTO 407                                          04010000
      CALL MYRD(ARRAY,2,23,5)                                           04020055
      CLD1=ARRAY(1)                                                     04030000
      BFC=ARRAY(2)                                                      04040000
      IF(BFC .GT.0.) BARFAC=BFC                                         04050000
      IF(CLD1.GT.0.)CLD=CLD1                                            04060000
 407  DO 408 I=1,IPOPS                                                  04070000
      POPN(I)=0.                                                        04080000
      POPP(I)=0.                                                        04090000
      POPA(I)=0.                                                        04100000
 408  POP(I)=0.                                                         04110000
      DO 410 ISO=1,NUMISO                                               04120000
 410  SIGISO(ISO)=0.                                                    04130000
      DO 411 I=1,JDIM                                                   04140000
 411  FLAND(I)=0.                                                       04150000
      DO 412 NROW=1,NZUP                                                04160005
      DO 412 NSTEP=1,NAUP                                               04170005
 412  XSECT(NROW,NSTEP)=0.                                              04180005
      DO 413 J=1,JDIM                                                   04190036
      PER(J)=1.                                                         04200041
      PFIS(J)=0.                                                        04210041
 413  YFISS(J)=0.                                                       04220036
      DO 4131 K=1,KEN                                                   04230066
 4131 SPECNT(K)=0.                                                      04240067
      DO 4132 K=1,KEP1                                                  04250066
 4132 SPECPT(K)=0.                                                      04260067
      DO 4133 K=1,KEA1                                                  04270066
 4133 SPECAT(K)=0.                                                      04280067
      DO 4134 K=1,KEF                                                   04290066
 4134 SPECFT(K)=0.                                                      04300067
      DO 4135 K=1,KEG                                                   04310066
 4135 SPECGT(K)=0.                                                      04320067
      DO 4136 K=1,KEG                                                   04330066
 4136 SPECQT(K)=0.                                                      04340067
C--------------------------------------------------------------------   04350000
C             GET COMPOUND NUCLEUS POPULATION                           04360000
C                                                                       04370000
      AP=MPROJ                                                          04380000
      AT=MTGT                                                           04390000
      ZT=IZT                                                            04400000
      ZP=IZP                                                            04410000
      CM=AT/(AP+AT)                                                     04420000
      IF(IEXC.EQ.1) GOTO 414                                            04430000
      ELAB=EQ                                                           04440000
      ELAB1=EQ1                                                         04450000
      EXCIT =ELAB*CM+QC                                                 04460000
      EXCIT1=ELAB1*CM+QC                                                04470000
      IF(EQ1.LE.0.)EXCIT1=0.                                            04480000
      GOTO 415                                                          04490000
 414  EXCIT=EQ                                                          04500000
      EXCIT1=EQ1                                                        04510000
      ELAB=(EXCIT-QC)/CM                                                04520000
      ELAB1=(EXCIT1-QC)/CM                                              04530000
      IF(EQ1.LE.0.)ELAB1=0.                                             04540000
C415  EXMAX2=EXCIT                                                      04550029
 415  H4=0.04783*AP*CM*CM*ELAB                                          04560029
      C6=31.42/H4                                                       04570000
C--------------------------------------                                 04580000
      IF(IFUS.EQ.9 ) GOTO 420                                           04590000
      IF(IFUS.EQ.1 .OR. IFUS.EQ.5) GOTO 417                             04600000
      IF(IFUS.EQ.8) GOTO 430                                            04610000
      IF(AP.GT.4.) GOTO 416                                             04620000
      KPART=0                                                           04630000
      IDUM=AP+0.01                                                      04640000
      IDUM1=ZP+0.01                                                     04650000
      IF(IDUM.EQ.4 .AND. IDUM1.EQ.2)KPART=3                             04660000
      IF(IDUM.EQ.1 .AND. IDUM1.EQ.0)KPART=1                             04670000
      IF(IDUM.EQ.1 .AND. IDUM1.EQ.1)KPART=2                             04680000
      IF(KPART.EQ.0) GOTO 416                                           04690000
      LMAX=15                                                           04700000
      IF(KPART.EQ.3) LMAX=31                                            04710000
      DO 418 L=1,50                                                     04720000
 418  POPN(L)=0.                                                        04730000
      IPOTS=0                                                           04740000
      IF(IFUS.EQ.7) IPOTS=1                                             04750000
      IOMPR=1                                                           04760000
      IF(LPRINT.GT.3) IOMPR=0                                           04770000
      CALL OM(KPART,ELAB,LMAX,AT,ZT,SIGABS,POPN(1),IPOTS,IOMPR)         04780000
      FLGRAZ=SIGABS/C6                                                  04790000
      FLGRAZ=SQRT(FLGRAZ)-0.5                                           04800000
      FLGRAZ=DMAX1(AMPAR1,FLGRAZ)                                       04810000
      GOTO 419                                                          04820000
 416  CALL PARAP(ELAB,AP,AT,ZP,ZT,FINERT,RATIOS,IFUS,LIMBAR,0,1,IDUM,   04830000
     1           POPN(1),POPN(401),FLGRAZ)                              04840000
 419  IF(LIMBAR.EQ.1) GOTO 424                                          04850000
 417  CALL FUSE (ELAB,AP,AT,ZP,ZT,POPN(1),JDIM,IFUS,FLGRAZ)             04860070
      GOTO 424                                                          04870000
C-------------------------------------------------- IFUS=9              04880000
 420  IRC=JUPPER-JLOWER+1                                               04890000
      CALL MYRD(POPN(JLOWER),IRC,23,5)                                  04900055
 421  FORMAT(9F8.3)                                                     04910000
C--------------------------------------------------                     04920000
 424  JDIM1=JDIM+1                                                      04930000
      DO 425 J=1,JDIM                                                   04940000
      JJ=JDIM1 -J                                                       04950000
      IF(POPN(JJ).GT.0.) GOTO 426                                       04960000
 425  CONTINUE                                                          04970000
 426  JMAX=JJ                                                           04980000
      IF(JUPPER.EQ.0)JUPPER=JMAX                                        04990000
      SIGFUS=0.                                                         05000000
      DO 427 J=1,JDIM                                                   05010000
      IF(J.LT.JLOWER)POPN(J)=0.                                         05020000
      IF(J.GT.JUPPER)POPN(J)=0.                                         05030000
      SIGFUS=SIGFUS+POPN(J)                                             05040000
 427  CONTINUE                                                          05050000
      DO 428 J=JDIM1,IPOPS                                              05060000
 428  POPN(J)=0.                                                        05070000
      M2=1                                                              05080000
      GOTO 460                                                          05090000
C-------------------------------------- IFUS=8                          05100000
 430  IF(IFUS1.EQ.0) CALL MYRD(ARRAY,9,23,5)                            05110055
      IF(IFUS1.EQ.1) CALL MYRD(ARRAY,10,23,5)                           05120055
      IF(IFUS1.EQ.2) CALL MYRD(ARRAY,14,23,5)                           05130055
      V0=ARRAY(1)                                                       05140000
      R0=ARRAY(2)                                                       05150000
      D=ARRAY(3)                                                        05160000
      Q2   =ARRAY(4)                                                    05170000
      IF(Q2.LT.1.) Q2=1.09*ZT*(AT**0.666667)*Q2*(1.+0.15577*Q2)         05180000
C     IF Q2 IS LESS THAN ONE IT IS INTERPRETED AS DEFO BETA             05190000
      CRED=ARRAY(5)                                                     05200000
      NOCURV=ARRAY(6)+0.01                                              05210000
      NOPROX=ARRAY(7)+0.01                                              05220000
      IOPT1=ABS(ARRAY(8))+0.01                                          05230000
      IF(ARRAY(8).LT.0.) IOPT1=-IOPT1                                   05240000
      IF(IOPT1.EQ.5) D=1.                                               05250032
      ITEST1=ABS(ARRAY(9))+0.01                                         05260000
      IF(ARRAY(9).LT.0.) ITEST1=-ITEST1                                 05270000
      IF(IFUS1.GT.0) SIGR0=ARRAY(10)                                    05280026
      CUTOF=2.5                                                         05290026
      XTH=0.                                                            05300026
      IF(IFUS1.LT.2) GOTO 436                                           05310026
      CUTOF=ARRAY(11)                                                   05320027
      IF(CUTOF.LE.0.) CUTOF=2.5                                         05330026
      DUM=CUTOF*SIGR0                                                   05340026
      IF(DUM.LT.50.) GOTO 434                                           05350026
      WRITE(6,432) SIGR0,CUTOF                                          05360026
 432  FORMAT(' UNREASONABLE PARAMETERS SIGR0 AND/OR CUTOF',2F8.3)       05370027
      GOTO 9991                                                         05380026
 434  XTH=ARRAY(12)                                                     05390026
      APUSH=ARRAY(13)                                                   05400026
      IF(APUSH.LE.0.) APUSH=12.                                         05410026
      FPUSH=ARRAY(14)                                                   05420026
      SPUSH=ARRAY(15)
      IF(FPUSH.LE.0.) FPUSH=0.75                                        05430026
 436  IF(CRED.LE.0.)CRED=1.                                             05440026
      IF(LIMBAR.EQ.1) LMAX=-1                                           05450000
C---------------------------------                                      05460000
      ELB=ELAB                                                          05470000
      IENG=1                                                            05480000
      M2=1                                                              05490000
      RNORM=1.                                                          05500000
      SIGFUS=0.                                                         05510000
      FLGRAZ=0.                                                         05520000
      IESTEP=1                                                          05530000
      JMAX=0                                                            05540000
      IF(EQ1.LE.0.) GOTO 444                                            05550000
      DUM=EXCIT-EXCIT1                                                  05560000
      M2=DUM+1.5                                                        05570000
 442  IDUM=M2/IESTEP                                                    05580000
      IF(IDUM.LT.20) GOTO 443                                           05590000
      IESTEP=IESTEP+1                                                   05600000
      GOTO 442                                                          05610000
 443  DELB=IESTEP/CM                                                    05620000
      RNORM=M2                                                          05630000
C                                                                       05640000
 444  DO 452 IE=1,M2,IESTEP                                             05650000
      LMAX=0                                                            05660000
      IF(LIMBAR.EQ.1)LMAX=-1                                            05670000
      FLGRZ=0.                                                          05680000
      R01=R0                                                            05690024
      IF(IOPT1.EQ.5) GOTO 4446                                          05700032
      IF(R01.GT.0.) GOTO 4448                                           05710025
 4446 AP3=AP**0.333333                                                  05720032
      AT3=AT**0.333333                                                  05730025
      RPS=1.28*AP3-0.76+0.8/AP3                                         05740025
      RTS=1.28*AT3-0.76+0.8/AT3                                         05750025
      IF(IOPT1.EQ.1 .OR. IOPT1.EQ.2 .OR.IOPT1.GT.3) GOTO 4445           05760032
      RPS=RPS-1./RPS                                                    05770025
      RTS=RTS-1./RTS                                                    05780025
 4445 R01=(RTS+RPS)/(AP3+AT3)                                           05790025
 4448 MPT=20                                                            05800031
C
C     I.Giese 2007
C     Hilfsvaribale statt Konstante 1, da Werte zurueckgegeben werden
      K_NORUTH=1
      K_ION=1
C
C     I.Giese 2007
C     Bemerkung: SIGFU nicht besetzt
      CALL FUSIO (ELB ,AP,ZP,AT,ZT,Q2,V0,R01,D,XTH,APUSH,FPUSH,         05810026
     1            SIGR0,CRED,K_NORUTH,NOPROX,NOCURV,                    05820025
     2        K_ION,POP ,  FLGRZ,LMAX,SIGFU ,ITEST1,IOPT1,MPT,CUTOF)    05830025
      LMAX=MIN0(LMAX,JDIM)                                              05840000
      JMAX=MAX0(JMAX,LMAX)                                              05850000
      JMAX1=LMAX                                                        05860000
      JMAX1=MIN0(JMAX1,JDIM)                                            05870000
      ICHECK=JMAX1*M2                                                   05880000
      IF(ICHECK.LE.IPOPS)GOTO 445                                       05890000
      WRITE(6,446) M2,JMAX1,IPOPS                                       05900000
 446  FORMAT(' DIMENSIONS TOO LARGE  M2*JMAX1>IPOPS',3I6)               05910000
      GOTO 9991                                                         05920000
 445  IF(JUPPER.GT.0)JMAX1=MIN0(LMAX,JUPPER)                            05930000
      FLGRAZ=DMAX1(FLGRAZ,FLGRZ)                                        05940000
      DO 448 J=JLOWER,JMAX1                                             05950000
      DO 447 IEE=1,IESTEP                                               05960000
      INDEX=IE-1+IEE                                                    05970000
      IF(INDEX.GT.M2) GOTO 447                                          05980000
      INDEX=INDEX+(J-1)*M2                                              05990000
      POPN(INDEX)=POP(J)/RNORM                                          06000000
      SIGFUS=SIGFUS+POPN(INDEX)                                         06010000
 447  CONTINUE                                                          06020000
 448  CONTINUE                                                          06030000
      ELB=ELB-DELB                                                      06040000
 452  CONTINUE                                                          06050000
C                                                                       06060000
      IF(EQ1.LE.0.) GOTO 450                                            06070000
      WRITE(6, 449) EXCIT,EXCIT1,M2,IESTEP,RNORM,DELB,SIGFUS            06080000
 449  FORMAT(' CROSS SECTION AVERAGED',2F6.1,2I6,3F8.3/)                06090000
 450  IF(JUPPER.EQ.0)JUPPER=JMAX                                        06100000
      DO 453 J=1,IPOPS                                                  06110000
 453  POP(J)=0.                                                         06120000
C---------------------------------                                      06130000
 460  IF(SIGFUS.GT.0.1E-5)GOTO 463                                      06140000
      WRITE(6,461) ELAB,EXCIT                                           06150000
 461  FORMAT(//'E_lab =',F8.3,' MeV',4X,'E* =',F8.3,' MeV',3X,
     1'nicht gerechnet, SIGFUS kleiner als 1 nb'//)                     06160000
      GO TO 400                                                         06180000
C---------------------------------                                      06190000
 463  N2=JUPPER                                                         06200000
      JUPYR=JUPPER+6                                                    06210000
      JUPYR=MIN0(JUPYR,JDIM-1)                                          06220000
      SIGMXZ =0.                                                        06230000
      SIGEVA=0.                                                         06240000
      DO 464 NROW=1,NZUP                                                06250000
 464  SIGZ(NROW)=0.                                                     06260000
      SIGFIS=0.                                                         06270000
 462  IU1=1-MOD(NZ1,2)                                                  06280000
      IU2=0                                                             06290000
      IDUM=(NZ1+1)/2                                                    06300000
      IF(MOD(IDUM,2).NE.0) IU2=1                                        06310000
      DO 465 I=1,NAUP2                                                  06320000
      DO 465 K=1,NZUP2                                                  06330000
 465  STORE(I,K)=0.                                                     06340000
      REWIND 10                                                         06350000
      REWIND 11                                                         06360000
      REWIND 12                                                         06370000
      REWIND 13                                                         06380000
      IF(NZ1.EQ.1) GOTO 469                                             06390000
 1001 FORMAT(18I4)                                                      06470000
 1002 FORMAT(3F8.2,4E12.4)                                              06480000
 1003 FORMAT(9F6.1)                                                     06490000
 1004 FORMAT(6E12.4)                                                    06500000
      NZ1=NROW+1                                                        06510000
      IF(MOD(NZ1,2).EQ.0) GOTO 469                                      06520000
      IUNIT=13-IU2                                                      06530000
      IREAD=1                                                           06540000
      DO 467 NSTEP=1,NA                                                 06550000
      IF(STORE(NZ1+1,NSTEP).LT.2) GOTO 467                              06560000
      CALL GETPUT(POP,IPOPS,IUNIT,IREAD,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)06570000
 467  CONTINUE                                                          06580000
 469  CONTINUE                                                          06590000
C     WRITE(6,860) IPAIR,CST,DELT,DELTA                                 06660000
C                                                                       06670079
 468  IDUM=EXCIT+1.                                                     06680079
      IF(IDUM.LT.IEDIM) GOTO 474                                        06690000
      WRITE(6,472) IEDIM                                                06700000
      GOTO 400                                                          06710000
 472  FORMAT(' die hoechste erlaubte Anregungsenergie ist ',I5,' MEV')  06720000
 474  CONTINUE                                                          06730000
C                                                                       06740000
C******************************* ROW LOOP (NZ) *************************06750000
      EXCIT0=EXCIT                                                      06760029
      IF(ISHELL.EQ.0) EXCIT=EXCIT-QC+QEFF                               06770029
      EXMAX2=EXCIT                                                      06780029
      DO 1000 NROW=NZ1,NZ                                               06790000
C                                                                       06800000
      IU1=1-IU1                                                         06810000
      IF(MOD(NROW,2).EQ.0)IU2=1-IU2                                     06820000
      KZ(1)=IZCN+1-NROW                                                 06830000
      ZPAR=KZ(1)                                                        06840000
      APAR =IACN+1-NROW                                                 06850000
      IN(1)=APAR-ZPAR                                                   06860000
      DO 476 NSTEP=2,NA                                                 06870000
      KZ(NSTEP)=KZ(1)                                                   06880000
 476  IN(NSTEP)=IN(1)+1-NSTEP                                           06890000
C                                                                       06900000
 480  DO 490 NSTEP=1,NAUP                                               06910000
      TOTAL(NSTEP)=0.                                                   06920000
      SUMEG(NSTEP)=0.                                                   06930004
      SUMJG(NSTEP)=0.                                                   06940006
      AVGJ(NSTEP)=0.                                                    06950000
      AVGE(NSTEP)=0.                                                    06960000
      AVGEF(NSTEP)=0.                                                   06970000
      AVGJF(NSTEP)=0.                                                   06980000
      SIG  (NSTEP)=0.                                                   06990000
      SIGNN(NSTEP)=0.                                                   07000000
      SIGP(NSTEP)=0.                                                    07010000
      SIGA(NSTEP)=0.                                                    07020000
      TRIM(NSTEP)=0.                                                    07030000
      RP (NSTEP)=0.                                                     07040000
      RLOSTG(NSTEP)=0.                                                  07050000
      FLOST (NSTEP)=0.                                                  07060000
      RLOSTN(NSTEP)=0.                                                  07070000
      RLOSTP(NSTEP)=0.                                                  07080000
      RLOSTA(NSTEP)=0.                                                  07090000
      RLOSTF(NSTEP)=0.                                                  07100000
      IND(NSTEP)=0                                                      07110000
      ISEQ(NSTEP)=0                                                     07120000
 490  CONTINUE                                                          07130000
      DO 492 NSTEP = 1,NA                                               07140000
 492  ISEQ(NSTEP)=1                                                     07150000
      DO 493 K=1,KEN                                                    07160000
 493  SPECN(K)=0.                                                       07170000
      DO 494 K=1,KEP1                                                   07180000
 494  SPECP(K)=0.                                                       07190000
      DO 495 K=1,KEA1                                                   07200000
 495  SPECA(K)=0.                                                       07210000
      DO 496 K=1,KEF                                                    07220000
 496  SPECF(K)=0.                                                       07230000
      DO 497 K=1,KEG                                                    07240000
 497  SPECG(K)=0.                                                       07250000
      DO 498 K=1,KEG                                                    07260000
 498  SPECQ(K)=0.                                                       07270000
      IF(NROW.EQ.NZ1) GOTO 510                                          07280000
      DO 499 K=1,IPOPS                                                  07290000
      POP (K)=0.                                                        07300000
      POPN(K)=0.                                                        07310000
      POPP(K)=0.                                                        07320000
      POPA(K)=0.                                                        07330000
      PENTRY(K)=0.                                                      07340000
 499  CONTINUE                                                          07350000
 510  SIGMAX=0.                                                         07360000
C                                                                       07400000
C*********************************************************************  07410000
C                                    SEQUENCE LOOP (NA)                 07420000
      DO 800 NSTEP = 1,NA                                               07430000
      DO 630 ISO=1,NUMISO                                               07440000
      IF(ISOA(ISO).EQ.NSTEP .AND. ISOZ(ISO).EQ.NROW) GOTO 631           07450000
 630  CONTINUE                                                          07460000
      ISO=0                                                             07470000
 631  CONTINUE                                                          07480000
C                                                                       07490000
      APAR=ZPAR+IN(NSTEP)                                               07500000
      IF(NSTEP.EQ.1) GOTO 635                                           07510000
      IF(MOD(NSTEP,3).NE.0 .OR. ITRANS.EQ.0) GOTO 438                   07520000
 635                            IUNIT=21                                07530000
      IF(LPRINT.LT.3)                                                   07540000
     1WRITE(6,636) APAR,ZPAR,SIGLOW,IUNIT,IOWKB,KEN,KEP,KEA,LN,LP,LA,DEL07550000
 636  FORMAT(2F6.0,F9.4,8I6/' DEL',5F8.1)                               07560000
             CALL TRANSM(APAR,ZPAR,SIGLOW,TLOW,IPTR,IUNIT,IOWKB,QT)     07570003
      IF(LPRINT.LT.3)                                                   07580000
     1WRITE(6,475) APAR,ZPAR,SIGLOW,IUNIT,DEL(2),DEL(3),DEL(4)          07590000
 475  FORMAT(' TRANSM  APAR,ZPAR,SIGLOW,UNIT,DEL',2F6.1,E10.2,I4,3F8.1) 07600000
C                         TRANSM CALLS OVER2(WHICH CALLS TLD) OR OWKB   07610083
C                                                                       07620000
 438                                   CALL SEQIN(NROW,NSTEP,LDBM,      07630000
     1                                         CST,CLD)                 07640000
C                                                                       07650000
C                                                                       07660000
 640  APAR=ZPAR+IN(NSTEP)                                               07670000
      DO 641 K=1,5                                                      07680000
 641  SHEL(K)=SHELK(K)                                                  07690000
                         CALL ROT (NSTEP,NROW,APAR,ZPAR,SHEL,ITSTRT)    07700000
C                                                                       07710000
 660  IF(NSTEP.GT.1)                                                    07720000
     1SIGMAX=DMAX1(SIGMAX,SIG(NSTEP-1))                                 07730000
      DUM2=FRACT2*SIGMAX                                                07740000
      SUMLOW=DMAX1(ABSMIN,DUM2)                                         07750000
                                       CALL GRIDS(NROW,NSTEP,CUTPOP)    07760000
C                                  GRIDS CALLS SUMPOP,GETPUT,CUTOFF,TOT 07770000
      IF(ISEQ(NSTEP).EQ.0) GO TO 700                                    07780000
                                       CALL DENSTY(NROW,NSTEP,LDBM,LEP, 07790000
     1                                                        LEVPRN)   07800000
C                                                                       07810000
      IF(NOLEP.EQ.1)IAMX=10000                                          07820000
                                       CALL EVA(NROW,NSTEP,CUTPOP,NOLEP,07840000
     1                                    NOLJI,NOLJF,LEVPRN,LDBM)      07850000
      IF(LPRINT.LT.3)                                                   07870000
     1WRITE(6,1006) NSTEP                                               07880000
 1006 FORMAT( ' EVA-,STEP',I3)                                          07890000
C                                      EVA CALLS POPUL                  07900000
C                                                                       07910000
                                       CALL SEQOUT(NROW,NSTEP,NOPUT)    07920000
                                       GOTO 800                         07990000
 700  IF(NSTEP.EQ.NA) GOTO 800                                          08000000
      NSTEP1=NSTEP+1                                                    08010000
      DO 710 NST=NSTEP1,NA                                              08020000
      IF(STORE(NROW,NST).GT.0.) GOTO 800                                08030000
 710  CONTINUE                                                          08040000
      GOTO 850                                                          08050000
 800  CONTINUE                                                          08060000
C                                                                       08070000
C******************** END SEQUENCE LOOP(NA) ****************************08080000
C                                                                       08090000
 850  CALL ROWOUT(NROW,LOGUN,SIGMXZ,NOPUT,EXCIT0)                       08100030
C 860  FORMAT(' IPAIR,CST,DELT,DELTA',I4,7F8.3)                         08110000
      IF(NOPUT.EQ.1) GOTO 990                                           08120000
      DUM=0.                                                            08130000
      IDUM=0                                                            08140000
      CALL GETPUT(POP,IPOPS,10,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)    08150000
      REWIND 10                                                         08160000
      CALL GETPUT(POP,IPOPS,11,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)    08170000
      REWIND 11                                                         08180000
      IF(MOD(NROW,2).EQ.0)GOTO 990                                      08190000
      CALL GETPUT(POP,IPOPS,12,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)    08200000
      REWIND 12                                                         08210000
      CALL GETPUT(POP,IPOPS,13,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)    08220000
      REWIND 13                                                         08230000
 990  SIGMXZ=DMAX1(SIGMXZ,SIGZ(NROW))                                   08240000
      DUM=FRACT2*SIGMXZ                                                 08250000
      IF(SIGZ(NROW).LT.DUM) GOTO 1500                                   08260000
 1000 CONTINUE                                                          08270000
C                                     END ROW LOOP                      08280000
C***********************************************************************08290000
 1500 IF(LPRINT .GT. 4) GOTO 1510                                       08300000
      WRITE(6,403)EQ,EQ1                                                08310000
 1510 CALL ENGOUT(LOGUN,IDISC,EXCIT0)                                   08320030
 2000 CONTINUE                                                          08330000
      GOTO 400                                                          08340000
C                                     END ELAB LOOP                     08350000
C***********************************************************************08360000
 3000 CALL CASOUT(LOGUN,IDISC)                                          08370030
      IF(IEOF.EQ.1) GOTO 9999                                           08380021
      GOTO 99                                                           08390000
C***********************************************************************08400000
C                                       EXIT                            08410000
 9999 WRITE(6,9992)                                                     08440000
 9992 FORMAT(///' Ende HIVAP ')                                         08450000
 9993 FORMAT('C ',18A4)                                                 08460000
 9991 CLOSE(10)                                                         08470076
      CLOSE(11)                                                         08490076
      CLOSE(12)                                                         08510076
      CLOSE(13)                                                         08530076
      CLOSE(22)                                                         08550076
      CLOSE(23)                                                         08570076
      CLOSE(5)
      CLOSE(6)
      CLOSE(9)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(19)
      CLOSE(21)
      CLOSE(25) 
      STOP                                                              08590076
      END                                                               08600000
      BLOCK DATA                                                        08610000
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     I.Giese 2007
C     Wegen vlro.f mit 2 Entries
      COMMON/TESTS/L0_ENTRY,L1_ENTRY,L2_ENTRY
C
      COMMON/ELEM/ELEMNT/S/BLANK,NU,SYMB                                08620000
      CHARACTER*4 ELEMNT(128)/                                          08630050
     1 ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',               08640051
     2 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',               08650051
     3 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',               08660051
     4 'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',               08670051
     5 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',               08680051
     6 'Sb','Te',' J','Xe','Cs','Ba','La','Ce','Pr','Nd',               08690051
     7 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',               08700051
     8 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',               08710051
     9 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',               08720051
     X 'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',               08730051
     Y 'Md','No','Lr','Rf','Ha','Sg','Ns','Hs','Mt','110/',             08740051
     Z '111/','112/','113/','114/','115/','116/','117/','118/','119/',  08750051
     A '120/','121/','122/','123/','124/','125/','126/','127/','128/'/  08760051
C                                                                       08770000
      CHARACTER*1 BLANK/' '/,       NU(15)/'1','2','3','4','5','6','7', 08780049
     1 '8','9','A','B','C','D','E','F'/,SYMB(7)/'G','N','P','A','F',    08790000
     2 '*','X'/                                                         08800000
C
C     I.Giese 2007
      DATA L0_ENTRY/0/,L1_ENTRY/0/,L2_ENTRY/0/
C
      END                                                               08810000
