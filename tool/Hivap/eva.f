      SUBROUTINE EVA(NROW,NSTEP,CUTPOP,NOLEP,NOLJI,NOLJF,LEVPRN,LDBM)   00020000
C
C  Berechnet ED-JD - Bevoelkerung der Tochter aus Mutter - ED-JD -
C  Bevoelkerung, ruft POPUL auf
C  Option zur Behinderung der Spaltung
C  Stand: 2.5.1994
C 
      INCLUDE 'common.f'                                                00030010
      COMMON/DELFIS/EDELFIS,ADELFIS,INOF
      CHARACTER*1 DOT/'.'/, OU/'O'/                                     00040010
      INTEGER*2 NULL                                                    00050000
      INTEGER*2 TABL(4096)                                              00060010
C-----------------------------------------------------------------------00070000
C     MAYBE UNINITIALIZED: JOSC1,IADD1,IADD2,IADD3,DUM,LMAX2,LMIN1
C
C     I.Giese 2007
      SAVE
C
      BELOW=0.                                                          00080008
      NULL=0                                                            00090000
      DO 120 I=1,IPOPS                                                  00100000
 120  TABLE(I)=BLANK                                                    00110000
      HCST=41.35                                                        00120000
C     PLANCKS CONSTANT (MEVSEC*10E22)                                   00130000
      DO 150 K=1,KEG                                                    00140000
      IF(EFIVE(KEG-K).GT.0.) GOTO 151                                   00150000
 150  CONTINUE                                                          00160000
      KE2MX=0                                                           00170000
      GOTO 152                                                          00180000
 151  KE2MX=K                                                           00190000
C--------------------------------                                       00200000
 152  IF(LPRINT.LT.3)                                                   00210000
     1WRITE(6,100) M1,N1,NOP,NOA,NOF,NOG,KE2MX,NOLEP,NOLJI,NOLJF,JFJI   00220000
 100  FORMAT(' EVA M1,N1=',2I6/' NOP,NOA,NOF,NOG,KE2MX=',5I2,           00230000
     1       '  NOLEP,NOLJI,NOLJF,JFJI=',4I2)                           00240000
C--------------------------------                                       00250000
      EWG=FLOAT(MWG)+EXMING                                             00260000
      EWN=FLOAT(MWN)+EXMINN                                             00270000
      EWP=FLOAT(MWP)+EXMINP                                             00280000
      EWA=FLOAT(MWA)+EXMINA                                             00290000
      EWF=FLOAT(MWF)+EXMINF                                             00300000
      MG=EXMAX1+2.-DEL(1)                                               00310000
      MN=EXMAX1+Q(2)+2.-DEL(2)                                          00320000
      MP=EXMAX1+Q(3)+2 -DEL(3)                                          00330000
      MA=EXMAX1+Q(4)+2.-DEL(4)                                          00340000
      MF=EXMAX1+2.-DEL(5)                                               00350000
      ICU=0                                                             00360000
      IOFF=0                                                            00370000
      CUTS=0.                                                           00380000
C--------------------                                                   00390000
      IF(JOSC1.NE.1)JOSC1=0                                             00400000
      IF(JOSC .NE.1)JOSC =0                                             00410000
      IF(IADD1.NE.1)IADD1=0                                             00420000
      IF(IADD2.NE.1)IADD2=0                                             00430000
      IF(IADD3.NE.1)IADD3=0                                             00440000
C---------------------------------LEP LJI                               00450000
      LEP=1                                                             00460000
      IF(NOLEP.EQ.0 .AND. IAMX.LT.2000) LEP=2                           00470000
      LJI=1                                                             00480000
      IF(NOLJI.EQ.0 .AND. AVGJ(NSTEP).GT.10.)LJI=2                      00490000
C                                                                       00500000
      IF(LEP.EQ.1 .AND. LJI.EQ.1) GOTO 195                              00510000
      IOSCE=MOD(NSTEP,2 )                                               00520000
      IF(LEP.EQ.1) IOSCE=1                                              00530000
      IF(LJI.EQ.1) JOSC=1                                               00540000
C                                                                       00550000
      DO 180 IEP=1,M1                                                   00560000
      IF(LEP.EQ.2) IOSCE=1-IOSCE                                        00570000
      E=IEP-1                                                           00580000
      EXCURR=EXMAX1-E                                                   00590000
      JPEP=IEP                                                          00600000
      IF(IOSCE.EQ.1)JOSC1=1-JOSC1                                       00610000
      JOSC=JOSC1                                                        00620000
      DO 170 JP1=1,N1                                                   00630000
      IF(LJI.EQ.2)JOSC=1-JOSC                                           00640000
      IEF =EXCURR-EJAY(JP1)                                             00650000
      TABL(JPEP)=1                                                      00660010
      IF(((IOSCE.EQ.0) .OR. (JOSC.EQ.0)).AND. (IEF .GT.20))             00670000
     1   TABL(JPEP)=NULL                                                00680010
      IF(POP(JPEP).LE.0.) TABL(JPEP)=NULL                               00690010
      JPEP=JPEP+M1                                                      00700000
 170  CONTINUE                                                          00710000
 172  FORMAT(1X,120I1)                                                  00720000
 180  CONTINUE                                                          00730000
C---------------------REPOP                                             00740000
      M1N1=M1*N1                                                        00750000
      DO 192 IEP=1,M1                                                   00760000
      JPEP=IEP-M1                                                       00770000
      DO 190 JP1=1,N1                                                   00780000
      JPEP=JPEP+M1                                                      00790000
      IF(TABL(JPEP).GT.NULL) GOTO 190                                   00800010
      IF(POP(JPEP).LE.0.)GOTO 190                                       00810000
      ISU2=0                                                            00820000
      L=JPEP-M1                                                         00830000
      IF(L.GT.0)ISU2=ISU2+TABL(L)                                       00840010
      L=JPEP+M1                                                         00850000
      IF(L.LE.M1N1) ISU2=ISU2+TABL(L)                                   00860010
      IF(IEP.EQ.M1) GOTO 1810                                           00870000
      L=JPEP+1                                                          00880000
      IF(L.LE.M1N1) ISU2=ISU2+TABL(L)                                   00890010
      L=JPEP-M1+1                                                       00900000
      IF(L.GT.0)ISU2=ISU2+TABL(L)                                       00910010
      L=JPEP+M1+1                                                       00920000
      IF(L.LE.M1N1)ISU2=ISU2+TABL(L)                                    00930010
 1810 IF(IEP.EQ.1) GOTO 183                                             00940000
      L=JPEP-1                                                          00950000
      IF(L.GT.0)ISU2=ISU2+TABL(L)                                       00960010
      L=JPEP-M1-1                                                       00970000
      IF(L.GT.0)ISU2=ISU2+TABL(L)                                       00980010
      L=JPEP+M1-1                                                       00990000
      IF(L.LE.M1N1)ISU2=ISU2+TABL(L)                                    01000010
C                                                                       01010000
 183  ISUM=ISU2                                                         01020000
      IF(ISUM.EQ.0)GOTO  190                                            01030000
      DUM=POP(JPEP)/FLOAT(ISUM)                                         01040000
      POP(JPEP)=0.                                                      01050000
 1830 FORMAT(5I6,E12.4)                                                 01060000
C                                                                       01070000
      L=JPEP-M1                                                         01080000
      IF(L.LE.0)GOTO 1840                                               01090000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01100010
 1840 L=JPEP+M1                                                         01110000
      IF(L.GT.M1N1) GOTO 185                                            01120000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01130010
 185  IF(IEP.EQ.M1) GOTO 1841                                           01140000
      L=JPEP+1                                                          01150000
      IF(L.GT.M1N1) GOTO 1842                                           01160000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01170010
 1842 L=JPEP-M1+1                                                       01180000
      IF(L.LE.0) GOTO 1844                                              01190000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01200010
 1844 L=JPEP+M1+1                                                       01210000
      IF(L.GT.M1N1) GOTO 1841                                           01220000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01230010
 1841 IF(IEP.EQ.1) GOTO 190                                             01240000
      L=JPEP-1                                                          01250000
      IF(L.LE.0) GOTO 1843                                              01260000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01270010
 1843 L=JPEP-M1-1                                                       01280000
      IF(L.LE.0)GOTO 1845                                               01290000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01300010
 1845 L=JPEP+M1-1                                                       01310000
      IF(L.GT.M1N1) GOTO 190                                            01320000
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM                             01330010
 190  CONTINUE                                                          01340000
 192  CONTINUE                                                          01350000
 195  NOGG=NOG                                                          01360000
C------------------------------                                         01370000
C                                                                       01380000
C******************************  EP LOOP  ******************************01390000
C                                                                       01400000
      DO 2000 IEP=1    ,M1                                              01410000
C------------------------------                                         01420000
      E=FLOAT(IEP-1)+0.0001                                             01430000
      EXCURR=EXMAX1-E                                                   01440000
      ILAND=NROW+NSTEP+IEP                                              01450000
C------------------------------------------------------------           01460000
      NOG=NOGG                                                          01470000
      MAXG=MIN0(MG-IEP,KEG)                                             01480000
      IF(MAXG.LT.0)  MAXG=0                                             01490000
      IF(NOG    .EQ.1) MAXG=0                                           01500000
      IGL=MAXG+MOD(MAXG,2)                                              01510000
      IGL=MAX0(IGL,4)                                                   01520000
      ITG=5*IGL                                                         01530000
C                                                                       01540000
      MAXN=MIN0(MN-IEP,KEN)                                             01550000
      IF(MAXN.LT.0)  MAXN=0                                             01560000
      IF(NON.EQ.1) MAXN=0                                               01570000
      INL=MAXN+MOD(MAXN+1,2)                                            01580000
      INL=MAX0(INL,5)                                                   01590000
      ITN=(2*LN -1)*INL                                                 01600000
      IF(JFJI.EQ.1)ITN=INL                                              01610000
C------------------------------------------------------------           01620000
      MAXP=MIN0(MP-IEP,KEP)                                             01630000
      IF(MAXP.LT.0)  MAXP=0                                             01640000
      IF(NOP.EQ.1) MAXP=0                                               01650000
      IPL=MAXP+MOD(MAXP,2)                                              01660000
      IPL=MAX0(IPL,4)                                                   01670000
      ITP=(2*LP -1)*IPL                                                 01680000
      IF(JFJI.EQ.1)ITP=IPL                                              01690000
C------------------------------------------------------------           01700000
      MAXA=MIN0(MA-IEP,KEA)                                             01710000
      IF(MAXA.LT.0)  MAXA=0                                             01720000
      IF(NOA     .EQ.1) MAXA=0                                          01730000
      IAL=MAXA+MOD(MAXA,2)                                              01740000
      IAL=MAX0(IAL,4)                                                   01750000
      ITA=(2*LA -1)*IAL                                                 01760000
      IF(JFJI.EQ.1)ITA=IAL                                              01770000
C-------------------------------------------------------------          01780000
      MAXF=MIN0(MF-IEP,KEF)                                             01790000
      IF(EXCURR.LT.BFS(NSTEP,NROW))  MAXF=0                             01800000
      IF(NOF.EQ.1) MAXF=0                                               01810000
C--------------------------------------------------------------         01820000
C                                                                       01830000
C************ READJUST LEVEL DENSITY ***********************************01840000
C                                                                       01850000
  200 IF(IEP.EQ.1) GOTO 230                                             01860000
  201 IF(MWG.LT.IDEN/NWG.AND.MWA.LT.IDNA/NWA.AND.MWN.LT.IDEN/NWN.AND.   01870000
     1    MWP.LT.IDEN/NWP.AND.MWF.LT.IDEN/NWF) GOTO 230                 01880000
      IDUM=0                                                            01890000
      IF(EXMAXG.LT.EWG)  GOTO 205                                       01900000
      LADJ(1)=1                                                         01910000
      IDUM=1                                                            01920000
 205  IF(EXMAXN.LT.EWN)  GOTO 210                                       01930000
      LADJ(2)=1                                                         01940000
      IDUM=1                                                            01950000
 210  IF(EXMAXP.LT.EWP)  GOTO 215                                       01960000
      LADJ(3)=1                                                         01970000
      IDUM=1                                                            01980000
 215  IF(EXMAXA.LT.EWA)  GOTO 220                                       01990000
      LADJ(4)=1                                                         02000000
      IDUM=1                                                            02010000
 220  IF(EXMAXF.LT.EWF)  GOTO 225                                       02020000
      LADJ(5)=1                                                         02030000
      IDUM=1                                                            02040000
 225  IF(IDUM .EQ.0) GO TO 230                                          02050000
      CALL DENSTY(NROW,NSTEP,LDBM,1   ,LEVPRN)                          02060000
C     GO TO 200                                                         02070000
C***********************************************************************02080000
  230 JPEPMX=(N1-1)*M1+IEP                                              02090000
C                                    OFFSETS FOR LEVDNS                 02100000
      KONSTG=EXMAXG-EXMAX1+E+DEL(1)                                     02110000
      KONSTN=EXMAXN-EXMAX1+E-Q(2)+DEL(2)                                02120000
      KONSTP=EXMAXP-EXMAX1+E-Q(3)+DEL(3)                                02130000
      KONSTA=EXMAXA-EXMAX1+E-Q(4)+DEL(4)                                02140000
      KONSTF=EXMAXF-EXMAX1+E+DEL(5)                                     02150000
C                                                                       02160000
C******************* JP (JI,J) LOOP ************************************02170000
C                                                                       02180000
      DO 1000 JPEP= IEP,JPEPMX, M1                                      02190000
C--------------------------------------FIRST CHECK CUTOFF-------        02200000
      IF(POP(JPEP).GT.0.)GOTO 350                                       02210000
      ICU=ICU+1                                                         02220000
      TABLE(JPEP)=BLANK                                                 02230000
      GOTO 1000                                                         02240000
 350  JI=(JPEP-IEP)/M1                                                  02250002
      IF(POP(JPEP).GT.CUTPOP )GOTO 360                                  02260002
 355  IOFF=IOFF+1                                                       02270000
      CUTS=CUTS+POP(JPEP)                                               02280000
      PENTRY(JPEP)=0.                                                   02290003
      TABLE(JPEP)=DOT                                                   02300000
C     TOTAL(NSTEP)=TOTAL(NSTEP)-POP(JPEP)                               02310000
C     POP(JPEP)=0.                                                      02320000
      GOTO 1000                                                         02330000
C-----------------------------------------------------------------------02340000
 360  JI1=JI+1                                                          02350002
      ETH=EXCURR-EJAY(JI1)                                              02360008
      IF(ETH.GT.0.) GOTO 364                                            02370008
      TABLE(JPEP)=OU                                                    02380008
      BELOW=BELOW+POP(JPEP)                                             02390008
      POP(JPEP)=0.                                                      02400008
      GOTO 1000                                                         02410008
 364  IEF=ETH + 0.5                                                     02420008
C     IF(TABL(JPEP).EQ.NULL .AND. IEF.GT.20) GOTO 355                   02430010
C     ------------------------------ISOMER                              02440000
      IF(ISO.EQ.0) GOTO 380                                             02450000
      IF(JI.EQ.ISOJ(ISO) .AND. DUM.LT.1.1) SIGISO(ISO)=SIGISO(ISO)      02460000
     1                                                 +POP(JPEP)       02470000
C     ------------------------------                                    02480000
 380  IF(IOPT.EQ.0) GOTO 400                                            02490000
      IF(DUM                .GT.STRIPE) GOTO 400                        02500000
      P(1)=-1.                                                          02510000
      GOTO 900                                                          02520000
C-----------------------------------------------NEUTRONS-------         02530000
  400 IF(MAXN.EQ.0) GOTO 490                                            02540000
      IF(EXMAX2.LE.0.) GOTO 490                                         02550000
      DO 404 K=1,ITN                                                    02560000
  404 TMTBN (K)=0.                                                      02570000
      DO 408 K=1,INL                                                    02580000
  408 TSPECN(K)=0.                                                      02590000
      IF(JFJI.EQ.1) GOTO 407                                            02600000
      JFMAX1=MIN0(NWN,JI+LN )                                           02610000
      JFMIN1=MAX0(1,JI+2-LN )                                           02620000
      GOTO 409                                                          02630000
 407  RJI=JI                                                            02640000
      DJJ2=DJ(2)                                                        02650004
 401  IF(DJJ2.LE.0.) GOTO 402                                           02660004
      JDUM=10.*DJJ2+0.5                                                 02670005
      IF(JI.GT.JDUM) GOTO 403                                           02680004
      DJJ2=DJJ2-1.                                                      02690004
      GOTO 401                                                          02700004
 402  DJJ2=0.                                                           02710004
 403  JFMIN1=RJI-DJJ2  +1.5                                             02720004
      JFMIN1=MAX0(1,JFMIN1)                                             02730000
      JFMAX1=JFMIN1                                                     02740000
      LMAX2=LN+1                                                        02750000
      LMIN1=LMAX2                                                       02760000
 409  LJF=1                                                             02770000
      JFMINN=JFMIN1                                                     02780000
      MMWN=MWN                                                          02790000
      IF(NOLJF.EQ.1 .OR. IEF.LT.20) GOTO 410                            02800000
      IDUM=JFMAX1-JFMIN1                                                02810000
      IF(IDUM.LT.6) GOTO 410                                            02820000
      LJF=2                                                             02830000
      MMWN=MWN+MWN                                                      02840000
      IADD1=1-IADD1                                                     02850000
      IDUM=JI1-JFMIN1+IADD1                                             02860000
      JFMIN1=JFMIN1+MOD(IDUM,2)                                         02870000
      JFMINN=-JFMIN1                                                    02880000
 410  KNDEX=0                                                           02890000
      JFMIN = JFMIN1-1                                                  02900000
      JROWS=JFMIN*MWN                                                   02910000
C                                                                       02920000
C************************************** JD LOOP  N    ******************02930000
C                                                                       02940000
      DO 440 JF1=JFMIN1,JFMAX1,LJF                                      02950000
      IF(JFJI.EQ.1)GOTO 420                                             02960000
      LMAX2=JI1+JF1                                                     02970000
      LMIN1=IABS(JI1-JF1)+1                                             02980000
C---------------------------------------------                          02990000
 420  DO 430 IE=1,MAXN                                                  03000000
      IW=KONSTN+IE                                                      03010000
      IF(IW.GT.MWN) GO TO 435                                           03020000
      IW=IW+JROWS                                                       03030000
      IF(OMEGN(IW).LE.0.)   GO TO 435                                   03040000
      B=TCOFN(IE,LMIN1)                                                 03050000
      IF(LMAX2.LE.LN ) B=B-TCOFN(IE,LMAX2)                              03060000
      B=B*OMEGN(IW)                                                     03070000
      TSPECN(IE)      =TSPECN(IE)      +B                               03080000
      INDEX=KNDEX+IE                                                    03090000
      TMTBN(INDEX)=B                                                    03100000
  430 CONTINUE                                                          03110000
C---------------------------------------------                          03120000
  435 KNDEX=KNDEX+MAXN                                                  03130000
      JROWS=JROWS+MMWN                                                  03140000
  440 CONTINUE                                                          03150000
C                                                                       03160000
C     IF(NSTEP.GT.1 .OR. JI.GT.7) GOTO 459                              03170000
C     WRITE(6,450) JFMIN1,JFMAX1,LJF,JFJI,MAXN,KONSTN,MWN,NWN,M2,N2     03180000
C450  FORMAT(' EVA',10I6)                                               03190000
C     IT1=1                                                             03200000
C     DO 456 JF1=JFMIN1,JFMAX1,LJF                                      03210000
C     IT2=IT1+MAXN-1                                                    03220000
C     WRITE(6,455) (TMTBN(I8),I8=IT1,IT2)                               03230000
C455  FORMAT(10E12.4)                                                   03240000
C     IT1=IT1+MAXN                                                      03250000
C456  CONTINUE                                                          03260000
C459  CONTINUE                                                          03270000
C**************************************** END JD LOOP  N   *************03280000
C                                                                       03290000
C------------------------------------NEUTRON INTEGRATIONS------         03300000
      P(2)=0.                                                           03310000
      SUM=0.                                                            03320000
C     NEUTRON INTEGRATIONS HAVE NO IMPLIED ZEROS                        03330000
      DO 470 K= 2 ,INL,2                                                03340000
  470 P(2)=P(2)+TSPECN(K)                                               03350000
      HIST(2)=P(2)                                                      03360000
      P(2)=P(2)*2.                                                      03370000
      DO 480 K= 3 ,INL,2                                                03380000
  480 SUM =SUM +TSPECN(K)                                               03390000
      HIST(2)=HIST(2)+SUM+TSPECN(1)                                     03400000
      P(2)=P(2)+SUM                                                     03410000
      P(2)=(2.*P(2)+TSPECN(1))/3.                                       03420000
      P(2)=P(2)*LJF                                                     03430000
      HIST(2)=HIST(2)*LJF                                               03440000
      GOTO 500                                                          03450000
  490 P(2)=0.                                                           03460000
      HIST(2)=0.                                                        03470000
C-----------------------------------------------PROTONS---------        03480000
  500 IF(MAXP.EQ.0) GOTO 590                                            03490000
      IF(EXMAX3.LE.0.) GOTO 590                                         03500000
      DO 504 K=1,IPL                                                    03510000
  504 TSPECP(K)=0.                                                      03520000
      DO 508 K=1,ITP                                                    03530000
  508 TMTBP (K)=0.                                                      03540000
      IF(JFJI.EQ.1) GOTO507                                             03550000
      JFMAX1=MIN0(NWP,JI+LP )                                           03560000
      JFMIN1=MAX0(1,JI+2-LP )                                           03570000
      GOTO 509                                                          03580000
 507  RJI=JI                                                            03590000
      DJJ3=DJ(3)                                                        03600004
 501  IF(DJJ3.LE.0.) GOTO 502                                           03610004
      JDUM=10.*DJJ3+0.5                                                 03620005
      IF(JI.GT.JDUM) GOTO 503                                           03630004
      DJJ3=DJJ3-1.                                                      03640004
      GOTO 501                                                          03650004
 502  DJJ3=0.                                                           03660004
 503  JFMIN1=RJI-DJJ3  +1.5                                             03670004
      JFMIN1=MAX0(1,JFMIN1)                                             03680000
      JFMAX1=JFMIN1                                                     03690000
      LMAX2=LP+1                                                        03700000
      LMIN1=LMAX2                                                       03710000
 509  LJF=1                                                             03720000
      JFMINP=JFMIN1                                                     03730000
      MMWP=MWP                                                          03740000
      IF(NOLJF.EQ.1 .OR. IEF.LT.20) GOTO 510                            03750000
      IDUM=JFMAX1-JFMIN1                                                03760000
      IF(IDUM.LT.6) GOTO 510                                            03770000
      LJF=2                                                             03780000
      MMWP=MWP+MWP                                                      03790000
      IADD2=1-IADD2                                                     03800000
      IDUM=JI1-JFMIN1+IADD2                                             03810000
      JFMIN1=JFMIN1+MOD(IDUM,2)                                         03820000
      JFMINP=-JFMIN1                                                    03830000
 510  KNDEX=0                                                           03840000
      JFMIN = JFMIN1-1                                                  03850000
      JROWS=JFMIN*MWP                                                   03860000
C                                                                       03870000
C************************************** JD LOOP  P    ******************03880000
C                                                                       03890000
      DO 540 JF1=JFMIN1,JFMAX1,LJF                                      03900000
      IF(JFJI.EQ.1)GOTO 520                                             03910000
      LMAX2=JI1+JF1                                                     03920000
      LMIN1=IABS(JI1-JF1)+1                                             03930000
C-----------------------------------------                              03940000
 520  DO 530 IE=1,MAXP                                                  03950000
      IW=KONSTP+IE                                                      03960000
      IF(IW.GT.MWP) GO TO 535                                           03970000
      IW=IW+JROWS                                                       03980000
      IF(OMEGP(IW).LE.0.)   GO TO 535                                   03990000
      B=TCOFP(IE,LMIN1)                                                 04000000
      IF(LMAX2.LE.LP ) B=B-TCOFP(IE,LMAX2)                              04010000
      B=B*OMEGP(IW)                                                     04020000
      TSPECP(IE)      =TSPECP(IE)      +B                               04030000
      INDEX=KNDEX+IE                                                    04040000
      TMTBP(INDEX)=B                                                    04050000
  530 CONTINUE                                                          04060000
C-----------------------------------------                              04070000
  535 KNDEX=KNDEX+MAXP                                                  04080000
      JROWS=JROWS+MMWP                                                  04090000
  540 CONTINUE                                                          04100000
C                                                                       04110000
C**************************************** END JD LOOP  P   *************04120000
C                                                                       04130000
C------------------------------------PROTON INTEGRATIONS---------       04140000
      P(3)=0.                                                           04150000
      SUM=0.                                                            04160000
C     ZEROS IMPLIED AT BOTH ENDS OF CHARGED PARTICLE SPECTRA            04170000
      DO 570 K= 1 ,IPL,2                                                04180000
      P(3)=P(3)+TSPECP(K)                                               04190000
  570 CONTINUE                                                          04200000
      HIST(3)=P(3)                                                      04210000
      P(3)=2.*P(3)                                                      04220000
      DO 580 K= 2 ,IPL,2                                                04230000
      SUM =SUM +TSPECP(K)                                               04240000
  580 CONTINUE                                                          04250000
      HIST(3)=HIST(3)+SUM                                               04260000
      P(3)=(P(3)+SUM)*0.6666667                                         04270000
      P(3)=P(3)*LJF                                                     04280000
      HIST(3)=HIST(3)*LJF                                               04290000
C     END PROTON EVAPORATION                                            04300000
      GOTO 600                                                          04310000
  590 P(3)=0.                                                           04320000
      HIST(3)=0.                                                        04330000
C---------------------------------------------ALPHAS---------------     04340000
 600  IF(MAXA.EQ.0) GOTO 690                                            04350000
      IF(EXMAX4.LE.0.) GOTO 690                                         04360000
      DO 604 K=1,IAL                                                    04370000
 604  TSPECA(K)=0.                                                      04380000
      DO 608  K=1,ITA                                                   04390000
 608  TMTBA(K)=0.                                                       04400000
      IF(JFJI.EQ.1) GOTO 607                                            04410000
      JFMIN1=MAX0(1,JI+2-LA )                                           04420000
      JFMAX1=MIN0(JI+LA ,NWA)                                           04430000
      GOTO 609                                                          04440000
 607  RJI=JI                                                            04450000
      DJJ4=DJ(4)                                                        04460004
 601  IF(DJJ4.LE.0.) GOTO 602                                           04470004
      JDUM=5.*DJJ4+0.5                                                  04480004
      IF(JI.GT.JDUM) GOTO 603                                           04490004
      DJJ4=DJJ4-1.                                                      04500004
      GOTO 601                                                          04510004
 602  DJJ4=0.                                                           04520004
 603  JFMIN1=RJI-DJJ4  +1.5                                             04530004
      JFMIN1=MAX0(1,JFMIN1)                                             04540000
      JFMAX1=JFMIN1                                                     04550000
      LMAX2=LA+1                                                        04560000
      LMIN1=LMAX2                                                       04570000
 609  LJF=1                                                             04580000
      JFMINA=JFMIN1                                                     04590000
      MMWA=MWA                                                          04600000
      IF(NOLJF.EQ.1 .OR. IEF.LT.20) GOTO 610                            04610000
      IDUM=JFMAX1-JFMIN1                                                04620000
      IF(IDUM.LT.6) GOTO 610                                            04630000
      LJF=2                                                             04640000
      MMWA=MWA+MWA                                                      04650000
      IADD3=1-IADD3                                                     04660000
      IDUM=JI1-JFMIN1+IADD3                                             04670000
      JFMIN1=JFMIN1+MOD(IDUM,2)                                         04680000
      JFMINA=-JFMIN1                                                    04690000
 610  KNDEX=0                                                           04700000
      JFMIN = JFMIN1-1                                                  04710000
      JROWS=JFMIN*MWA                                                   04720000
C                                                                       04730000
C************************************* JD LOOP ALPHAS ******************04740000
C                                                                       04750000
      DO 640 JF1=JFMIN1,JFMAX1,LJF                                      04760000
      IF(JFJI.EQ.1)GOTO 620                                             04770000
      LMAX2=JI1+JF1                                                     04780000
      LMIN1=IABS(JI1-JF1)+1                                             04790000
C------------------------------------------------                       04800000
 620  DO 630 KE=1,MAXA                                                  04810000
      IW=KONSTA+KE                                                      04820000
      IF(IW.GT.MWA) GO TO 635                                           04830000
      IW=IW+JROWS                                                       04840000
      IF(OMEGA(IW).LE.0.)   GO TO 635                                   04850000
      B=TCOFA(KE,LMIN1)                                                 04860000
      IF(LMAX2.LE.LA )B=B-TCOFA(KE,LMAX2)                               04870000
      B=B*OMEGA(IW)                                                     04880000
      TSPECA(KE)=TSPECA(KE)+B                                           04890000
      INDEX=KNDEX+KE                                                    04900000
      TMTBA(INDEX)=B                                                    04910000
  630 CONTINUE                                                          04920000
C------------------------------------------------                       04930000
  635 KNDEX=KNDEX+MAXA                                                  04940000
      JROWS=JROWS+MMWA                                                  04950000
  640 CONTINUE                                                          04960000
C                                                                       04970000
C**************************************** END JD LOOP ALPHAS ***********04980000
C----------------------------------ALPHA INTEGRATIONS-------            04990000
      P(4)=0.                                                           05000000
      SUM=0.                                                            05010000
      DO 670 KE= 1 ,IAL,2                                               05020000
      P(4)=P(4)+TSPECA(KE)                                              05030000
  670 CONTINUE                                                          05040000
      HIST(4)=P(4)                                                      05050000
      P(4)=2.*P(4)                                                      05060000
      DO 680 KE=2  ,IAL,2                                               05070000
      SUM =SUM +TSPECA(KE)                                              05080000
  680 CONTINUE                                                          05090000
      HIST(4)=HIST(4)+SUM                                               05100000
      P(4)=(P(4)+SUM)*0.6666667                                         05110000
      P(4)=P(4)*LJF                                                     05120000
      HIST(4)=HIST(4)*LJF                                               05130000
      GOTO 700                                                          05140000
  690 P(4)=0.                                                           05150000
      HIST(4)=0.                                                        05160000
C------------------------------------------ GAMMAS -------------------- 05170000
C                                                                       05180000
  700 IF(NOG.EQ.1) GOTO 790                                             05190000
      IF(MAXG.LE.0) GO TO 790                                           05200000
      DO 704 KE=1,IGL                                                   05210000
      TSPECQ(KE)=0.                                                     05220000
  704 TSPECG(KE)=0.                                                     05230000
      DO 708 K=1,ITG                                                    05240000
  708 TMTBG(K)=0.                                                       05250000
      DUMQ=0.                                                           05260000
      LGD=2                                                             05270000
      KE2M=KE2MX                                                        05280000
C     IF(IEF.GT.20) KE2M=0                                              05290000
      IF(KE2M.EQ.0) LGD=1                                               05300000
      JFMIN=MAX0(JI-LGD,0)                                              05310000
      JFMIN1=JFMIN+1                                                    05320000
      JFMAX1=MIN0(JI+1+LGD,NWG)                                         05330000
      GAMJ=1.                                                           05340000
      IF(JFACTR.EQ.1)GAMJ=JI+JI+1                                       05350000
      IF(JFMAX1.LT.JFMIN1) GOTO 790                                     05360000
      JROWS=JFMIN*MWG                                                   05370000
      KNDEX=0                                                           05380000
C************************************** JD LOOP GAMMAS **************   05390000
C                                                                       05400000
      DO 740 JF1=JFMIN1,JFMAX1                                          05410000
      JDG=IABS(JI1-JF1)                                                 05420000
C----------------------------------------------                         05430000
      DO 730 KE=1,MAXG                                                  05440000
      IW=KONSTG+KE                                                      05450000
      IF(IW.GT.MWG)GOTO 735                                             05460000
      IW=IW+JROWS                                                       05470000
      IF(OMEGG(IW).LE.0.)GOTO 735                                       05480000
      DUMD=0.                                                           05490000
      IF(JDG.LE.1) DUMD=ECUBE(KE)*OMEGG(IW)*GAMJ                        05500000
      TSPECG(KE)=TSPECG(KE)+DUMD                                        05510000
      IF(KE.GT.KE2M) GOTO 720                                           05520000
      DUMQ=EFIVE(KE)*OMEGG(IW)                                          05530000
      TSPECQ(KE)=TSPECQ(KE)+DUMQ                                        05540000
      DUMD=DUMD+DUMQ                                                    05550000
 720  INDEX=KNDEX+KE                                                    05560000
      TMTBG(INDEX)=DUMD                                                 05570000
 730  CONTINUE                                                          05580000
C----------------------------------------------                         05590000
 735  KNDEX=KNDEX+MAXG                                                  05600000
      JROWS=JROWS+MWG                                                   05610000
 740  CONTINUE                                                          05620000
C*********************************** END JD LOOP GAMMAS ************    05630000
C                                                                       05640000
C--------------------------- GAMMA INTEGRATIONS-----------              05650000
      P(1)=0.                                                           05660000
      SUM=0.                                                            05670000
      DO 770 KE=1,IGL,2                                                 05680000
 770  P(1)=P(1)+TSPECQ(KE)+TSPECG(KE)                                   05690000
      HIST(1)=P(1)                                                      05700000
      P(1)=2.*P(1)                                                      05710000
      DO 780 KE=2,IGL,2                                                 05720000
 780  SUM =SUM +TSPECG(KE)+TSPECQ(KE)                                   05730000
      HIST(1)=HIST(1)+SUM                                               05740000
      P(1)=(P(1)+SUM)*0.6666667                                         05750000
      GO TO 800                                                         05760000
 790  P(1)=0.                                                           05770000
      HIST(1)=0.                                                        05780000
C----------------------------------------------- FISSION-------------   05790000
C                                                                       05800000
 800  IF(MAXF.EQ.0)GOTO 890                                             05810000
      IF(JI.GE.NWF) GOTO 890                                            05820000
      DO 804 KE=1,KEF                                                   05830000
 804  TSPECF(KE)=0.                                                     05840000
      JROWS=JI*MWF                                                      05850000
      B=0.                                                              05860000
      DO 830 KE=1,MAXF                                                  05870000
      IW=KONSTF+KE                                                      05880000
      IF(IW.GT.MWF) GOTO 835                                            05890000
      IW=IW+JROWS                                                       05900000
      IF(OMEGF(IW).LE.0.)GOTO 835                                       05910000
      TSPECF(KE)=OMEGF(IW)                                              05920000
 830  CONTINUE                                                          05930000
C------------------------------------                                   05940000
 835  SUM=0.                                                            05950000
      SUM1=0.                                                           05960000
      IF(MAXF.EQ.1) GOTO 885                                            05970000
      DO 870 KE=2,MAXF,2                                                05980000
 870  SUM=SUM+TSPECF(KE)                                                05990000
      IF(MAXF.EQ.2) GOTO 885                                            06000000
      DO 880 KE=3,MAXF,2                                                06010000
 880  SUM1=SUM1+TSPECF(KE)                                              06020000
 885  P(5)=0.333333*TSPECF(1)+1.333333*SUM+0.6666667*SUM1               06030000
C
C Behinderung der Spaltung, fuer E* > EDELFIS Spaltung behindert
C
      IF(INOF.NE.1) GOTO 889
      IF(EXMAX1.LT.EDELFIS) GOTO 889
      BEHIF1=0.6931
      EXPFAK=BEHIF1*(EXMAX1-EDELFIS)/ADELFIS
      BEHIFAK=1./EXP(EXPFAK)
      P(5)=P(5)*BEHIFAK
C   
 889  HIST(5)=SUM+SUM1+TSPECF(1)                                        06040000
      GOTO 900                                                          06050000
 890  P(5)=0.                                                           06060000
      HIST(5)=0.                                                        06070000
C********************************************************************** 06080000
C                                                                       06090000
 900  IF(M1.LE.0) WRITE (6,125) JPEP,IEP,NSTEP,JI                       06100000
      CALL POPUL(JPEP,IEP,NSTEP,JI)                                     06110000
      IF(M1.LE.0) WRITE (6,125) JPEP,IEP,NSTEP,JI                       06120000
 125  FORMAT(' POPUL ',10I10)                                           06130000
C                                                                       06140000
      IF(ILAND.NE.3 .OR.IOPT.NE.0) GOTO 1000                            06150000
      PALL=P(1)+P(2)+P(3)+P(4)+P(5)                                     06160000
      IF(PALL.LE.0.) GOTO 1000                                          06170000
      IW=IEP+JI*MWG                                                     06180000
      IDUM=IW+DEL(1)+0.5                                                06190000
      IW1=IW                                                            06200000
      IF(MWG.GE.IDUM)IW1=IDUM                                           06210000
      IF(OMEGG(IW1).LE.0.)GOTO 1000                                     06220000
      DUM=OMEGG(IW)/OMEGG(IW1)                                          06230000
      FLAND(JI1)=HCST*DUM*OMEGG(IW)/ PALL                               06240000
 1000 CONTINUE                                                          06250000
C                                                                       06260000
C********************************* END JP LOOP *************************06270000
C                                                                       06280000
 2000 CONTINUE                                                          06290000
C                                                                       06300000
C********************************* END EP LOOP *************************06310000
C                                                                       06320000
      IF(RLOSTF(NSTEP).GT.0.) AVGJF(NSTEP)=AVGJF(NSTEP)/RLOSTF(NSTEP)   06330000
      IF(RLOSTF(NSTEP).GT.0.) AVGEF(NSTEP)=AVGEF(NSTEP)/RLOSTF(NSTEP)   06340000
      IF(SSP(2)       .GT.0.) AVESP(2    )=AVESP(2    )/SSP(2)          06350000
      IF(SSP(3)       .GT.0.) AVESP(3    )=AVESP(3    )/SSP(3)          06360000
      IF(SSP(4)       .GT.0.) AVESP(4    )=AVESP(4    )/SSP(4)          06370000
      EAVN(NSTEP)=AVESP(2)                                              06380000
      EAVP(NSTEP)=AVESP(3)                                              06390000
      EAVA(NSTEP)=AVESP(4)                                              06400000
      MN=M1*N1                                                          06410000
      TRIM(NSTEP)=TRIM(NSTEP)+CUTS                                      06420000
      IF(LPRINT.LT.3)                                                   06430000
     1WRITE(6,5000) MN,ICU,IOFF,CUTPOP,CUTS,BELOW                       06440008
 5000 FORMAT(/' EVA  M1*N1',I6,4X,'ZEROS',I6,4X,'CUTS',I6,4X,'CUTPOP',  06450000
     1        E12.4,4X,'SUM CUTS',E12.4 /' LOSS BELOW YRAST LINE',E12.4)06460008
      NOG=NOGG                                                          06470000
      RETURN                                                            06480000
      END                                                               06490000
