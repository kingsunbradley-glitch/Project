C        KC02.GRO.FORT(POPUL)                                           00010075
      SUBROUTINE POPUL(JPEP,IEP,NSTEP,JI)                               00020053
C
C  Wird von 'EVA' aufgerufen ( innnerhalb EP,JP - Schleife )
C  Stand: 3.5.1994
C
       INCLUDE 'common.f'                                                00040076
C     MAYBE UNINITIALIZED: PER0
C-----------------------------------------------------------------------00050045
C
C     I.Giese 2007
      SAVE
C                                                                       00060005
      RJI=JI                                                            00070021
      JP1=JI+1                                                          00080068
      E=IEP-1                                                           00090058
      EXCURR=EXMAX1-E                                                   00100058
C                                                                       00110021
      IF(P(1).GE.-0.5) GOTO 40                                          00120013
      FLAND(JP1)=FLAND(JP1)+POP(JPEP)                                   00130008
      PALL=0.                                                           00140008
      GOTO 45                                                           00150008
 40   PALL=P(2)+P(3)+P(4)+P(5)                                          00160015
      PER0=P(1)+P(2)+P(3)+P(4)                                          00170075
      IF(P(1).GT.0.) GOTO 50                                            00180015
C                                                                       00190021
 45   IF(NOG.EQ.0) GOTO 47                                              00200021
      IF(PALL.GT.0.)GOTO 50                                             00210021
 47   TABLE(JPEP)=SYMB(6)                                               00220021
C                           SYMB(6) *   SYMB(7) X                       00230024
      IF((EXCURR+.1).LT.EJAY(JI+1)) GOTO 48                             00240065
C     IF(P(5).GT.0.) GOTO 48                                            00250062
      SIG(NSTEP)=SIG(NSTEP)+POP(JPEP)                                   00260002
      IF(JP1.LT.3) GOTO 46                                              00270069
      JP11=JP1+2                                                        00280069
      DO 43 JPP=3,JP1,2                                                 00290069
      JP11=JP11-2                                                       00300069
      IEG=EJAY(JP11)-EJAY(JP11-2)+1.                                    00310069
      IF(IEG.LT.11) SPECQ(IEG)=SPECQ(IEG)+POP(JPEP)                     00320069
 43   CONTINUE                                                          00330069
C 
C     I.Giese 2007
C     Abfrage auf ISO.NE.0
C
 46   IF (ISO.NE.0 .AND. JI.GT.ISOJ(ISO))
     1                    SIGISO(ISO)=SIGISO(ISO)+POP(JPEP)             00340069
C
C     I.Giese 2007: Index 0 taucht etwa 900 Mal auf !!!
C     COMMON/ISOM/NUMISO,ISO,ISOA(5),ISOZ(5),ISOJ(5),SIGISO(5)
C
C     Nur Ausgabe wegen moeglicher Fehlerquelle
C     IF(ISO.LT.1.OR.ISO.GT.5) WRITE(6,5400) ISO,JI,ISOJ(ISO),
C    1                                    SIGISO(ISO),ISOJ(4),ISOJ(5)
C5400 FORMAT('POPUL-ISO-INDEX: ',3(I12,1X),E12.2,1X,2(I12,1X))
C
      IF(PALL.EQ.0.)RETURN                                              00350003
 48   RP(NSTEP)=RP(NSTEP)+POP(JPEP)                                     00360058
      TABLE(JPEP)=SYMB(7)                                               00370008
      PENTRY(JPEP)=0.                                                   00380069
      IF(P(5).GT.0.) GOTO 67                                            00390060
      RETURN                                                            00400060
C---------------------------------------------------------------        00410024
 50   PMAX=0.                                                           00420057
      DO 60 K=1,5                                                       00430003
 60   PMAX=DMAX1(PMAX,P(K))                                             00440015
      DO 65 K=1,5                                                       00450003
      IF(PMAX.GT.P(K)) GOTO 65                                          00460015
      TABLE(JPEP)=SYMB(K)                                               00470008
 65   CONTINUE                                                          00480003
      IF(PALL.GT.0.) PENTRY(JPEP) = PENTRY(JPEP) * P(1)/(PALL+P(1))     00490069
 67   PALL=PALL+P(1)                                                    00500057
      FACT=POP(JPEP)/PALL                                               00510002
C--------------------------------------------------------------------   00520000
C                                         GAMMAS                        00530000
 100  IF(HIST(1).LE.0.) GOTO 200                                        00540000
      FACTOR=FACT*P(1)/HIST(1)                                          00550015
      JFMAX1=MIN0(NWG,JI+LGD+1)                                         00560009
      JF1UP=JFMAX1                                                      00570044
      JFMIN1=MAX0(1,JI+1-LGD)                                           00580008
      IF(JFMIN1.GT.JF1UP) GOTO 130                                      00590009
      MUP=M1-IEP                                                        00600012
      IF(MUP.EQ.0  .AND.  TABLE(JPEP).EQ.SYMB(1)) GOTO 47               00610057
      KEUP=MIN0(MAXG,MUP)                                               00620009
      IF(KEUP.LT.1) GOTO 130                                            00630035
      IPOP=(JFMIN1-1)*M1+IEP                                            00640009
      ITM=0                                                             00650009
C                                                                       00660000
      DO 120 JF1=JFMIN1,JF1UP                                           00670000
      DO 110 KE=1,KEUP                                                  00680000
      IPO =IPOP+KE                                                      00690009
      IT =ITM+KE                                                        00700009
      POP (IPO )=POP (IPO )+TMTBG(IT )*FACTOR                           00710009
      TMTBG(IT )=0.                                                     00720009
 110  CONTINUE                                                          00730000
      ITM=ITM+MAXG                                                      00740009
      IF(JF1.LT.N1) IPOP=IPOP+M1                                        00750044
 120  CONTINUE                                                          00760000
C                                                                       00770000
 130  SUM=0.                                                            00780009
      DO 140 ITM=1,ITG                                                  00790008
 140  SUM=SUM+TMTBG(ITM)                                                00800000
      RLOSTG(NSTEP)=RLOSTG(NSTEP)+SUM*FACTOR                            00810002
      IE=DEL(1)                                                         00820047
C                                                                       00830000
      DO 150 KE=1,MAXG                                                  00840000
      IE=IE+1                                                           00850047
      SPECG(IE)=SPECG(IE)+ TSPECG(KE)*FACT                              00860003
      SPECQ(IE)=SPECQ(IE)+ TSPECQ(KE)*FACT                              00870047
 150  CONTINUE                                                          00880047
C--------------------------------------------------------------------   00890000
C                                         NEUTRONS                      00900000
 200  IF(HIST(2).LE.0.) GOTO 300                                        00910000
      FACTOR=FACT*P(2)/HIST(2)                                          00920015
      LJF=1                                                             00930028
      ITWO=0                                                            00940034
      IF(JFJI.EQ.1) GOTO 204                                            00950021
      JFMAX1=MIN0(NWN,JI+LN )                                           00960004
      JF1UP=MIN0(N2,JFMAX1)                                             00970009
      JFMIN1=MAX0(1,JI+2-LN)                                            00980034
      IF(JFMINN.LT.0)LJF=2                                              00990022
      IF(LJF.EQ.2) ITWO=1                                               01000034
      IDUM=IABS(JFMINN)                                                 01010051
      IF(IDUM.GT.JFMIN1) ITWO=-1                                        01020051
      JFMIN1=IDUM                                                       01030051
      IF(JFMIN1.GT.JF1UP) GOTO 230                                      01040011
      GOTO 206                                                          01050021
 204  IDUM=DJJ2                                                         01060070
      FDUM=IDUM                                                         01070042
      FDUM=DJJ2-FDUM                                                    01080070
C     IF(FDUM.GT.0.25 .AND. FDUM.LT.0.75) ITWO=1                        01090077
      JF1=RJI-DJJ2+1.5                                                  01100070
      JFMIN1=MAX0(1,JF1)                                                01110021
      JF1UP=JFMIN1                                                      01120021
 206  MUP=M2-IEP+1                                                      01130021
      KEUP=MIN0(MAXN,MUP)                                               01140009
      IF(KEUP.LT.1) GOTO 230                                            01150035
      IPOP=(JFMIN1-1)*M2 +IEP-1                                         01160009
      ITM=0                                                             01170009
      IDUM=JI+1-LN                                                      01180019
      MM2=M2                                                            01190022
      IF(LJF.EQ.2)MM2=M2+M2                                             01200022
C                                                                       01210030
      DO 220 JF1=JFMIN1,JF1UP,LJF                                       01220051
      JID=JF1-IDUM                                                      01230019
      IF(JF1.EQ.N2)ITWO=0                                               01240034
      DO 210 KE=1,KEUP                                                  01250027
      IPO =IPOP+KE                                                      01260009
      IT =ITM+KE                                                        01270009
      DUM=TMTBN(IT)*FACTOR                                              01280017
      POPN(IPO )= POPN(IPO )+DUM                                        01290017
      DJN(JID)=DJN(JID)+DUM                                             01300019
      TMTBN(IT )=0.                                                     01310009
      IF(ITWO) 205,210,207                                              01320051
 207  IPO1=IPO+M2                                                       01330051
      JID1=JID+1                                                        01340052
      GOTO 208                                                          01350051
 205  IPO1=IPO-M2                                                       01360051
      JID1=MAX0(1,JID-1)                                                01370052
 208  POPN(IPO1)=POPN(IPO1)+DUM                                         01380051
      DJN(JID1)=DJN(JID1)+DUM                                           01390051
 210  CONTINUE                                                          01400000
      ITM=ITM+MAXN                                                      01410009
      IPOP=IPOP+MM2                                                     01420022
 220  CONTINUE                                                          01430000
C                                                                       01440000
 230  SUM=0.                                                            01450011
      IDUM=ITN/LJF                                                      01460022
      DO 240 ITM=1,IDUM                                                 01470022
 240  SUM=SUM+TMTBN(ITM)                                                01480000
      RLOSTN(NSTEP)=RLOSTN(NSTEP)+SUM*FACTOR*FLOAT(LJF)                 01490034
C                                                                       01500000
      FACT1=FACT*FLOAT(LJF)                                             01510027
      DUM=0.                                                            01520021
      DO 250 KE=1,MAXN                                                  01530027
      FIE= DEL(2)+KE -1                                                 01540038
      IE=FIE+1                                                          01550038
      DUM=DUM+TSPECN(KE)*FIE                                            01560021
 250  SPECN(IE)=SPECN(IE)+TSPECN(KE)*FACT1                              01570024
      AVESP(2)=AVESP(2)+DUM*FACT1                                       01580024
      SSP(2)=SSP(2)+HIST(2)*FACT                                        01590021
C--------------------------------------------------------------------   01600000
C                                         PROTONS                       01610000
 300  IF(HIST(3).LE.0.) GOTO 400                                        01620000
      FACTOR=FACT*P(3)/HIST(3)                                          01630015
      LJF=1                                                             01640028
      ITWO=0                                                            01650034
      IF(JFJI.EQ.1) GOTO 304                                            01660021
      JFMAX1=MIN0(NWP,JI+LP )                                           01670004
      JF1UP=MIN0(N3,JFMAX1)                                             01680009
      JFMIN1=MAX0(1,JI+2-LP)                                            01690034
      IF(JFMINP.LT.0)LJF=2                                              01700023
      IF(LJF.EQ.2)ITWO=1                                                01710034
      IDUM=IABS(JFMINP)                                                 01720051
      IF(IDUM.GT.JFMIN1) ITWO=-1                                        01730051
      JFMIN1=IDUM                                                       01740051
      IF(JFMIN1.GT.JF1UP) GOTO 330                                      01750011
      GOTO 306                                                          01760021
 304  IDUM=DJJ3                                                         01770070
      FDUM=IDUM                                                         01780042
      FDUM=DJJ3-FDUM                                                    01790070
C     IF(FDUM.GT.0.25 .AND. FDUM.LT.0.75) ITWO=1                        01800077
      JF1=RJI-DJJ3+1.5                                                  01810070
      JFMIN1=MAX0(1,JF1)                                                01820042
      JF1UP=JFMIN1                                                      01830042
 306  MUP=M3-IEP+1                                                      01840021
      KEUP=MIN0(MAXP,MUP)                                               01850009
      IF(KEUP.LT.1) GOTO 330                                            01860035
      IPOP=(JFMIN1-1)*M3+IEP-1                                          01870009
      ITM=0                                                             01880009
      IDUM=JI+1-LP                                                      01890020
      MM3=M3                                                            01900023
      IF(LJF.EQ.2)MM3=M3+M3                                             01910023
C                                                                       01920000
      DO 320 JF1=JFMIN1,JF1UP,LJF                                       01930023
      JID=JF1-IDUM                                                      01940020
      IF(JF1.EQ.N3) ITWO=0                                              01950034
      DO 310 KE=1,KEUP                                                  01960027
      IPO =IPOP+KE                                                      01970009
      IT =ITM+KE                                                        01980009
      DUM=TMTBP(IT)*FACTOR                                              01990017
      POPP(IPO )= POPP(IPO )+DUM                                        02000017
      DJP(JID)=DJP(JID)+DUM                                             02010020
      TMTBP(IT )=0.                                                     02020009
      IF(ITWO) 305,310,307                                              02030051
 307  IPO1=IPO+M3                                                       02040051
      JID1=JID+1                                                        02050052
      GOTO 308                                                          02060051
 305  IPO1=IPO-M3                                                       02070051
      JID1=MAX0(1,JID-1)                                                02080052
 308  POPP(IPO1)=POPP(IPO1)+DUM                                         02090051
      DJP(JID1)=DJP(JID1)+DUM                                           02100051
 310  CONTINUE                                                          02110000
      ITM=ITM+MAXP                                                      02120009
      IPOP=IPOP+MM3                                                     02130023
 320  CONTINUE                                                          02140000
C                                                                       02150000
 330  SUM=0.                                                            02160011
      IDUM=ITP/LJF                                                      02170023
      DO 340 ITM=1,IDUM                                                 02180023
 340  SUM=SUM+TMTBP(ITM)                                                02190000
      RLOSTP(NSTEP)=RLOSTP(NSTEP)+SUM*FACTOR*FLOAT(LJF)                 02200034
C                                                                       02210000
      FACT1=FACT*FLOAT(LJF)                                             02220027
      DUM=0.                                                            02230021
      DO 350 KE=1,MAXP                                                  02240027
      FIE= DEL(3)+KE-1                                                  02250038
      IE=FIE+1                                                          02260038
      DUM=DUM+TSPECP(KE)*FIE                                            02270021
 350  SPECP(IE)=SPECP(IE)+TSPECP(KE)*FACT1                              02280024
      AVESP(3)=AVESP(3)+DUM*FACT1                                       02290024
      SSP(3)=SSP(3)+HIST(3)*FACT                                        02300021
C--------------------------------------------------------------------   02310000
C                                         ALPHAS                        02320000
 400  IF(HIST(4).LE.0.) GOTO 500                                        02330000
      FACTOR=FACT*P(4)/HIST(4)                                          02340015
      LJF=1                                                             02350028
      ITWO=0                                                            02360034
      IF(JFJI.EQ.1) GOTO 404                                            02370021
      JFMAX1=MIN0(NWA,JI+LA )                                           02380004
      JF1UP=MIN0(N4,JFMAX1)                                             02390009
      JFMIN1=MAX0(1,JI+2-LA)                                            02400034
      IF(JFMINA.LT.0)LJF=2                                              02410023
      IF(LJF.EQ.2)ITWO=1                                                02420034
      IDUM=IABS(JFMINA)                                                 02430051
      IF(IDUM.GT.JFMIN1) ITWO=-1                                        02440051
      JFMIN1=IDUM                                                       02450051
      IF(JFMIN1.GT.JF1UP) GOTO 430                                      02460011
      GOTO 406                                                          02470021
 404  IDUM=DJJ4                                                         02480070
      FDUM=IDUM                                                         02490042
      FDUM=DJJ4-FDUM                                                    02500070
C     IF(FDUM.GT.0.25 .AND. FDUM.LT.0.75) ITWO=1                        02510077
      JF1=RJI-DJJ4+1.5                                                  02520070
      JFMIN1=MAX0(1,JF1)                                                02530042
      JF1UP=JFMIN1                                                      02540042
 406  MUP=M4-IEP+1                                                      02550021
      KEUP=MIN0(MAXA,MUP)                                               02560009
      IF(KEUP.LT.1) GOTO 430                                            02570035
      IPOP=(JFMIN1-1)*M4+IEP-1                                          02580009
      ITM=0                                                             02590009
      IDUM=JI+1-LA                                                      02600020
      MM4=M4                                                            02610023
      IF(LJF.EQ.2)MM4=M4+M4                                             02620023
C                                                                       02630000
      DO 420 JF1=JFMIN1,JF1UP,LJF                                       02640023
      JID=JF1-IDUM                                                      02650020
      IF(JF1.EQ.N4) ITWO=0                                              02660034
      DO 410 KE=1,KEUP                                                  02670027
      IPO =IPOP+KE                                                      02680009
      IT =ITM+KE                                                        02690009
      DUM=TMTBA(IT)*FACTOR                                              02700017
      POPA(IPO )= POPA(IPO )+DUM                                        02710017
      DJA(JID)=DJA(JID)+DUM                                             02720039
      TMTBA(IT )=0.                                                     02730009
      IF(ITWO) 405,410,407                                              02740051
 407  IPO1=IPO+M4                                                       02750051
      JID1=JID+1                                                        02760052
      GOTO 408                                                          02770051
 405  IPO1=IPO-M4                                                       02780051
      JID1=MAX0(1,JID-1)                                                02790052
 408  POPA(IPO1)=POPA(IPO1)+DUM                                         02800051
      DJA(JID1)=DJA(JID1)+DUM                                           02810051
 410  CONTINUE                                                          02820000
      ITM=ITM+MAXA                                                      02830009
      IPOP=IPOP+MM4                                                     02840023
 420  CONTINUE                                                          02850000
C                                                                       02860000
 430  SUM=0.                                                            02870011
      IDUM=ITA/LJF                                                      02880023
      DO 440 ITM=1,IDUM                                                 02890023
 440  SUM=SUM+TMTBA(ITM)                                                02900000
      RLOSTA(NSTEP)=RLOSTA(NSTEP)+SUM*FACTOR *FLOAT(LJF)                02910034
C                                                                       02920000
      FACT1=FACT*FLOAT(LJF)                                             02930027
      DUM=0.                                                            02940021
      DO 450 KE=1,MAXA                                                  02950027
      FIE= DEL(4)+KE -1                                                 02960038
      IE=FIE+1                                                          02970038
      DUM=DUM+TSPECA(KE)*FIE                                            02980021
 450  SPECA(IE)=SPECA(IE)+TSPECA(KE)*FACT1                              02990024
      AVESP(4)=AVESP(4)+DUM*FACT1                                       03000024
      SSP(4)=SSP(4)+HIST(4)*FACT                                        03010021
C---------------------------------------------------------------------- 03020001
C                                      FISSION                          03030001
 500  IF(      P(5).LE.0.) GOTO 600                                     03040015
      DUM=POP   (JPEP)*   P(5)/PALL                                     03050015
      RLOSTF(NSTEP)=RLOSTF(NSTEP)+DUM                                   03060010
      JI1=JI+1                                                          03070036
      YFISJ(JI1)=YFISJ(JI1)+DUM                                         03080036
      YPOPJ(JI1)=YPOPJ(JI1)+POP(JPEP)                                   03090071
      YERJ(JI1)=YERJ(JI1)+     POP(JPEP) * PER0/PALL                    03100075
      AVGJF(NSTEP)=AVGJF(NSTEP)+DUM*FLOAT(JI)                           03110017
      AVGEF(NSTEP)=AVGEF(NSTEP)+DUM*(EXMAX1+1-IEP)                      03120017
C                                                                       03130001
      DO 550 KE=1,MAXF                                                  03140001
      IE=KE+DEL(5)                                                      03150003
 550  SPECF(IE)=SPECF(IE)+TSPECF(KE)*FACT                               03160003
C---------------------------------------------------------------------- 03170001
 600  CONTINUE                                                          03180001
      RETURN                                                            03190001
      END                                                               03200001
