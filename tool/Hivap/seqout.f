      SUBROUTINE SEQOUT(NROW,NSTEP,NOPUT)                               00020008
C
C  Druckt Ausgabe innerhalb  einer Sequenz(N) - Schleife
C  letzte Modifikation 20.4.1990
C  Stand: 3.5.1994
C
      INCLUDE 'common.f'                                                00030073
C--------------------------------------------------------------------   00040076
C
C     I.Giese 2007
      SAVE
C
C                                                                       00050008
      LOGUN=6                                                           00060040
      MZ=KZ(NSTEP)                                                      00070012
      MAA=MZ+IN(NSTEP)                                                  00080012
      IF(LPRINT.LE.1)                                                   00090023
     1WRITE(LOGUN,5 ) NROW,NSTEP,MAA,ELEMNT(MZ),ELAB,EXCIT              00100040
 5    FORMAT(' SEQOUT NROW,NSTEP',2I4,9X,'PARENT ',I3,A4,4X,'ELAB',F8.1,00110015
     1 4X,'EXCIT',F8.1/)                                                00120012
      IF(NOPUT.EQ.1) GOTO 6                                             00130008
C                                                                       00140076
      IF(LPRINT.LT.3) THEN                                              00150076
         DO 31 K=1,KEG                                                  00160076
         SPECGT(K)=SPECGT(K)+SPECG(K)                                   00170076
 31      SPECQT(K)=SPECQT(K)+SPECQ(K)                                   00180076
         DO 32 K=1,KEN                                                  00190076
 32      SPECNT(K)=SPECNT(K)+ SPECN(K)                                  00200076
         DO 33 K=1,KEP1                                                 00210076
 33      SPECPT(K)=SPECPT(K)+SPECP(K)                                   00220076
         DO 34 K=1,KEA1                                                 00230076
 34      SPECAT(K)=SPECAT(K)+SPECA(K)                                   00240076
         DO 35 K=1,KEF                                                  00250076
 35      SPECFT(K)=SPECFT(K)+SPECF(K)                                   00260076
      ENDIF                                                             00270076
C                                                                       00280076
C                                     PUT P                             00290076
C                                                                       00300076
      MN=M3*N3                                                          00310008
      SUMP=0.                                                           00320008
      DO 190 I=1,MN                                                     00330008
 190  SUMP=SUMP+POPP(I)                                                 00340008
      IF(LPRINT.LE.1)                                                   00350023
     1WRITE(LOGUN,195) SUMP,SUMLOW                                      00360040
 195  FORMAT(' SUMP,SUMLOW',2E14.4)                                     00370008
      IF(  SUMP     .LT.SUMLOW) GOTO 200                                00380008
       STORE(NROW+1,NSTEP)  = STORE(NROW+1,NSTEP)  +1.0001              00390009
      IUNIT=11-IU1                                                      00400008
      IREAD=0                                                           00410008
      I   =MZ-1                                                         00420019
      II   =MAA-1                                                       00430019
      CALL GETPUT(POPP,IPOPS,IUNIT,IREAD,M3,N3,EXMAX3,I,II,ACN,ZCN)     00440019
      IF(LPRINT.LE.1)                                                   00450023
     1WRITE(LOGUN,196) IUNIT,II ,ELEMNT(I)                              00460040
 196  FORMAT(' WROTE POPP ON UNIT',I4,15X,I3,A4)                        00470015
C                                                                       00480008
C                                     PUT A                             00490008
C                                                                       00500008
 200  MN=M4*N4                                                          00510008
      SUMA=0.                                                           00520008
      DO 204 I=1,MN                                                     00530069
 204  SUMA=SUMA+POPA(I)                                                 00540069
      IF(LPRINT.LE.1)                                                   00550023
     1WRITE(LOGUN,208) SUMA,SUMLOW                                      00560069
 208  FORMAT(' SUMA,SUMLOW',2E14.4)                                     00570069
      IF(SUMA       .LT.SUMLOW) GOTO   6                                00580008
       STORE(NROW+2,NSTEP+2)= STORE(NROW+2,NSTEP+2)+2.0001              00590009
      IUNIT=13-IU2                                                      00600020
      IREAD=0                                                           00610008
      IDUM=MAA-4                                                        00620019
      I=MZ-2                                                            00630019
      CALL GETPUT(POPA,IPOPS,IUNIT,IREAD,M4,N4,EXMAX4,I,IDUM,ACN,ZCN)   00640019
      IF(LPRINT.LE.1)                                                   00650023
     1WRITE(LOGUN,212) IUNIT,IDUM,ELEMNT(I)                             00660069
 212  FORMAT(' WROTE POPA ON UNIT',I4,15X,I3,A4)                        00670069
C---------------------------------------------------------------------- 00680004
C                                                                       00690004
 6    DO 216 J=1,N1                                                     00700069
      YFISS(J)=YFISS(J)+YFISJ(J)                                        00710069
      DUM=1.                                                            00720071
      IF(YPOPJ(J).GT.0.1E-20) DUM=YERJ(J)/YPOPJ(J)                      00730077
      IF(PER(J).GT..1E-20) PER(J)=PER(J)*DUM                            00740078
      DUMF=0.                                                           00750069
      IF(YPOPJ(J).GT.0.1E-20) DUMF=YFISJ(J)/YPOPJ(J)                    00760077
      PFIS(J)=PFIS(J)+(1.-PFIS(J))*DUMF                                 00770069
 216  CONTINUE                                                          00780069
      IF(JFIS.EQ.0) GOTO 11                                             00790069
      IF(NROW.NE.1 .OR. NSTEP.NE.1) GOTO 226                            00800064
      IACN=ACN+0.01                                                     00810064
      IZCN=ZCN+0.01                                                     00820064
      IF(JFIS.LT.4) WRITE(19,223) IACN,ELEMNT(IZCN),TOTAL(1),ELAB,N1    00830072
 223  FORMAT('C CNJ',I7,A4,E12.4,'(MB)',4X,                             00840064
     1       'ELAB',F8.2,I6/'H: X-1   Y')                               00850068
      IF(JFIS.LT.4) WRITE(19,231)(J,SPECJ(J),J=1,N1)                    00860072
 226  IF(JFIS.EQ.1 .OR. JFIS.EQ.3) GOTO 228                             00870067
      GOTO 11                                                           00880067
 228  AMAA=MAA                                                          00890067
      ZMZ=MZ                                                            00900044
      AN=AMAA-ZMZ                                                       00910044
C
C     I.Giese 2007
C     0.0D0
C
      CALL FISROT(AMAA,ZMZ,AN,1.0D0,DELR,DELSP,ERO,0.0D0 )              00920043
      CALL FISROT(AMAA,ZMZ,AN,0.0D0,DELR,DELS0,ERO,0.0D0 )              00930043
      DUMFF=DELSP-DELS0                                                 00940043
      WRITE(19,230) MAA,ELEMNT(MZ),RLOSTF(NSTEP),N1,ELAB,DUMFF          00950041
 230  FORMAT('C FISSJ',I4,A4,E12.4,'(MB)',3X,'JMAX=',I4,3X,             00960061
     1       'ELAB',F8.2,4X,E12.4/'H: X-1   Y')                         00970064
      WRITE(19,231)(J,YFISJ(J),J=1,N1)                                  00980062
 231  FORMAT(4(I4,E12.4,';'))                                           00990061
C--------------------------------------                                 01000050
C                             CALCULATE APROXIMATE TOTAL GAMMA ENERGY   01010050
 11   SUM=0.                                                            01020050
      SUME=0.                                                           01030050
      SUMJ=0.                                                           01040051
      SUMEG(NSTEP)=0.                                                   01050050
      SUMJG(NSTEP)=0.                                                   01060051
      IF(NOG.EQ.1) GOTO 250                                             01070050
      DO 240 IEP=1,M1                                                   01080050
      E=IEP-1                                                           01090050
      EXCURR=EXMAX1-E                                                   01100050
      JPEPMX=(N1-1)*M1+IEP                                              01110050
      RL=0.                                                             01120053
      DO 236 JPEP=IEP,JPEPMX,M1                                         01130050
C     IF(PENTRY(JPEP).GT.0.) GOTO 234                                   01140065
      IF(TABLE(JPEP).EQ.SYMB(1) .OR. TABLE(JPEP).EQ.SYMB(6)) GOTO 234   01150065
      GOTO 235                                                          01160051
 234  SUM=SUM+PENTRY(JPEP)                                              01170050
      SUME=SUME+PENTRY(JPEP)*EXCURR                                     01180050
      SUMJ=SUMJ+PENTRY(JPEP)*RL                                         01190053
 235  RL=RL+1.                                                          01200054
 236  CONTINUE                                                          01210050
 240  CONTINUE                                                          01220050
      IF(SUM.GT.0.) SUMEG(NSTEP)=SUME/SUM                               01230050
      IF(SUM.GT.0.) SUMJG(NSTEP)=SUMJ/SUM                               01240051
C-----------------------------------------------                        01250050
 250  IF(LPRINT.GT.1) GOTO 22                                           01260061
      WRITE(LOGUN,12) TOTAL(NSTEP)                                      01270040
 12   FORMAT(' TOTAL',E14.4)                                            01280003
      WRITE(LOGUN,14) FLOST(NSTEP),TRIM(NSTEP),RLOSTG(NSTEP),           01290040
     1 RLOSTN(NSTEP),                                                   01300040
     2 RLOSTP(NSTEP),RLOSTA(NSTEP),RP(NSTEP)                            01310040
 14   FORMAT(       ' FLOST(SUMPOP),TRIM(CUTOFF),RLOST(POPUL),RP '/     01320016
     1 7E13.4)                                                          01330003
C----------------------------------------------------------------       01340060
      WRITE(LOGUN,20) MAA,ELEMNT(MZ)                                    01350040
 20   FORMAT(/' J SPECTRUM OF PARENT',I4,A4)                            01360017
      ZERO=0.                                                           01370031
C
C     I.Giese 2007
C     1.0D0
C
      CALL OUT1(SPECJ,N1,ZERO,1.0D0,6)                                  01380031
      WRITE(LOGUN,18) MAA,ELEMNT(MZ)                                    01390040
 18   FORMAT(/' J SPECTRUM OF FISS  ',I4,A4)                            01400036
      ZERO=0.                                                           01410036
      CALL OUT1(YFISJ,N1,ZERO,1.0D0,6)                                  01420036
 19   ZERO=EXMAX1+1-MIN0(IEDIM,M1)                                      01430046
      WRITE(LOGUN,21) MAA,ELEMNT(MZ)                                    01440040
 21   FORMAT(/' E SPECTRUM OF PARENT',I4,A4)                            01450023
      IDUM=MIN0(IEDIM,M1)                                               01460046
      CALL OUT1(SPECE,IDUM,ZERO,1.0D0,6)                                01470046
C----------------------------------------------  ENTRYJ,GSPEC  JFIS=2,3 01480061
 22   IF(JFIS.EQ.2 .OR. JFIS.EQ.3) GOTO 305                             01490075
      GOTO 330                                                          01500067
 305  CALL TOT(N1,M1,PENTRY,SPECJ,AMX,AVG,SIGGAM,SPECE,AVE,EXMAX1)      01510067
      WRITE(19,310) MAA,ELEMNT(MZ),AVG,AVE,SIGGAM,SIG(NSTEP)            01520061
 310  FORMAT('C ENTRYJ',I4,A4,1X,'AVJ=',F6.1,3X,'AVE=',F6.1,2X,         01530063
     1 'SIGENTRY',E10.3,2X,'SIGER',E10.3/'H: X-1   Y')                  01540063
      WRITE(19,231)(J,SPECJ(J),J=1,N1)                                  01550062
      DO 320 IE=1,KEG                                                   01560061
 320  SPECJ(IE)=SPECG(IE)+SPECQ(IE)                                     01570061
      WRITE(19,325)                                                     01580061
 325  FORMAT('C GSPEC  1,2 ETC MEV'/'H: X   Y,H0')                      01590061
      WRITE(19,231)(IE,SPECJ(IE),IE=1,16)                               01600061
C----------------------------------------------  DJ-SPECS               01610061
 330  ILOG=6                                                            01620061
      IF(LPRINT.GT.1) ILOG=0                                            01630031
      ZERO=1-LN                                                         01640027
      IF(LPRINT.LT.2)WRITE(LOGUN,23)                                    01650061
 23   FORMAT(/' DJN SPEC')                                              01660027
      IDUM=2.*LN+1                                                      01670027
      CALL OUT1(DJN,IDUM,ZERO,1.0D0,ILOG)                               01680031
      AVDJN(NSTEP)=ZERO                                                 01690031
      ZERO=1-LP                                                         01700027
      IF(LPRINT.LT.2)WRITE(LOGUN,25)                                    01710061
 25   FORMAT(/' DJP SPEC')                                              01720027
      IDUM=2.*LP+1                                                      01730027
      CALL OUT1(DJP,IDUM,ZERO,1.0D0,ILOG)                               01740031
      AVDJP(NSTEP)=ZERO                                                 01750031
      ZERO=1-LA                                                         01760027
      IF(LPRINT.LT.2)WRITE(LOGUN,27)                                    01770061
 27   FORMAT(/' DJA SPEC')                                              01780027
      IDUM=2.*LA+1                                                      01790027
      CALL OUT1(DJA,IDUM,ZERO,1.0D0,ILOG)                               01800031
      AVDJA(NSTEP)=ZERO                                                 01810031
C--------------------------------------------                           01820061
      IF(LPRINT.GT.2)GOTO 65                                            01830032
      IF(IOPT.EQ.0) GOTO 36                                             01840076
      WRITE(LOGUN,30) STRIPE                                            01850040
 30   FORMAT(/' J SPECTRUM(LINEAR) OF G RESIDUAL IN YRAST STRIPE',      01860007
     1        F6.1,' MEV'/)                                             01870007
      ZERO=0.                                                           01880031
      CALL OUT1(FLAND,N1,ZERO,1.0D0,6)                                  01890031
C
C     I.Giese 2007
C     0.0D0
C
      CALL OUT2(FLAND,N1,1,0,0,0,0.0D0,10,0,1)                          01900007
 36   CONTINUE                                                          01910076
C-------------------------------------------- INFO VARIOUS BRANCHES     01920060
      WRITE(LOGUN,60)                                                   01930040
 60   FORMAT(/3X,'N   M   EXMAX  AV.EX. NW  MW  EXMAXW',6X,'Q',6X,      01940023
     1 'AL   DELTA',5X,'AVJ',6X,'DJ',8X,'TOTAL  BRANCH')                01950010
      DELTAJ=0.                                                         01960024
      IRANCH=10000                                                      01970001
      WRITE(LOGUN,100)N1,M1,EXMAX1,AVGE(NSTEP),NWG,MWG,EXMAXG,          01980040
     1 Q(1),AL(1),                                                      01990040
     2   DELTA(1),AVGJ(NSTEP),DELTAJ,TOTAL(NSTEP),IRANCH,BLANK          02000040
C-N-----------------                                                    02010060
 65   CALL TOT(N2,M2,POPN,SPECJ,AMX,AVG,SIGNN(NSTEP),SPECE,AVE,EXMAX2)  02020032
      IRANCH=SIGNN(NSTEP)*10000./TOTAL(NSTEP) +0.5                      02030045
      DELTAJ=AVGJ(NSTEP)-AVG                                            02040024
      IF(LPRINT.LT.3)                                                   02050032
     1WRITE(LOGUN,100)N2,M2,EXMAX2,AVE,NWN,MWN,EXMAXN,Q(2),AL(2),       02060040
     2DELTA(2),                                                         02070040
     3            AVG        ,DELTAJ,SIGNN(NSTEP),IRANCH,SYMB(2)        02080040
C-P-----------------                                                    02090060
      CALL TOT(N3,M3,POPP,SPECJ,AMX,AVG,SIGP (NSTEP),SPECE,AVE,EXMAX3)  02100023
      IRANCH=SIGP(NSTEP)*10000./TOTAL(NSTEP) +0.5                       02110045
      DELTAJ=AVGJ(NSTEP)-AVG                                            02120024
      IF(LPRINT.LT.3)                                                   02130032
     1WRITE(LOGUN,100)N3,M3,EXMAX3,AVE,NWP,MWP,EXMAXP,Q(3),AL(3),       02140040
     2            DELTA(3),                                             02150040
     3            AVG        ,DELTAJ,SIGP (NSTEP),IRANCH,SYMB(3)        02160040
C-A-----------------                                                    02170060
      CALL TOT(N4,M4,POPA,SPECJ,AMX,AVG,SIGA (NSTEP),SPECE,AVE,EXMAX4)  02180023
      IF(LPRINT.GT.2) RETURN                                            02190032
      IRANCH=SIGA(NSTEP)*10000./TOTAL(NSTEP) +0.5                       02200045
      DELTAJ=AVGJ(NSTEP)-AVG                                            02210024
      WRITE(LOGUN,100)N4,M4,EXMAX4,AVE,NWA,MWA,EXMAXA,Q(4),AL(4),       02220040
     2            DELTA(4),                                             02230040
     1            AVG        ,DELTAJ,SIGA (NSTEP),IRANCH,SYMB(4)        02240032
C-F-----------------                                                    02250060
      IRANCH=SIGF(NSTEP)*10000./TOTAL(NSTEP) +0.5                       02260045
      QF=-BFS(NSTEP,NROW)                                               02270028
      WRITE(LOGUN,101) AVGEF(NSTEP),NWF,MWF,EXMAXF,QF ,AL(5),DELTA(5),  02280040
     1           AVGJF(NSTEP),RLOSTF  (NSTEP),IRANCH,SYMB(5)            02290023
C-G-----------------                                                    02300060
      IRANCH=SIG(NSTEP)*10000./TOTAL(NSTEP) +0.5                        02310045
      WRITE(LOGUN,102) SIG(NSTEP),IRANCH,SYMB(1)                        02320040
 102  FORMAT(80X,E14.3,I6,2X,A2)                                        02330029
 100  FORMAT(2I4,2F8.1,2I4,F8.1,5F8.2,E14.3,I6,2X,A2)                   02340023
 101  FORMAT(16X, F8.1,2I4,F8.1,4F8.2,8X,E14.3,I6,2X,A2)                02350025
      WRITE(LOGUN,105) (AVESP(K),K=2,4),(SSP(K),K=2,4)                  02360040
 105  FORMAT(/' SPECTRA N,P,A  AVE',3F8.2,8X,' SUMS',3E12.4)            02370028
C-------------------------------------------------------------          02380060
C PRINT-PLOTS (LPRINT=0)                                                02390060
      IF(LPRINT.GT.0) RETURN                                            02400023
      WRITE(LOGUN,120) NSTEP                                            02410040
 120  FORMAT(/' SPECTRA G(E1,E2),N,P,A,F   FACTOR 2.  NSTEP=',I3)       02420074
      WRITE(LOGUN,125) (SYMB(6),I=1,10)                                 02430040
 125  FORMAT(20(A1,4X))                                                 02440013
      FACTOR=2.                                                         02450002
      CALL OUT2(SPECG,KEG,1,1,0,0,FACTOR,10,0,0)                        02460014
      WRITE(LOGUN,127) SYMB(1)                                          02470040
C
C     I.Giese 2007
C     2.0D0
C
      IF(NOE2.EQ.0)CALL OUT2(SPECQ,KEG,1,1,0,0,2.0D0,10,0,0)            02480014
      CALL OUT2(SPECN,KEN,1,0,0,0,FACTOR,10,0,0)                        02490002
      WRITE(LOGUN,127) SYMB(2)                                          02500040
      CALL OUT2(SPECP,KEP1,1,0,0,0,FACTOR,10,0,0)                       02510002
      WRITE(LOGUN,127) SYMB(3)                                          02520040
      CALL OUT2(SPECA,KEA1,1,0,0,0,FACTOR,10,0,0)                       02530002
      WRITE(LOGUN,127) SYMB(4)                                          02540040
      IF(NOF.EQ.0)   CALL OUT2(SPECF,KEF,1,0,0,0,FACTOR ,10,0,0)        02550002
      IF(NOF.EQ.0)   WRITE(LOGUN,127) SYMB(5)                           02560040
 127  FORMAT(1H+,55X,A2)                                                02570014
C                                                                       02580002
      IDUM=MAA-1                                                        02590017
      WRITE(LOGUN,140) EXMAX2,NSTEP,IDUM,ELEMNT(MZ),MAA,ELEMNT(MZ)      02600040
 140  FORMAT(/6X,'POPN   FACTOR=2.',6X,'EXMAX=',F8.1,6X,'NSTEP=',I3,    02610017
     1       4X,'DAUGHTER',I4,A4,4X,'PARENT',I4,A4)                     02620034
      M0=-EXMAX2                                                        02630002
      CALL OUT2(POPN,M2,N2,M0,0,1,FACTOR,10,1,0)                        02640002
C                                                                       02650002
      IF(N1.GT.M1) GOTO 162                                             02660002
      IF(M1.GT.127) GOTO 180                                            02670016
      WRITE(LOGUN,145) EXMAX1,NSTEP,MAA,ELEMNT(MZ)                      02680040
 145  FORMAT(/6X,'EMISSION TABLE    EXCIT HOR',4X,'EXMAX=',F8.1,        02690014
     1        6X,'NSTEP=',I3,4X,'PARENT',I4,A4)                         02700015
      WRITE(LOGUN,150) (SYMB(6),I=1,20)                                 02710040
 150  FORMAT(5X,20(A1,4X))                                              02720016
      IREP=0                                                            02730033
      DO 160 J=1,N1                                                     02740002
      JJ=    J*M1+1                                                     02750013
      J0=J-1                                                            02760002
      WRITE(LOGUN,170) J0,(TABLE(JJ-K),K=1,M1)                          02770040
      DO 155 K=1,M1                                                     02780028
      IF(TABLE(JJ-K).NE.BLANK) GOTO 160                                 02790028
 155  CONTINUE                                                          02800028
      IREP=IREP+1                                                       02810033
      IF(IREP.GE.2)GOTO 180                                             02820033
 160  CONTINUE                                                          02830028
      GOTO 180                                                          02840002
 162  JMAX=N1-1                                                         02850013
      IF(N1.GT.127) GO TO 180                                           02860016
      WRITE(LOGUN,164) JMAX,NSTEP,MAA,ELEMNT(MZ)                        02870040
 164  FORMAT(/6X,'EMISSION TABLE    SPINS HOR',4X,'JMAX=',I4,4X,        02880034
     1 'NSTEP=',I3,4X,'PARENT'  ,I4,A4)                                 02890035
      WRITE(LOGUN,125) (SYMB(6),I=1,20)                                 02900040
      IREP=0                                                            02910033
      DO 168 K=1,M1                                                     02920028
      K0=EXMAX1+2-K                                                     02930018
      WRITE(LOGUN,170) K0,(TABLE(K+(J -1)*M1),J=1,N1)                   02940040
      I=K-M1                                                            02950028
      DO 165 J=1,N1                                                     02960028
      I=I+M1                                                            02970028
      IF(TABLE(I).NE.BLANK) GOTO 167                                    02980034
 165  CONTINUE                                                          02990028
      IREP=IREP+1                                                       03000033
      IF(IREP.GE.2)GOTO 180                                             03010033
      GOTO 168                                                          03020034
 167  IREP=0                                                            03030034
 168  CONTINUE                                                          03040028
 170  FORMAT(I4,1X,127A1)                                               03050016
 180  CONTINUE                                                          03060002
      RETURN                                                            03070002
      END                                                               03080000
