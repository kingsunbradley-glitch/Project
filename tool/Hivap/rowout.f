      SUBROUTINE  ROWOUT(NROW,LOGUN,SIGMXZ,NOPUT,EXCIT0)                00020053
C
C  Steuert Datenausgabe
C  Stand: 3.5.1994
C
      INCLUDE 'common.f'                                                00030056
C
C     I.Giese 2007
      SAVE
C
                     LOGUN=6                                            00040038
      IF(LPRINT.GT.3)LOGUN=22                                           00050037
      SIGZ(NROW)=0.                                                     00060024
      DO 30 NSTEP=1,NA                                                  00070022
      XSECT(NROW,NSTEP)=SIG(NSTEP)                                      00080024
      SIGFIS=SIGFIS+RLOSTF(NSTEP)                                       00090024
      SIGEVA=SIGEVA+SIG(NSTEP)                                          00100024
 30   SIGZ(NROW)=SIGZ(NROW)+SIG(NSTEP)                                  00110024
      SIGMXZ=DMAX1(SIGMXZ,SIGZ(NROW))                                   00120024
C                                                                       00130005
      IF(LPRINT.GE.3) THEN                                              00140061
         DO 31 K=1,KEG                                                  00150061
         SPECGT(K)=SPECGT(K)+SPECG(K)                                   00160061
 31      SPECQT(K)=SPECQT(K)+SPECQ(K)                                   00170061
         DO 32 K=1,KEN                                                  00180061
 32      SPECNT(K)=SPECNT(K)+ SPECN(K)                                  00190061
         DO 33 K=1,KEP1                                                 00200061
 33      SPECPT(K)=SPECPT(K)+SPECP(K)                                   00210061
         DO 34 K=1,KEA1                                                 00220061
 34      SPECAT(K)=SPECAT(K)+SPECA(K)                                   00230061
         DO 35 K=1,KEF                                                  00240061
 35      SPECFT(K)=SPECFT(K)+SPECF(K)                                   00250061
      ENDIF                                                             00260061
C                                                                       00270059
      MZ=KZ(NROW)                                                       00280006
      IF(LPRINT.GT.3)GOTO 116                                           00290028
      WRITE(LOGUN,50) NROW,ELEMNT(MZ),SIGZ(NROW)                        00300024
 50   FORMAT(1H1,4X,'ROW',I3,6X,'ELEMENT ',A4,4X,'X-SECTION',E12.3)     00310020
 55   WRITE (LOGUN,101)(TITLE(I),I=1,18)                                00320022
 101  FORMAT   ( 18A4)                                                  00330005
      WRITE(LOGUN,110) ELAB,QC,EXCIT0                                   00340053
 110  FORMAT(' ELAB(MEV)',F9.3,5X,'Q-VALUE' ,F9.3,4X,'EXCIT',F9.3)      00350005
      WRITE (LOGUN,114) SIGFUS,SIGEVA,SIGFIS,SIGMXZ                     00360024
 114  FORMAT(' SIGFUS,SIGEVA,SIGFIS,SIGMXZ',4E12.4)                     00370024
C     FOR PLT02                                                         00380005
 116                  GOTO 160                                          00390032
 120  FORMAT( 9F8.3)                                                    00400005
 130  FORMAT(18I4)                                                      00410005
C                                                                       00420005
 160  IF(LPRINT.GT.3) GOTO 340                                          00430028
      WRITE(LOGUN,170)IZT,MTGT,IZP,MPROJ                                00440027
 170  FORMAT(' ZT=',I3,4X,'AT',I4,4X,'ZP',I3,4X,'AP',I4//)              00450006
C---------------------------------------------------------------        00460041
      NA1=1                                                             00470042
      NA2=MIN0(12,NA)                                                   00480042
 171  DO 175 NSTEP=1,NA                                                 00490041
 175  IND(NSTEP)=MZ+IN(NSTEP)                                           00500006
      WRITE(LOGUN,176) (IND(NSTEP),ELEMNT(MZ),NSTEP=NA1,NA2)            00510041
 176  FORMAT(11X,12(2X,I3,A4,1X))                                       00520034
      WRITE(LOGUN,180)(SIG(NSTEP) ,NSTEP=NA1 ,NA2)                      00530041
 180  FORMAT(' SIG',7X,12E10.3)                                         00540005
      DO 182 NSTEP=NA1,NA2                                              00550041
      IND(NSTEP)=0                                                      00560033
      IF(TOTAL(NSTEP).GT.0.)IND(NSTEP)=SIG  (NSTEP)*10000./TOTAL(NSTEP) 00570033
     1                                 +.5                              00580033
 182  CONTINUE                                                          00590033
      WRITE(LOGUN,183)(IND(NSTEP) ,NSTEP=NA1 ,NA2)                      00600041
 183  FORMAT(' B-G',6X,12I10  )                                         00610034
      DO 185 NSTEP=NA1,NA2                                              00620041
      IND(NSTEP)=0                                                      00630033
      IF(TOTAL(NSTEP).GT.0.)IND(NSTEP)=SIGNN(NSTEP)*10000./TOTAL(NSTEP) 00640033
     1                                 +.5                              00650033
 185  CONTINUE                                                          00660033
      WRITE(LOGUN, 190)(IND  (NSTEP),NSTEP=NA1 ,NA2)                    00670041
 190  FORMAT('  R-N',5X,12I10)                                          00680058
      DO 195 NSTEP=NA1,NA2                                              00690041
      IND(NSTEP)=0                                                      00700033
      IF(TOTAL(NSTEP).GT.0.)IND(NSTEP)=SIGP (NSTEP)*10000./TOTAL(NSTEP) 00710033
     1                                 +.5                              00720033
 195  CONTINUE                                                          00730033
      WRITE(LOGUN,200)( IND (NSTEP),NSTEP=NA1,NA2)                      00740041
 200  FORMAT(' BR-P',5X,12I10  )                                        00750034
      DO 205 NSTEP=NA1,NA2                                              00760041
      IND(NSTEP)=0                                                      00770033
      IF(TOTAL(NSTEP).GT.0.)IND(NSTEP)=SIGA (NSTEP)*10000./TOTAL(NSTEP) 00780033
     1                                 +.5                              00790033
 205  CONTINUE                                                          00800033
      WRITE(LOGUN,210)( IND (NSTEP),NSTEP=NA1,NA2)                      00810041
 210  FORMAT(' BR-A'  ,5X,12I10  )                                      00820034
      DO 209 NSTEP=NA1,NA2                                              00830041
      IND(NSTEP)=0                                                      00840033
      IF(TOTAL(NSTEP).GT.0.)IND(NSTEP)=SIGF (NSTEP)*10000./TOTAL(NSTEP) 00850033
     1                                 +.5                              00860033
 209  CONTINUE                                                          00870033
      IF(NOF.EQ.0) WRITE(LOGUN,211) (IND   (NSTEP),NSTEP=NA1,NA2)       00880041
 211  FORMAT(' BR-F'   ,5X,12I10  )                                     00890034
      IF(NOF.EQ.0) WRITE(LOGUN,212) (RLOSTF(NSTEP),NSTEP=NA1,NA2)       00900041
 212  FORMAT(' SIGF'   ,6X,12E10.3)                                     00910033
      WRITE(LOGUN,230)( TOTAL(NSTEP),NSTEP=NA1,NA2)                     00920041
 230  FORMAT(' TOTAL',5X,12E10.3)                                       00930005
      WRITE(LOGUN,231)(  AVGE(NSTEP),NSTEP=NA1,NA2)                     00940041
 231  FORMAT(' EXCP',5X,12F10.1)                                        00950033
      WRITE(LOGUN, 232)( AVGJ(NSTEP),NSTEP=NA1,NA2)                     00960041
 232  FORMAT(' SPINP',4X,12F10.1/)                                      00970033
      WRITE(LOGUN,233)( AVGEF(NSTEP),NSTEP=NA1,NA2)                     00980041
 233  FORMAT(' EXCF',5X,12F10.1)                                        00990033
      WRITE(LOGUN, 234)(AVGJF(NSTEP),NSTEP=NA1,NA2)                     01000041
 234  FORMAT(' SPINF',4X,12F10.1/)                                      01010033
      WRITE(LOGUN,235)( AVDJN(NSTEP),NSTEP=NA1,NA2)                     01020041
 235  FORMAT(' DJN'  ,6X,12F10.2)                                       01030033
      WRITE(LOGUN,236)( AVDJP(NSTEP),NSTEP=NA1,NA2)                     01040041
 236  FORMAT(' DJP'  ,6X,12F10.2)                                       01050033
      WRITE(LOGUN,237)( AVDJA(NSTEP),NSTEP=NA1,NA2)                     01060041
 237  FORMAT(' DJA'  ,6X,12F10.2)                                       01070033
      WRITE(LOGUN,238)( EAVN (NSTEP),NSTEP=NA1,NA2)                     01080041
 238  FORMAT(' AV.E N',3X,12F10.2)                                      01090057
      WRITE(LOGUN,239)( EAVP (NSTEP),NSTEP=NA1,NA2)                     01100041
 239  FORMAT(' AV.E P',3X,12F10.2)                                      01110057
      WRITE(LOGUN,240)( EAVA (NSTEP),NSTEP=NA1,NA2)                     01120041
 240  FORMAT(' AV.E A',3X,12F10.2)                                      01130057
      WRITE(LOGUN,251)(SUMEG (NSTEP),NSTEP=NA1,NA2)                     01140045
 251  FORMAT(' SUM GAM',2X,12F10.1)                                     01150046
      WRITE(LOGUN,252)(SUMJG (NSTEP),NSTEP=NA1,NA2)                     01160049
 252  FORMAT(' SPIN GAM',1X,12F10.1)                                    01170051
C                                                                       01180045
      WRITE(LOGUN,241) (TRIM(NSTEP),NSTEP=NA1,NA2)                      01190041
 241  FORMAT(/ ' TRIM',6X,12E10.3)                                      01200045
      WRITE(LOGUN,245) (FLOST (NSTEP),NSTEP=NA1,NA2)                    01210041
 245  FORMAT(' FLOST ',4X,12E10.3)                                      01220057
      WRITE(LOGUN,250) (RP(NSTEP),NSTEP=NA1,NA2)                        01230041
 250  FORMAT(' RP',8X,12E10.3)                                          01240006
      WRITE(LOGUN,260) (RLOSTG(NSTEP),NSTEP=NA1,NA2)                    01250041
 260  FORMAT(' RLOSTG',4X,12E10.3)                                      01260057
      WRITE(LOGUN,270) (RLOSTN(NSTEP),NSTEP=NA1,NA2)                    01270041
 270  FORMAT(' RLOSTN',4X,12E10.3)                                      01280057
      WRITE(LOGUN,280) (RLOSTP(NSTEP),NSTEP=NA1,NA2)                    01290041
 280  FORMAT(' RLOSTP',4X,12E10.3)                                      01300057
      WRITE(LOGUN,290) (RLOSTA(NSTEP),NSTEP=NA1,NA2)                    01310041
 290  FORMAT(' RLOSTA',4X,12E10.3)                                      01320057
      WRITE(LOGUN,295) (IME   (NSTEP),NSTEP=NA1,NA2)                    01330041
 295  FORMAT(' TIME', 5X,12I10//)                                       01340042
C-----------------------------------------------------------------      01350041
      IF(NA2.EQ.NA) GOTO 303                                            01360042
      NA1=NA1+12                                                        01370042
      NA2=NA2+12                                                        01380042
      NA2=MIN0(NA,NA2)                                                  01390042
      GOTO 171                                                          01400041
 303  IF(LPRINT.LT.3) GOTO 340                                          01410041
      ROUND=DEL(1)-AINT(DEL(1))                                         01420031
      WRITE(LOGUN,311) SYMB(1)                                          01430031
C
C     I.Giese 2007
C     1.0D0
C
      CALL OUT1(SPECG,20,ROUND,1.0D0,LOGUN)                             01440031
      IF(NOE2.EQ.1) GOTO 304                                            01450031
      ROUND=DEL(1)-AINT(DEL(1))                                         01460047
      WRITE(LOGUN,310) SYMB(1)                                          01470031
      CALL OUT1(SPECQ,20,ROUND,1.0D0,LOGUN)                             01480031
 304  ROUND=DEL(2)-AINT(DEL(2))                                         01490031
      WRITE(LOGUN,310) SYMB(2)                                          01500031
      CALL OUT1(SPECN,25,ROUND,1.0D0,LOGUN)                             01510031
      ROUND=DEL(3)-AINT(DEL(3))                                         01520031
      WRITE(LOGUN,310) SYMB(3)                                          01530031
      CALL OUT1(SPECP,45,ROUND,1.0D0,LOGUN)                             01540031
      ROUND=DEL(4)-AINT(DEL(4))                                         01550031
      WRITE(LOGUN,310) SYMB(4)                                          01560031
      CALL OUT1(SPECA,80,ROUND,1.0D0,LOGUN)                             01570031
      IF(NOF.EQ.1) GOTO 340                                             01580031
      ROUND=DEL(5)-AINT(DEL(5))                                         01590031
      WRITE(LOGUN,310) SYMB(5)                                          01600031
      CALL OUT1(SPECF,20,ROUND,1.0D0,LOGUN)                             01610031
 310  FORMAT(' SPECTRUM ',A1)                                           01620056
 311  FORMAT(/' SPECTRUM ',A1)                                          01630056
 340  IF(LPRINT.GT.2) GOTO 999                                          01640062
      WRITE (LOGUN,350)                                                 01650021
 350  FORMAT(/' STORAGE')                                               01660021
      WRITE(LOGUN,355)(STORE(NROW+1,NSTEP),NSTEP=1,NA)                  01670024
      WRITE(LOGUN,355)(STORE(NROW+2,NSTEP),NSTEP=1,NA)                  01680024
 355  FORMAT(15F4.0)                                                    01690024
 600  FORMAT(6E12.4)                                                    01700005
 999  IF(NOPUT.EQ.1) GOTO 9999                                          01710025
C-------------------------------------------------------------------    01720024
 1001 FORMAT(18I4)                                                      01730024
 1002 FORMAT(3F8.2,4E12.4)                                              01740024
 1003 FORMAT(9F6.1)                                                     01750024
 1004 FORMAT(6E12.4)                                                    01760024
 9999 RETURN                                                            01770040
      END                                                               01780000
