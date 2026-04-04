      SUBROUTINE ENGOUT(LOGUN,IDISC,EXCIT0)                             00020028
C
C  Steuert  Ausgabe
C  Stand: 3.5.1994
C
      INCLUDE 'common.f'                                                00030034
C                                                                       00040024
      DIMENSION SIGMAS(50)                                              00050021
C     DIMENSION IS NAUP+NZUP                                            00060021
C     MAYBE UNINITIALIZED: J1
C
C     I.Giese 2007
      SAVE
C
      NEXC=NEXC+1                                                       00070009
                    LOGUN=6                                             00080015
      IF(LPRINT. GT.4)LOGUN=22                                          00090014
      IACN =ACN                                                         00100002
      IZCN =ZCN                                                         00110002
      I0=IACN-(IACN/10)*10+1                                            00120003
      DO 20 I=1,10                                                      00130003
      I1=I0-I                                                           00140003
      IF(I1.LT.0)I1=I1+10                                               00150003
 20   IND(I)=I1                                                         00160003
      NANZ=NA+NZ                                                        00170004
      DO 30 N=1,NANZ                                                    00180002
 30   SIGMAS(N)=0.                                                      00190002
      DO 40  NROW=1,NZ                                                  00200002
      DO 35  NSTEP=1,NA                                                 00210002
      N=NSTEP+NROW-1                                                    00220003
 35   SIGMAS(N)=SIGMAS(N)+XSECT(NROW,NSTEP)                             00230002
 40   CONTINUE                                                          00240002
C---------------------------------------------------------              00250030
      WRITE (LOGUN,50)                                                  00260002
 50   FORMAT(/   ' HIVAP  ENGOUT'/)                                     00270017
      WRITE(LOGUN,55)(TITLE(I),I=1,18)                                  00280002
 55   FORMAT(18A4)                                                      00290002
      IF(ELAB1.LE.0.) GOTO 58                                           00300017
      WRITE(LOGUN,57) ELAB1,ELAB,EXCIT1,EXCIT0,IACN,ELEMNT(IZCN)        00310028
 57   FORMAT(' AVERAGED CALCULATION  RANGE(LAB):',2F8.1,4X,'OR(EXCIT)', 00320017
     1        2F8.1/' CN=',I4,A4)                                       00330017
      GOTO 59                                                           00340017
 58   WRITE (LOGUN,60) ELAB,EXCIT0,IACN,ELEMNT(IZCN)                    00350028
 59   WRITE (LOGUN,65) SIGFUS,FLGRAZ,SIGEVA,SIGFIS                      00360017
 60   FORMAT(4X,'ELAB=',F8.1,6X,'EXCIT=',F8.1,6X,'CN=',I4,A4)           00370002
 65   FORMAT(4X,'FUSION X-SECT=',E12.4,6X,'LCRIT=',F6.1/                00380002
     1 4X,'EVAP SO FAR',E12.4,6X,'FISSION',E12.4)                       00390002
C-----------------------------------------------------                  00400030
      ELABE(NEXC)=ELAB                                                  00410009
      ELABE1(NEXC)=ELAB1                                                00420018
      SIGFUE(NEXC)=SIGFUS                                               00430009
      SIGEVE(NEXC)=SIGEVA                                               00440009
      SIGFIE(NEXC)=SIGFIS                                               00450009
      FLGRAE(NEXC)=FLGRAZ                                               00460009
C--------------------------------------------------------               00470030
      WRITE (LOGUN,105)                                                 00480002
 105  FORMAT(/' MASS DISTRIBUTION')                                     00490002
      WRITE(LOGUN,106)(IND(I),I=1,10)                                   00500003
 106  FORMAT(I9,4I12,6X,5I12)                                           00510003
      DO 120 N=1,NANZ,10                                                00520002
      N1=MIN0(N+9,NANZ)                                                 00530002
      MAS=IACN -N +1                                                    00540002
      WRITE (LOGUN,115) MAS,(SIGMAS(N2),N2=N,N1)                        00550002
 115  FORMAT(I4,5E12.4,6X,5E12.4)                                       00560002
 120  CONTINUE                                                          00570002
C--------------------------------------------------------               00580030
      WRITE(LOGUN,125)                                                  00590002
 125  FORMAT(/' Z-DISTIBUTION')                                         00600002
      NZ0=MIN0(NZ,12)                                                   00610002
      MZ0=ZCN-NZ                                                        00620002
      WRITE(LOGUN,130) (ELEMNT(MZ0+NROW),NROW=1,NZ0)                    00630002
 130  FORMAT(12(4X,A4,3X))                                              00640002
      WRITE (LOGUN,135) (SIGZ(NZ+1-NROW),NROW=1,NZ0)                    00650002
 135  FORMAT(12E11.3)                                                   00660002
C--------------------------------------------------------               00670030
      WRITE (LOGUN,150)                                                 00680002
 150  FORMAT(/'  CROSS SECTIONS (MB)')                                  00690002
      NA0=MIN0(NA,12)                                                   00700002
      ZAV=0.                                                            00710022
      AAV=0.                                                            00720022
      ZDUM=0.                                                           00730022
      ADUM=0.                                                           00740022
      SUM=0.                                                            00750022
      DO 170 NROW=1,NZ                                                  00760002
      MZ=IZCN -NROW+1                                                   00770002
      ZZ=MZ                                                             00780022
      IDUM=IACN -NROW+2                                                 00790002
      DO 155 NSTEP=1,NA                                                 00800002
      XSECTE(NROW,NSTEP,NEXC)=XSECT(NROW,NSTEP)                         00810009
      IND(NSTEP)=IDUM-NSTEP                                             00820022
      AA=IND(NSTEP)                                                     00830022
      XS=XSECT(NROW,NSTEP)                                              00840022
      ZAV=ZAV+ZZ*XS                                                     00850022
      AAV=AAV+AA*XS                                                     00860022
      ZDUM=ZDUM+ZZ*ZZ*XS                                                00870022
      ADUM=ADUM+AA*AA*XS                                                00880022
      SUM=SUM+XS                                                        00890022
 155  CONTINUE                                                          00900022
      WRITE(LOGUN,160)(IND(NA+1-NSTEP),ELEMNT(MZ),NSTEP=1,NA0)          00910002
      WRITE(LOGUN,165)(XSECT(NROW,NA+1-NSTEP),NSTEP=1,NA0)              00920002
 160  FORMAT(12(I5,A4,2X))                                              00930002
 165  FORMAT(12E11.3)                                                   00940002
 170  CONTINUE                                                          00950002
C----------------------------------------------------------             00960030
      IF(SUM.LE.0.) GOTO 176                                            00970022
      ZAV=ZAV/SUM                                                       00980022
      AAV=AAV/SUM                                                       00990022
      ZDUM=ZDUM/SUM                                                     01000022
      ADUM=ADUM/SUM                                                     01010022
      DELTAZ=IZCN-ZAV                                                   01020023
      DELTAA=IACN-AAV                                                   01030023
      FWHMZ=0.                                                          01040023
      FWHMA=0.                                                          01050023
      DUMZ=ZDUM-ZAV*ZAV                                                 01060026
      DUMA=ADUM-AAV*AAV                                                 01070026
      IF(DUMZ.LT.0. .OR. DUMA.LT.0.) GOTO 300                           01080026
      IF(DELTAZ.GT.0.01) FWHMZ=2.36*SQRT(ZDUM-ZAV*ZAV)                  01090026
      IF(DELTAA.GT.0.01) FWHMA=2.36*SQRT(ADUM-AAV*AAV)                  01100023
      WRITE(LOGUN,172) ZAV,DELTAZ,FWHMZ,AAV,DELTAA,FWHMA                01110023
 172  FORMAT(/' AVERAGE Z,DELTAZ,FWHMZ SO FAR', 3F8.2/                  01120022
     1        ' AVERAGE A,DELTAA,FWHMA SO FAR', 3F8.2/)                 01130022
C------------------------------------------------------------           01140030
      ROUND=DEL(1)-AINT(DEL(1))                                         01150035
      WRITE(LOGUN,311) SYMB(1)                                          01160035
C
C     I.Giese 2007
C     1.0D0
C
      CALL OUT1(SPECGT,20,ROUND,1.0D0,LOGUN)                            01170035
      IF(NOE2.EQ.1) GOTO 304                                            01180035
      ROUND=DEL(1)-AINT(DEL(1))                                         01190035
      WRITE(LOGUN,310) SYMB(1)                                          01200035
      CALL OUT1(SPECQT,20,ROUND,1.0D0,LOGUN)                            01210035
 304  ROUND=DEL(2)-AINT(DEL(2))                                         01220035
      WRITE(LOGUN,310) SYMB(2)                                          01230035
      CALL OUT1(SPECNT,25,ROUND,1.0D0,LOGUN)                            01240035
      ROUND=DEL(3)-AINT(DEL(3))                                         01250035
      WRITE(LOGUN,310) SYMB(3)                                          01260035
      CALL OUT1(SPECPT,45,ROUND,1.0D0,LOGUN)                            01270035
      ROUND=DEL(4)-AINT(DEL(4))                                         01280035
      WRITE(LOGUN,310) SYMB(4)                                          01290035
      CALL OUT1(SPECAT,80,ROUND,1.0D0,LOGUN)                            01300035
      IF(NOF.EQ.1) GOTO 176                                             01310036
      ROUND=DEL(5)-AINT(DEL(5))                                         01320035
      WRITE(LOGUN,310) SYMB(5)                                          01330035
      CALL OUT1(SPECFT,20,ROUND,1.0D0,LOGUN)                            01340035
 310  FORMAT(' SPECTRUM ',A1)                                           01350035
 311  FORMAT(/' SPECTRUM ',A1)                                          01360035
C------------------------------------------------------------           01370035
 176  IF(IOPT.NE.0) GOTO 180                                            01380022
      WRITE(LOGUN,173)                                                  01390012
 173  FORMAT(/' HALFLIFE FIRST STEP IN UNITS OF 10**(-22) SEC AS FUNTION01400012
     1 OF J')                                                           01410012
      WRITE(LOGUN,175)(FLAND(I),I=1,JUPYR)                              01420012
 175  FORMAT(10E12.3)                                                   01430012
C----------------------------------------------------------             01440030
 180  IF(NUMISO.EQ.0) GOTO 200                                          01450030
      WRITE(LOGUN,185)                                                  01460012
 185  FORMAT(/' ISOMER CROSS SECTIONS')                                 01470012
      DO 190 ISO=1,NUMISO                                               01480012
      MZ=IZCN-ISOZ(ISO)+1                                               01490011
      IDUM=IACN-ISOZ(ISO)-ISOA(ISO)+2                                   01500011
      WRITE(LOGUN,187) IDUM,ELEMNT(MZ),SIGISO(ISO)                      01510012
 187  FORMAT(I5,A4,2X,E12.4)                                            01520012
 190  CONTINUE                                                          01530012
C----------------------------------------------------------             01540030
 200  IF(JFIS.NE.0) THEN                                                01550037
         DO 205 J=1,JDIM                                                01560037
         J1=JDIM+1-J                                                    01570037
         IF(YFISS(J1).GT.0.) GOTO 206                                   01580037
 205     CONTINUE                                                       01590037
 206     N1MX=J1                                                        01600037
         IF(JFIS.EQ.5) THEN                                             01610037
            DO 208 J=1,N1MX                                             01620037
 208        YFISS(J)=1.00000-YFISS(J)                                   01630037
            WRITE(19,210) ELAB,N1MX,JFIS                                01640037
 210        FORMAT('C 1-YFISS',4X,'ELAB=',F8.3,I6,I4/'H: X-1    Y')     01650037
         ELSE                                                           01660037
            WRITE(19,211) ELAB,N1MX,JFIS                                01670037
 211        FORMAT('C YFISS',4X,'ELAB=',F8.3,I6,I4/'H: X-1    Y')       01680037
         ENDIF                                                          01690037
         WRITE(19,215) (J,YFISS(J),J=1,N1MX)                            01700037
 215     FORMAT(4(I4,E12.4,';'))                                        01710037
      ENDIF                                                             01720037
C----------------------------------------------------------             01730030
 250  IF(NEXC.GE.NENMX) CALL CASOUT(LOGUN,IDISC)                        01740039
      RETURN                                                            01750002
C----------------------------------------------------------             01760030
 300  WRITE(LOGUN,301) ZDUM,ADUM,ZAV,AAV ,SUM                           01770026
 301  FORMAT(' ENGOUT NEG SQRT  ZDUM,ADUM,ZAV,AAV=',4E12.4,' SUM=',E12.401780026
     1      )                                                           01790026
      GOTO 176                                                          01800026
      END                                                               01810002
