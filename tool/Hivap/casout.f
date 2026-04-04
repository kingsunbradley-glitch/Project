      SUBROUTINE CASOUT(LOGUN,IDISC)                                    00010012
C 
C  Steuerung der Datenausgabe
C  Stand: 5. 5. 1994
C
      INCLUDE 'common.f'                                                00030033
      COMMON/DELFIS/EDELFIS,ADELFIS,INOF
      DIMENSION XSUM(30),EXCITE(NENMX),EXCITE1(NENMX)                   00040037
      DIMENSION ZAV(NENMX),AAV(NENMX),SUMZN(NENMX),ZSTDEV(NENMX)        00050034
      DIMENSION ASTDEV(NENMX),SUMN(NENMX,27)                            00060034
C ------------------------------------------------------------          00070033
C
C     I.Giese 2007
      SAVE
C
      AMPAR1=0.
      LOGUN=6                                                           00080033
      IF(NEXC.EQ.0) RETURN                                              00090018
      AP=MPROJ                                                          00100027
      AT=MTGT                                                           00110027
      WRITE(LOGUN,50)                                                   00120002
 50   FORMAT(/' HIVAP  CASOUT')                                         00130002
      WRITE(LOGUN,55)(TITLE(I),I=1,18)                                  00140002
 55   FORMAT(/18A4/)                                                    00150002
      IF(INOF.NE.1) GOTO 57
      WRITE(LOGUN,56) EDELFIS,ADELFIS
 56   FORMAT(/2X,'Rechnung mit Behinderung der Spaltung; Parameter:',
     11X,'E_Beh=',F5.1,2X,'DE_Beh=',F5.1)      
 57   WRITE(LOGUN,58)                                                   00160002
 58   FORMAT(/5X,'E_lab  A*MeV   E*/MeV',5X,'SI_fus',6X,'SI_VR ',4X,    00170002
     1       'SI_Spalt  L_krit',5X,'b_Spalt'/)                          00180004
C                                                                       00190033
      DO 70 I=1,NEXC                                                    00200002
      ZAV(I)=0.                                                         00210033
      AAV(I)=0.                                                         00220033
      SUMZN(I)=0.                                                       00230034
      DO NROW=1,NZ                                                      00240036
        SUMN(I,NROW)=0.                                                 00250036
      ENDDO                                                             00260034
      ZSTDEV(I)=0.                                                      00270033
      ASTDEV(I)=0.                                                      00280033
      DUM=ELABE(I)*FLOAT(MTGT )/ACN +  QC                               00290003
      DUM0=ELABE1(I)*FLOAT(MTGT )/ACN +  QC                             00300010
      DUM1=ELABE(I)/FLOAT(MPROJ)                                        00310002
      DUM11=ELABE1(I)/FLOAT(MPROJ)                                      00320010
      DUM2=0.                                                           00330004
      IF(SIGFUE(I).GT.0.) DUM2=SIGFIE(I)*100./SIGFUE(I)                 00340004
      WRITE(LOGUN,60)ELABE(I),DUM1,DUM,SIGFUE(I),SIGEVE(I),SIGFIE(I),   00350002
     1               FLGRAE(I),DUM2                                     00360004
      IF(ELABE1(I).GT.0.)WRITE(LOGUN,60) ELABE1(I),DUM11,DUM0           00370010
      EXCITE(I)=DUM                                                     00380037
C
C     I.Giese 2007
      EXCITE1(I)=0.0
C
      IF(ELABE1(I).GT.0.) EXCITE1(I)=DUM0                               00390037
 60   FORMAT(F10.3,F8.2,F8.1,3E12.4,F8.1,F12.2)                         00400004
 70   CONTINUE                                                          00410002
C---------------------------------------------                          00420002
      IZCN=ZCN+0.01                                                     00430002
      IACN=ACN+0.01                                                     00440002
      WRITE (LOGUN,745)
745   Format(/'     Querschnitte / mbarn')                              00450033
      NA0=MIN0(NA,12)                                                   00460002
C                                                                       00470033
      DO 170 NROW=1,NZ                                                  00480002
      MZ=IZCN -NROW+1                                                   00490002
      IDUM=IACN -NROW+2                                                 00500002
      ZZ=MZ-IZCN                                                        00510038
C                                                                       00520033
      DO 155 NSTEP=1,NA                                                 00530002
      IND(NSTEP)=IDUM-NSTEP                                             00540033
      AA=IND(NSTEP)-IACN                                                00550038
      DO IE =1,NEXC                                                     00560033
        XS=XSECTE(NROW,NSTEP,IE)                                        00570034
        ZAV(IE)=ZAV(IE)+ZZ*XS                                           00580034
        AAV(IE)=AAV(IE)+AA*XS                                           00590034
        ZSTDEV(IE)=ZSTDEV(IE)+ZZ*ZZ*XS                                  00600034
        ASTDEV(IE)=ASTDEV(IE)+AA*AA*XS                                  00610034
        SUMZN(IE)=SUMZN(IE)+XS                                          00620034
        SUMN(IE,NROW)=SUMN(IE,NROW)+XS                                  00630034
      ENDDO                     ! IE                                    00640034
 155  CONTINUE                  ! NSTEP                                 00650033
C                                                                       00660034
      WRITE(LOGUN,'(/''  E*/MeV '',12(I5,A4,1X))')                       00670033
     1            (IND(NA+1-NSTEP),ELEMNT(MZ),NSTEP=1,NA0)              00680033
      DO     IE=1,NEXC                                                  00690033
         WRITE(LOGUN,'(F8.1,12E10.3)')                                  00700033
     1            EXCITE(IE),(XSECTE(NROW,NA+1-NSTEP,IE),NSTEP=1,NA0)   00710037
         IF(EXCITE1(IE).GT.0.) WRITE(LOGUN,'(F8.1)') EXCITE1(IE)        00720037
      ENDDO   ! IE                                                      00730033
 170  CONTINUE                                                          00740002
C                                                                       00750033
      IF(NZ.GT.2) THEN                                                  00760039
        WRITE(LOGUN,'(/'' Z-Verteilung  '')')                           00770039
        NZ0=MIN(NZ,12)                                                  00780039
        MZ0=ZCN-NZ                                                      00790039
        WRITE(LOGUN,'('' E*/MeV '',12(4X,A4,3X))')                      00800039
     1            (ELEMNT(MZ0+NROW),NROW=1,NZ0)                         00810034
        DO IE=1,NEXC                                                    00820039
         IF(SUMZN(IE).GT.0.) THEN                                       00830034
            WRITE(LOGUN,'(F7.1,E11.3,11E10.3)') EXCITE(IE),             00840038
     1           (SUMN(IE,NZ+1-NROW),NROW=1,NZ0)                        00850034
            IF(EXCITE1(IE).GT.0.) WRITE(LOGUN,'(F7.1)') EXCITE1(IE)     00860038
            ZAV(IE)=ZAV(IE)/SUMZN(IE)                                   00870034
            AAV(IE)=AAV(IE)/SUMZN(IE)                                   00880034
            ZSTDEV(IE)=ZSTDEV(IE)/SUMZN(IE)                             00890034
            DUM=DMAX1(AMPAR1,ZSTDEV(IE)-ZAV(IE)**2)                     00900034
            ZSTDEV(IE)=SQRT(DUM)                                        00910034
            ZAV(IE)=ZAV(IE)+ZCN                                         00920038
            ASTDEV(IE)=ASTDEV(IE)/SUMZN(IE)                             00930034
            DUM=DMAX1(AMPAR1,ASTDEV(IE)-AAV(IE)**2)                     00940034
            ASTDEV(IE)=SQRT(DUM)                                        00950034
            AAV(IE)=AAV(IE)+ACN                                         00960038
         ELSE                                                           00970034
            ZAV(IE)=0.                                                  00980034
            AAV(IE)=0.                                                  00990034
            ZSTDEV(IE)=0.                                               01000034
            ASTDEV(IE)=0.                                               01010034
         ENDIF                                                          01020034
        ENDDO             ! IE                                          01030039
C                                                                       01040034
        WRITE(LOGUN,'(/'' Z und A Mittelwerte    ''/                    01050039
     1''   E*/MeV     MW.Z  St.Abw.Z      MW.A  St.Abw.A'') ')          01060038
        DO IE=1,NEXC                                                    01070039
          WRITE(LOGUN,'(F8.1,4F10.2)')                                  01080039
     1       EXCITE(IE),ZAV(IE),ZSTDEV(IE),AAV(IE),ASTDEV(IE)           01090037
        ENDDO                                                           01100039
      ENDIF   ! NZ>3                                                    01110039
C ----------------------------------------------------------------      01120034
      IF(IDISC.EQ.0) GOTO 300                                           01130014
      IF(IDISC.EQ.2) GOTO 210                                           01140029
      DO 200 NROW=1,NZ                                                  01150012
      MZ=IZCN-NROW+1                                                    01160012
      DO 185 NSTEP=1,NA                                                 01170012
      DO 171 I=1,NEXC                                                   01180014
      IF(XSECTE(NROW,NSTEP,I).GT.0.) GOTO 174                           01190014
 171  CONTINUE                                                          01200014
      GOTO 185                                                          01210014
 174  MA=IACN-NROW-NSTEP+2                                              01220014
      IF(MA.GT.99) WRITE(25,172) MA,ELEMNT(MZ),MA,ELEMNT(MZ)            01230030
      IF(MA.LE.99) WRITE(25,173) MA,ELEMNT(MZ),MA,ELEMNT(MZ)            01240030
 172  FORMAT('H:  X  Y  ''',I3,A4,''''/'C ',I3,A4)                      01250030
 173  FORMAT('H:  X  Y  ''',I2,A4,''''/'C ',I2,A4)                      01260030
      WRITE(25,175) (EXCITE(I),XSECTE(NROW,NSTEP,I),I=1,NEXC)           01270037
 175  FORMAT(4(F6.1,E11.3,';'))                                         01280027
 185  CONTINUE                                                          01290012
 200  CONTINUE                                                          01300012
C-----------------------------                                          01310029
      IF(IDISC.EQ.3) GOTO 210                                           01320027
      GOTO 230                                                          01330029
 210  WRITE(25,212)                                                     01340029
 212  FORMAT('H:   X Y  ''SIGFUS'''/'C EXCIT VS SIGFUS(MB)')            01350031
      WRITE(25,175) (EXCITE(I), SIGFUE(I),I=1,NEXC)                     01360037
      WRITE(25,214)                                                     01370023
 214  FORMAT('H:   X Y  ''SIGER'''/'C EXCIT VS SIGEVA(MB)')             01380031
      WRITE(25,175)(EXCITE(I),SIGEVE(I),I=1,NEXC)                       01390037
      WRITE(25,217)                                                     01400022
 217  FORMAT('H:   X Y ''SIGFISS'''/'C EXCIT VS SIGFISS(MB)')           01410031
      WRITE(25,175)(EXCITE(I),SIGFIE(I),I=1,NEXC)                       01420037
C --------                                                              01430031
      DO 226 NROW=1,NZ                                                  01440031
      MZ=IZCN-NROW+1                                                    01450031
      WRITE(25,220) ELEMNT(MZ)                                          01460031
 220  FORMAT('H: X Y  ''',A4,'''')                                      01470031
      DO 224 I=1,NEXC                                                   01480031
         XSUM(I)=0.                                                     01490031
         DO 222 NSTEP=1,NA                                              01500031
         XSUM(I)=XSUM(I)+XSECTE(NROW,NSTEP,I)                           01510032
 222     CONTINUE                                                       01520031
 224  CONTINUE                                                          01530031
      WRITE(25,175) (EXCITE(I),XSUM(I),I=1,NEXC)                        01540037
 226  CONTINUE                                                          01550031
C-------------                                                          01560029
 230  DUM=FLOAT(MTGT)/ACN                                               01570029
      WRITE(25 ,231) QC,DUM,NEXC                                        01580029
 231  FORMAT('C   QC=',F8.2,4X,'AT/ACN=',F9.4,4X,'NEXC=',I4)            01590029
C                                                                       01600029
      DUM0=-QC                                                          01610029
      DUM1=1.                                                           01620029
      WRITE(25,232) DUM0,DUM1                                           01630029
 232  FORMAT('C RESC  X0',F9.3,2X,'X1',F4.0,6X,' FROM 1',               01640029
     1       4X,'SHIFT TO ECM')                                         01650029
C                                                                       01660029
      DUM0=-QC*ACN/AT                                                   01670029
      DUM1=ACN/AT                                                       01680029
      WRITE(25,234) DUM0,DUM1                                           01690029
 234  FORMAT('C RESC  X0',F9.3,2X,'X1',F10.5,' FROM 1',                 01700029
     1      4X,'SHIFT TO ELAB')                                         01710029
C                                                                       01720029
      DUM1=ACN/(AP*AT)                                                  01730029
      DUM0=DUM1*(-QC)                                                   01740029
      WRITE(25,236) DUM0,DUM1                                           01750029
 236  FORMAT('C RESC  X0',F9.3,2X,'X1',F10.6,' FROM 1',4X,'SHIFT TO MEV/01760029
     1U')                                                               01770029
 300  NEXC=0                                                            01780012
      RETURN                                                            01790002
      END                                                               01800000
