      SUBROUTINE OWKB(K,KEMAX,LMAX,TLE,APAR,ZPAR,TLIM,SIGLOW,TLOW,      00010000
     1                KDIM,Q2,IPRNT)                                    00020000
C    
C  Berechnet Transmissionskoeefizienten fuer n,p,alphas fuer Abdampfung
C  nach der WKB - Methode, ruft FUSION auf
C  Stand: 3.5.1994
C
C     CALCULATES TRANSMISSION COEFF. N,P,A FOR EVAP.THEORY WITH WKB     00040000
C     JE=1 CORRESPONDS TO TLIM(LOWEST ENERGY)                           00050000
C     K=1 NEUTRONS,=2 PROTONS,=3 ALPHAS                                 00060000
C     APAR,ZPAR PARENT NUCLEUS                                          00070000
C     KEMAX=MAXIMUM NUMBER OF ENERGIES FOR WHICH TRANSM.IS CALCULATED   00080000
C     LMAX=MAX.ORBITAL L FOR EMITTED N,P,A                              00090000
C---------------------------------------------------------------------  00100000
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TLE(1),SIG(100),SIGML(60)                               00110000
      COMMON/POTS/ VRN,R0RN,ADIFRN,VIN,RIN,ADIFIN,                      00120000
     1             VRP,R0RP,ADIFRP,VIP,RIP,ADIFIP,RCLMBP,CBFACP,        00130000
     2             VRA,R0RA,ADIFRA,VIA,RIA,ADIFIA,RCLMBA,CBFACA,QQ2     00140000
      COMMON/FUS/VBFUS,RBFUS,HWFUS                                      00150005
C
C     I.Giese 2007
      SAVE
C
      K=IABS(K)                                                         00160000
      HBAR=.65819E-21                                                   00170000
C
C     I.Giese 2007
C     Wert 1...E-44 nur fuer RISC6000-Architektur
C     PC/Linux: 1.1754994E-38 kleinste sinnvolle Zahl
C
      AMU=1.03641E-44                                                   00180000
      R0FUS=0.                                                          00190000
      CRED=1.                                                           00200000
      NORUTH=1                                                          00210000
      NOPROX=1                                                          00220000
      NOCURV=1                                                          00230000
      ION=1                                                             00240000
      ITEST=0                                                           00250000
      IOPT=-1                                                           00260000
C---------------------------                                            00270000
      GOTO (10,20,30),K                                                 00280000
C                                   NEUTRONS PEREY                      00290001
 10   AP=1.                                                             00300000
      ZP=0.                                                             00310000
      AT=APAR-AP                                                        00320000
      V0=VRN                                                            00330000
      R0=R0RN                                                           00340000
      D=ADIFRN                                                          00350000
      IF(V0.LE.0.) V0=47.01                                             00360002
      IF(R0.LE.0.) R0=1.322-7.6E-4* AT+4.E-6* AT* AT-8.E-9* AT**3.      00370000
      R0T=R0                                                            00380000
      IF(D.LE.0.)  D=0.660                                              00390000
      IF(IPRNT.NE.0) WRITE(6,15) V0,R0,D                                00400001
 15   FORMAT(/' OWKB NEUTRONS',4X,'V0=',F7.2,4X,'R0=',F7.3,4X,'D=',F7.3)00410001
      GOTO 35                                                           00420000
C-------------------PROTONS OWN FIT                                     00430001
 20   AP=1.                                                             00440000
      ZP=1.                                                             00450000
      AT=APAR-AP                                                        00460000
      V0=VRP                                                            00470000
      R0=R0RP                                                           00480000
      D=ADIFRP                                                          00490000
C     PEREY 1963 POT                                                    00500000
C     IF(V0.LE.0.) V0=53.3+27.*( AT-2.*ZT)/AT +0.4*ZT/AT**0.33333       00510001
C     IF(R0.LE.0.) R0=1.25                                              00520001
C     FOLLOWING PARAMETERS FOR TEST PURPOSES: FITTED TO DATA HOFMANN    00530001
C      V0=60                                                            00540000
C      R0=1.40                                                          00550000
C     WKB POT FITTED TO KURCEWICZ DATA TH232+P(1981)                    00560001
      IF(V0.LE.0.) V0=62                                                00570001
      IF(R0.LE.0.) R0=1.254                                             00580001
      R0T=R0                                                            00590000
      IF(D.LE.0.)  D=0.750                                              00600001
      IF(IPRNT.NE.0) WRITE(6,25)                                        00610000
 25   FORMAT(/' OWKB  PROTONS')                                         00620000
      GOTO 35                                                           00630000
C--ALPH-MODIFIED  SATCHLER (IGO USE V0=1100.;R0=1.17-1.6/AT3;D=0.574)   00640001
C             FITTED TO ALPHA DATA SEP 81                               00650000
 30   AP=4.                                                             00660000
      ZP=2.                                                             00670000
      AT =APAR-AP                                                       00680000
      AT3=AT**0.333333                                                  00690000
      V0=VRA                                                            00700000
      R0=R0RA                                                           00710000
      D=ADIFRA                                                          00720000
      IF(V0.LE.0.) V0=50.2                                              00730000
      IF(R0.LE.0.) R0=1.2067                                            00740000
      R0T=R0+1.6/AT3                                                    00750000
      IF(D.LE.0.)  D=0.564                                              00760000
      IF(IPRNT.NE.0) WRITE(6,34)                                        00770000
 34   FORMAT(/' OWKB  ALPHAS ')                                         00780000
C---------------------------                                            00790000
 35   Q2=QQ2                                                            00800000
C---------------------------                                            00810000
 45   ELAB=TLIM -1.                                                     00820000
      AT=APAR-AP                                                        00830000
      ZT=ZPAR-ZP                                                        00840000
      ARED=AMU*AP*AT/(AP+AT)                                            00850000
      HW0=HBAR*HBAR/ARED                                                00860000
      KE=0                                                              00870000
C-----------------------------                                          00880000
 50   ELAB=ELAB+1.                                                      00890000
      IF(ELAB.GT.150) GOTO 1999                                         00900003
      ECM=ELAB*AT/(AP+AT)                                               00910000
      SIG0=15.708*HW0/ECM                                               00920000
      LMX =-IABS(LMAX)                                                  00930000
C
C     I.Giese 2007
C     0.0D0
C
      CALL FUSION(ELAB,AP,ZP,AT,ZT,Q2,V0,R0,D,0.0D0,0.0D0,0.0D0,        00940004
     1            CRED,NORUTH,                                          00950005
     2            NOPROX,NOCURV,ION,SIGML,FLGRAZ,LMX ,SIGF,ITEST,IOPT)  00960004
      IF(ITEST.EQ.99) GOTO 1999                                         00970003
      IF(SIGF.GT.SIGLOW) GOTO 55                                        00980000
      TLIM=TLIM+1.                                                      00990000
      GOTO 50                                                           01000000
 55   KE=KE+1                                                           01010000
      DO 60 L=1,LMX                                                     01020000
      IDUM=KE + KDIM *(L-1)                                             01030000
      FL=L+L-1                                                          01040000
      TLE(IDUM)  =SIGML(L)/(SIG0*FL)                                    01050000
      IF(TLE(IDUM).LT.TLOW) TLE(IDUM)=0.                                01060000
      SIG(KE)=SIGF                                                      01070000
 60   CONTINUE                                                          01080000
      IF(KE.LT.KEMAX) GOTO 50                                           01090000
C------------------------------                                         01100000
      IF(IPRNT.NE.0) WRITE(6,315) K,LMAX,APAR,ZPAR,TLIM,V0,R0,D,Q2,R0T  01110000
 315  FORMAT(/' OWKB K,LMAX,A,Z,EMIN',2I3,3F8.1/  ' V0,R0,D,Q2=', F8.2, 01120000
     1 2F9.4,F8.1,4X,'R0T=',F9.4)                                       01130000
 999  RETURN                                                            01140000
1999  WRITE(6,2000) ELAB,ITEST                                          01150003
2000  FORMAT(' CALCULATION STOPS BECAUSE PROBLEM WITH WKB TRANSMISSION'/01160003
     1       ' ELAB,ITEST=',F10.1,I4)                                   01170003
      STOP                                                              01180003
      END                                                               01190000
