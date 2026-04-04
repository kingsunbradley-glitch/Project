      SUBROUTINE OVER2(K,KEMAX,LMAX,TLE,AMASS,ZEE,TLIM,SIGLOW,TLOW,     00010000
     1                 KDIM,IPRNT)                                      00020000
C
C  Berechnet Transmissionskoeffizienten fuer n,p,alphas fuer Abdampfung
C  nach dem optischen Modell, ruft TDL auf
C  Stand: 3.5,1994
C  
C     MOD 2 NOV 81 TO FEED IN POTS FROM MAIN                            00030000
C     MOD 20 JUN 90 SATCHLER POT FOR ALPHAS IS DEFAULT                  00040007
C     CALCULATES TRANSMISSION COEFFICIENTS FOR HIVAP                    00050000
C     FOR DEUTERON AND TRITON EVAPORATION THIS MUST BE COMPLETED BY INSE00060000
C     RELEVANT OPTICAL MODEL PARAMETERS                                 00070000
C     JE=1 CORRESPONDS TO E=TLIM                                        00080000
C     K=1 NEUTRONS,=2 PROTONS,=3 ALPHAS,=4 DEUTERONS,=5 TRITONS         00090000
C     KEMAX=MAXIMUM NUMBER OF ENERGIES FOR WHICH TRANSMISSION IS CALCUL.00100000
C     IF 'KEMAX' IS NEGATIVE,WILL INTERPRETE 'SIGLOW' AS ENERGY STEP(MEV00110000
C     ELSE STEP IS 1MEV AND SIGLOW IS CUTOFF ABSORB. CROSS SECTION(MB)  00120000
C     'TLIM' IS LOWEST ENERGY,'LMAX' IS MAX ORBITAL ANG.MOMENTUM        00130000
C                                                                       00140000
      IMPLICIT REAL*8 (A-H,O-Z)
C     UNUSED: ARRAY
      DIMENSION TLE (1) ,ARRAY(10)                                      00150000
      COMMON/TLJ/ T(2,31),V(15),V1(3)                                   00160000
      COMMON/POTS/ VRN,R0RN,ADIFRN,VIN,R0IN,ADIFIN,                     00170004
     1             VRP,R0RP,ADIFRP,VIP,R0IP,ADIFIP,RCLMBP,CBFACP,       00180004
     2             VRA,R0RA,ADIFRA,VIA,R0IA,ADIFIA,RCLMBA,CBFACA,QQ2    00190004
C      DOUBLE PRECISION T,V,V1,H4,H1,Z,W1,STPLTH                         00200000
C-------------------------------------------------------------------    00210000
C     MAYBE UNINITIALIZED: IDEFN,IDEFW,IDEFP,IDEFPW,CBFAC
C
C     I.Giese 2007
      SAVE
C
 1001 FORMAT(6E10.4)                                                    00220000
 1002 FORMAT(1I1)                                                       00230000
 1003 FORMAT(1H ,8F10.3  )                                              00240000
 1004 FORMAT(6F10.0)                                                    00250000
 1008 FORMAT(9F8.3)                                                     00260000
      DDE=1.                                                            00270000
      IF(KEMAX.GT.0.) GOTO 21                                           00280000
      KEMAX=IABS(KEMAX)                                                 00290000
      DDE=SIGLOW                                                        00300000
      SIGLOW=0.                                                         00310000
 21   IDUM=KEMAX*LMAX                                                   00320000
      DO 20 I=1,IDUM                                                    00330000
 20   TLE(I)=0.                                                         00340000
C--------------------------------------------------------------------   00350000
C                                      PARTICLE LOOP                    00360000
C                                                                       00370000
      K=IABS(K)                                                         00380000
      E=TLIM                                                            00390000
      JE1=TLIM+1.                                                       00400000
      E1=0.                                                             00410000
      E2=0.                                                             00420000
      E3=0.                                                             00430000
      TSUM1=0.                                                          00440000
      TSUM2=0.                                                          00450000
      TSUM3=0.                                                          00460000
      TSUM = 0.                                                         00470000
      XLM=LMAX                                                          00480000
      GOTO (1,2,3,320,340),K                                            00490000
C-----------------------------------   NEUTRON  --------------          00500000
C     DEFAULT IS    WILMORE-HODGSON NP55(64)673,ALSO PEREY NDTA10(72)53900510007
 1    XJP=0.50                                                          00520000
      XMP=1.00                                                          00530000
      XMT=AMASS-1.0                                                     00540000
      ZP=0.0                                                            00550000
      ZT=ZEE                                                            00560000
      VR   =VRN                                                         00570007
      IDEFN=0        ! USER POTENTIALS NOT ENERGY DEPENDENT             00580007
      RV   =R0RN                                                        00590000
      AV   =ADIFRN                                                      00600000
      W    =VIN                                                         00610000
      IDEFW=0        ! USER POTENTIALS NOT ENERGY DEPENDENT             00620007
      RW   =R0IN                                                        00630000
      AW   =ADIFIN                                                      00640000
      CBFAC=1.                                                          00650000
C                                                                       00660007
 200  IF(VR.LE.0.) IDEFN=1     ! ENERGY DEPENDENCE IN ENERGY LOOP       00670007
      IF(VR.LE.0.) VR=47.01                                             00680007
      IF(RV.LE.0.) RV=1.322  -7.6E-4*XMT  +4.E-6*XMT*XMT  -8.E-9*XMT**3.00690007
      IF(AV.LE.0.) AV=.660                                              00700000
      IF(W.LE.0.) IDEFW=1       ! ENERGY DEPENDENCE IN ENERGY LOOP      00710007
      IF(W.LE.0.) W=9.520                                               00720000
      IF(RW.LE.0.) RW=1.266-3.7E-4*XMT  +2.E-6*XMT*XMT  -4.E-9*XMT **3. 00730000
      IF(AW.LE.0.) AW=0.48                                              00740000
      S=1.00         ! SURFACE ABSORPTION FLAG                          00750007
C                                                                       00760007
 205  VSO=7.00       ! SPIN ORBIT                                       00770007
      RSO=RW                                                            00780000
      ASO=AW                                                            00790000
      RCLMB=.001                                                        00800000
      GO TO 4                                                           00810000
C  ---------------------------------   PROTON ---------------------     00820007
C      DEFAULT:    PROTON POT BECCHETTI-GREENLEES PR182(1969)1190       00830002
 2    XJP=0.50                                                          00840000
      XMP=1.00                                                          00850000
      XMT=AMASS-1.0                                                     00860000
      ZP=1.0                                                            00870000
      ZT=ZEE-1.0                                                        00880000
      VR   =VRP                                                         00890007
      IDEFP=0                                                           00900000
      RV   =R0RP                                                        00910000
      AV   =ADIFRP                                                      00920000
      W    =VIP                                                         00930000
      IDEFPW=0                                                          00940002
      RW   =R0IP                                                        00950000
      AW   =ADIFIP                                                      00960000
      RCLMB=RCLMBP                                                      00970000
      CBFAC =CBFACP                                                     00980000
 210  IF(VR.LE.0.) IDEFP=1       ! ENERGY DEPENDENCE IN ENERGY LOOP     00990007
      IF(VR.LE.0.)  VR=54.  +24.*(XMT -2.*ZT)/XMT  +0.4*ZT/XMT**0.33333 01000007
      IF(RV.LE.0.)  RV=1.17                                             01010002
      IF(AV.LE.0.)  AV=.750                                             01020002
      IF(W.LE.0.)   IDEFPW=1      ! ENERGY DEPENDENCE IN ENERGY LOOP    01030007
      IF(W.LE.0.)   W=11.8 + 12.*(XMT-2.*ZT)/XMT                        01040002
      IF(RW.LE.0.)  RW=1.32                                             01050002
      IF(AW.LE.0.)  AW=.51 + 0.7*(XMT-2.*ZT)/XMT                        01060002
      IF(RCLMB.LE.0.) RCLMB=1.16                                        01070002
      IF(CBFAC.LE.0.) CBFAC=1.                                          01080000
      S=1.00                        ! SURFACE ABSORPTION FLAG           01090007
C                                                                       01100007
 215  VSO=6.2                       ! SPIN ORBIT                        01110007
      RSO=1.01                                                          01120007
      ASO=0.75                                                          01130002
      GO TO 4                                                           01140000
C------------------------------------  ALPHA ---------------            01150000
C     DEFAULT                          SATCHLER,NP70(1965)177           01160007
 3    XJP=0.0                                                           01170000
      XMP=4.0                                                           01180000
      XMT=AMASS-4.0                                                     01190000
      ZP=2.00                                                           01200000
      ZT=ZEE-2.0                                                        01210000
      VR   =VRA                                                         01220007
      RV   =R0RA                                                        01230004
      AV   =ADIFRA                                                      01240000
      W    =VIA                                                         01250000
      RW   =R0IA                                                        01260000
      AW   =ADIFIA                                                      01270000
      RCLMB=RCLMBA                                                      01280000
      CBFAC=CBFACA                                                      01290000
C                                                                       01300007
      IF(VR.LE.0.)    VR=50.2                                           01310007
      IF(RV.LE.0.)    RV=1.2  +1.5/XMT**0.33333                         01320007
      IF(AV.LE.0.)    AV=.564                                           01330000
      IF(W.LE.0.)     W=12.30                                           01340000
      IF(RW.LE.0.)    RW=1.2  +1.5/XMT**0.33333                         01350007
      IF(AW.LE.0.)    AW=AV                                             01360000
      IF(RCLMB.LE.0.) RCLMB=1.30                                        01370000
      IF(CBFAC.LE.0.) CBFAC=1.                                          01380000
      S=0.0           ! VOLUME ABSORPTION FLAG                          01390007
      RSO=1.                                                            01400007
      VSO=.001                                                          01410000
      ASO=1.0                                                           01420000
      GOTO 4                                                            01430000
C------------------------------------  DEUTERON ----------------        01440000
 320  CONTINUE                                                          01450000
      GOTO 4                                                            01460000
C------------------------------------  TRITON ---------------           01470000
 340  CONTINUE                                                          01480000
C-----------------------------------------------------------------------01490000
 4    IF(IPRNT.NE.1) GOTO 350                                           01500000
      WRITE (6,1005)                                                    01510000
 1005 FORMAT( /' PROJ. SPIN  A(PROJ.) A(TARGET)  Z(PROJ.) Z(TARGET)')   01520005
      WRITE(6,1003)XJP,XMP,XMT,ZP,ZT                                    01530000
      WRITE (6,1006)                                                    01540000
 1006 FORMAT( /  '    V(REAL)   R(REAL)   A(REAL)  W(IMAG.)  R(IMAG.)  A01550005
     1(IMAG.) ')                                                        01560000
      WRITE(6,1003)VR,RV,AV,W,RW,AW                                     01570007
      WRITE (6,1007)                                                    01580000
 1007 FORMAT (/  '      V(SO)     R(SO)     A(SO)   R(COUL)        S   L01590005
     1(LIMIT) ')                                                        01600000
      WRITE(6,1003)VSO,RSO,ASO,RCLMB,S,XLM                              01610000
 350  CONTINUE                                                          01620000
      XF=XMT**.3333                                                     01630000
      XM=XMT/(XMT + XMP)                                                01640000
      XJT=0.0                                                           01650000
      IC=0                                                              01660000
C-----------------------------------------------------------------------01670000
C                                      ENERGY LOOP                      01680000
C                                                                       01690000
      KE=0                                                              01700000
      DO 130 JE=JE1,80                                                  01710000
 370  V(1)=RV*XF                                                        01720000
      V(2)=AV                                                           01730000
      V(3)=RW*XF                                                        01740000
      V(4)=AW                                                           01750000
      V(5)=S                                                            01760000
      V(6)=VSO                                                          01770000
      V(9)=RCLMB*XF                                                     01780000
      V(10)=RSO*XF                                                      01790000
      V(11)=ASO                                                         01800000
      V1(1)=XJP                                                         01810000
      V1(2)=+0.0                                                        01820000
      W1=0.04783258*XMP*XM                                              01830000
      STPLTH=0.1                                                        01840000
      V1(3)=XLM                                                         01850000
      V(7)=VR                                                           01860007
      IF(K.EQ.1 .AND. IDEFN.EQ.1) V(7)=VR-0.267*E   ! ENERGY DEP        01870007
      IF(K.EQ.2 .AND. IDEFP.EQ.1) V(7)=VR-0.32 *E                       01880007
      V(8)=W                                                            01890000
      IF(K.EQ.1 .AND. IDEFW.EQ.1)V(8)=W-0.053*E                         01900000
      IF(K.EQ.2 .AND. IDEFPW.EQ.1)V(8)=W-0.25*E                         01910002
      H4=0.04783*XMP*(XM**2)*E                                          01920000
      H1=DSQRT(H4)                                                      01930000
      Z=CBFAC*(0.03478*ZP*ZT*XMP*XM)/H1                                 01940000
      C6=31.42/H4                                                       01950000
      DO 5 J=1,2                                                        01960000
      DO 5 I=1,31                                                       01970000
 5    T(J,I)=0.0                                                        01980000
C--------------------------------------------------------------------   01990000
      CALL TLD(H4,N1,J1,H1,Z,W1,STPLTH)                                 02000000
C--------------------------------------------------------------------   02010000
      TSUM=0.0                                                          02020000
      AVSUM=0.                                                          02030000
      INTRPO=2.*XJP                                                     02040000
C                                                                       02050000
      DO 7 I=1,LMAX                                                     02060000
      AL=I+I-1                                                          02070000
      IF(INTRPO.EQ.0) GOTO 6                                            02080000
      TK=I-1                                                            02090000
      FL=1./(2.*TK+1)                                                   02100000
      T(1,I)=FL*(T(1,I)*TK + T(2,I)*(TK+1.))                            02110000
 6    TSUM=TSUM + T(1,I)*AL*C6                                          02120000
 7    CONTINUE                                                          02130000
C                                                                       02140000
      IF(TSUM.GE.SIGLOW  ) GOTO 70                                      02150000
      TLIM=TLIM+1.                                                      02160000
      GOTO 10                                                           02170000
 70   KE=KE+1                                                           02180000
      DO 8 I=1,LMAX                                                     02190000
      IDUM=KE + KDIM *(I-1)                                             02200000
      TLE(IDUM)=T(1,I)                                                  02210000
      IF(TLE(IDUM).LT.TLOW) TLE(IDUM)=0.                                02220000
 8    CONTINUE                                                          02230000
      IF(KE.GE.KEMAX) GOTO 14                                           02240000
 10   CONTINUE                                                          02250000
      E=E+DDE                                                           02260000
 130  CONTINUE                                                          02270000
C                                                                       02280000
C                                      STOP ENERGY LOOP                 02290000
C                                                                       02300000
 14   CONTINUE                                                          02310000
C                                                                       02320000
C                                      STOP PARTICLE LOOP               02330000
C                                                                       02340000
 15   RETURN                                                            02350000
      END                                                               02360000
