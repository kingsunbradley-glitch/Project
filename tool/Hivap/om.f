      SUBROUTINE OM (K,ELAB,LMAX,ATARG,ZTARG,SIGTOT,SIGML,IPOT,IPRNT)   00020007
C 
C  Berechnet Absorptionsquerschnitte fuer n,p,a nach dem optischen 
C  Modell
C  Stand: 3.5.1994
C
C     CALCULATES ABSORPTION CROSS SECTIONS FOR NEUTRONS,PROTONS,ALPHAS  00040001
C     FOR DEUTERON AND TRITON  THIS MUST BE COMPLETED BY INSERTING      00050001
C     RELEVANT OPTICAL MODEL PARAMETERS                                 00060000
C     K=1 NEUTRONS,=2 PROTONS,=3 ALPHAS,=4 DEUTERONS,=5 TRITONS         00070000
C     ELAB    =PROJECTILES ENERGY IN MEV                                00080007
C     LMAX    MAXIMUM+1 ORBITAL ANGULAR MOMENTUM                        00090007
C     ATARG,ZTARG  MASS AND Z OF TARGET NUCLEUS                         00100001
C     SIGTOT  TOTAL ABSORPTION CROSS SECTION IN MILLIBARN(OUTPUT)       00110001
C     SIGML   PARTIAL CROSS SECTIONS(AVERAGED FOR SPIN ORBIT) (OUTPUT)  00120007
C     IPOT    =0  TAKE DEFAULT POTENTIAL                                00130007
C             =1  READ IN NEW POTENTIAL                                 00140007
C             =2  KEEP OLD POTENTIAL                                    00150007
C     IPRNT   PRINT FLAG                                                00160007
C     T       TRANSMISSION COEFFICIENTS (OUTPUT)                        00170007
C---------------------------------------------------------------------- 00180007
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIGML (1) ,ARRAY(10)                                    00190001
      COMMON/TLJ/ T(2,31),V(15),V1(3)                                   00200000
      DOUBLE PRECISION T,V,V1,H4,H1,Z,W1,STPLTH                         00210000
C-------------------------------------------------------------------    00220000
C     MAYBE UNINITIALIZED: KOLD
C
C     I.Giese 2007
      SAVE
C
 1001 FORMAT(6E10.4)                                                    00230000
 1002 FORMAT(1I1)                                                       00240000
 1003 FORMAT(1H ,8F10.3  )                                              00250000
 1009 FORMAT(1H ,8F10.3//)                                              00260004
 1004 FORMAT(6F10.0)                                                    00270000
 1008 FORMAT(9F8.3)                                                     00280000
      IF(LMAX.GT.31) LMAX=31                                            00290001
      DO 20 I=1,LMAX                                                    00300001
 20   SIGML(I)=0.                                                       00310001
C--------------------------------------------------------------------   00320000
C ------------------------------------ PARTICLE SWITCH                  00330007
C                                                                       00340000
      K=IABS(K)                                                         00350000
      ISAME=0                                                           00360007
      IF(K.EQ.KOLD) ISAME=1                                             00370007
      KOLD=K                                                            00380007
      IPOT1=IPOT                                                        00390007
      IF(ISAME.EQ.0 .AND.IPOT.EQ.2) IPOT1=0                             00400007
      E=ELAB                                                            00410001
      TSUM = 0.                                                         00420000
      XLM=LMAX                                                          00430000
      GOTO (1,2,3,320,340),K                                            00440000
C------------------------------------- NEUTRON                          00450007
 1    XJP=0.50                                                          00460000
      XMP=1.00                                                          00470000
      XMT=ATARG                                                         00480001
      ZP=0.0                                                            00490000
      ZT=ZTARG                                                          00500001
      IF(IPOT1.EQ.0) GOTO 200                                           00510007
      IF(IPOT1.EQ.2) GOTO 205                                           00520007
C     READ(5,1008) P,RV,AV,W,RW,AW                                      00530000
      IF(IPRNT.NE.13) CALL MYRD(ARRAY,7,6,5)                            00540013
      IF(IPRNT.EQ.13) CALL MYRD(ARRAY,7,13,5)                           00550013
      P    =ARRAY(1)                                                    00560000
      RV   =ARRAY(2)                                                    00570000
      AV   =ARRAY(3)                                                    00580000
      W    =ARRAY(4)                                                    00590000
      RW   =ARRAY(5)                                                    00600000
      AW   =ARRAY(6)                                                    00610000
C     RCLMB=ARRAY(7)                                                    00620000
      GOTO 205                                                          00630000
 200  P=47.01 - 0.267*E                                                 00640011
      RV=1.322-7.6E-4*XMT  +4.E-6*XMT*XMT      -8.E-9* XMT **3.         00650000
      AV=.660                                                           00660000
      W=9.520 -0.053*E                                                  00670011
      RW=1.266-3.7E-4*XMT  +2.E-6*XMT  *XMT  -4.E-9*XMT **3.            00680000
      AW=0.48                                                           00690000
 205  VSO=7.00                                                          00700000
      RSO=RW                                                            00710000
      ASO=AW                                                            00720000
      RCLMB=.001                                                        00730000
      S=1.00                 ! SURFACE IMAG                             00740014
      GO TO 4                                                           00750000
C ------------------------------------ PROTON                           00760007
 2    XJP=0.50                                                          00770000
      XMP=1.00                                                          00780000
      XMT=ATARG                                                         00790001
      ZP=1.0                                                            00800000
      ZT=ZTARG                                                          00810001
      IF(IPOT1.EQ.0) GO TO 210                                          00820007
      IF(IPOT1.EQ.2) GO TO 215                                          00830007
C     READ(5,1008) P,RV,AV,W,RW,AW,RCLMB                                00840000
      IF(IPRNT.NE.13) CALL MYRD(ARRAY,7,6,5)                            00850013
      IF(IPRNT.EQ.13) CALL MYRD(ARRAY,7,13,5)                           00860013
      P    =ARRAY(1)                                                    00870000
      RV   =ARRAY(2)                                                    00880000
      AV   =ARRAY(3)                                                    00890000
      W    =ARRAY(4)                                                    00900000
      RW   =ARRAY(5)                                                    00910000
      AW   =ARRAY(6)                                                    00920000
      RCLMB=ARRAY(7)                                                    00930000
      GOTO 215                                                          00940000
C     BECCHETTI-GREENLEES PR182(1969)1190                               00950012
 210  P=54.        +24.*(XMT  -2.*ZT)/XMT  +0.4*ZT/XMT  **0.33333       00960012
     1             -0.32*E                                              00970012
      RV=1.17                                                           00980012
      AV=.750                                                           00990012
      W=11.8 +12.*(XMT - 2.*ZT)/XMT-0.25*E                              01000012
      RW=1.32                                                           01010012
      AW=.51 + 0.7*(XMT - 2.*ZT)/XMT                                    01020012
      RCLMB=1.16                                                        01030012
C     SPIN ORBIT PEREY 1963                                             01040011
 215  RSO=1.01                                                          01050012
      VSO=6.2                                                           01060012
      ASO=.75                                                           01070012
      S=1.00             ! SURFACE IMAG                                 01080014
      GO TO 4                                                           01090000
C--------------------------------------ALPHA                            01100007
 3    XJP=0.0                                                           01110000
      XMP=4.0                                                           01120000
      XMT=ATARG                                                         01130001
      ZP=2.00                                                           01140000
      ZT=ZTARG                                                          01150001
      IF(IPOT1.EQ.0) GO TO 220                                          01160007
      IF(IPOT1.EQ.2) GO TO 225                                          01170007
C     READ(5,1008) P,RV,AV,W,RW,AW,RCLMB                                01180000
      IF(IPRNT.NE.13) CALL MYRD(ARRAY,7,6,5)                            01190013
      IF(IPRNT.EQ.13) CALL MYRD(ARRAY,7,13,5)                           01200013
      P    =ARRAY(1)                                                    01210000
      RV   =ARRAY(2) +1.5/XMT**0.33333                                  01220000
      AV   =ARRAY(3)                                                    01230000
      W    =ARRAY(4)                                                    01240000
      RW   =ARRAY(5)                                                    01250000
      AW   =ARRAY(6)                                                    01260000
      RCLMB=ARRAY(7)                                                    01270000
      GOTO 225                                                          01280000
 220  P=50.2                                                            01290000
C     RV=1.2+1.5/XMT**0.33333          SATCHLER                         01300000
      RV=1.275+1.5/XMT**0.33333                                         01310000
      AV=.564                                                           01320000
      W=12.30                                                           01330000
C     RW=RV                            SATCHLER                         01340000
      RW=1.25                                                           01350000
      AW=.564                                                           01360012
      RCLMB=1.30                                                        01370000
 225  RSO=1.                                                            01380000
      VSO=.001                                                          01390000
      ASO=1.0                                                           01400000
      S=0.0              ! VOLUME IMAG                                  01410014
      GOTO 4                                                            01420000
C                                      DEUTERON                         01430000
 320  CONTINUE                                                          01440000
      GOTO 4                                                            01450000
C                                      TRITON                           01460000
 340  CONTINUE                                                          01470000
C------------------------------------------------------------------     01480003
 4    IF(IPRNT.EQ.0) GOTO 350                                           01490003
      WRITE (6,1005)                                                    01500000
 1005 FORMAT(/' PROJ. SPIN  A(PROJ.) A(TARGET)  Z(PROJ.) Z(TARGET)')    01510004
      WRITE(6,1003)XJP,XMP,XMT,ZP,ZT                                    01520000
      WRITE (6,1006)                                                    01530000
 1006 FORMAT( '    V(REAL)   R(REAL)   A(REAL)  W(IMAG.)  R(IMAG.)  A   01540003
     1(IMAG.) ')                                                        01550000
      WRITE(6,1003)P,RV,AV,W,RW,AW                                      01560000
      WRITE (6,1007)                                                    01570000
 1007 FORMAT ('      V(SO)     R(SO)     A(SO)   R(COUL)        S   L   01580003
     1(LIMIT) ')                                                        01590000
      WRITE(6,1009)VSO,RSO,ASO,RCLMB,S,XLM                              01600004
C---------------------------------------------------------------------  01610003
      IF(IPRNT.NE.13) GOTO 350                                          01620003
      WRITE (13,1005)                                                   01630003
      WRITE(13,1003)XJP,XMP,XMT,ZP,ZT                                   01640003
      WRITE (13,1006)                                                   01650003
      WRITE(13,1003)P,RV,AV,W,RW,AW                                     01660003
      WRITE(13,1007)                                                    01670003
      WRITE(13,1009)VSO,RSO,ASO,RCLMB,S,XLM                             01680004
C---------------------------------------------------------------------  01690003
 350  CONTINUE                                                          01700000
      XF=XMT**.3333                                                     01710000
      XM=XMT/(XMT + XMP)                                                01720000
      XJT=0.0                                                           01730000
      IC=0                                                              01740000
C                                                                       01750000
C--------------------------------------------------------------------   01760001
C                                                                       01770000
      V(1)=RV*XF                                                        01780002
      V(2)=AV                                                           01790000
      V(3)=RW*XF                                                        01800000
      V(4)=AW                                                           01810000
      V(5)=S       ! IS ONE FOR SURFACE IMAG                            01820014
      V(6)=VSO                                                          01830000
      V(9)=RCLMB*XF                                                     01840000
      V(10)=RSO*XF                                                      01850000
      V(11)=ASO                                                         01860000
      V1(1)=XJP                                                         01870000
      V1(2)=+0.0                                                        01880000
      W1=0.04783258*XMP*XM                                              01890000
      STPLTH=0.1                                                        01900000
      V1(3)=XLM                                                         01910000
      V(7)=P                                                            01920000
      V(8)=W                                                            01930000
      H4=0.04783*XMP*(XM**2)*E                                          01940000
      H1=DSQRT(H4)                                                      01950000
      Z=(0.03478*ZP*ZT*XMP*XM)/H1                                       01960000
      C6=31.42/H4                                                       01970000
      DO 5 J=1,2                                                        01980000
      DO 5 I=1,31                                                       01990000
 5    T(J,I)=0.0                                                        02000000
C--------------------------------------------------------------------   02010000
      CALL TLD(H4,N1,J1,H1,Z,W1,STPLTH)                                 02020000
C--------------------------------------------------------------------   02030000
      TSUM=0.0                                                          02040000
      AVSUM=0.                                                          02050000
      INTRPO=2.*XJP                                                     02060000
C                                                                       02070000
      DO 7 I=1,LMAX                                                     02080000
      AL=I+I-1                                                          02090000
      IF(INTRPO.EQ.0) GOTO 6                                            02100000
      TK=I-1                                                            02110000
      FL=1./(2.*TK+1)                                                   02120000
      T(1,I)=FL*(T(1,I)*TK + T(2,I)*(TK+1.))                            02130000
 6    SIGML(I)=T(1,I)*AL*C6                                             02140001
      TSUM=TSUM+SIGML(I)                                                02150001
 7    CONTINUE                                                          02160000
      SIGTOT=TSUM                                                       02170001
C                                                                       02180000
C                                                                       02190000
C------------------------------------------------------------------     02200001
C                                                                       02210000
      RETURN                                                            02220001
      END                                                               02230000
