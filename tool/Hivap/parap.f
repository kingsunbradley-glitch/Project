      SUBROUTINE PARAP(ELAB,AP,AT,ZP,ZT,FINERT,RATIOS,IFUS,LIMBAR,JCAL, 00020006
     1                NOPRNT,  MX,SIGML,CRSROT,FLGRAZ)                  00030003
C
C  Berechnet Fusionsquerschnitte mit Parabelmethode (Hill-Wheeler)
C  Stand: 3.5.1994
C
C     QUICK CALCULATION OF NUCLEAR CROSS SECTIONS VIA OPTICAL MODEL     00040000
C     WITH PARABOLIC BARRIER - SEE THOMAS, PHYS REV 116,703 (1959) -    00050000
C                                                                       00060003
C     STANDARD V0,R0,D = 67.,1.17,0.574 (BLANN-PLASIL)                  00070003
C     IFUS=5   V0,R0,D = 70.,1.25,0.44 (SIKKELAND)                      00080003
C     IFUS=7     READS V0,R0,D FROM CARDS                               00090003
C     LIMBAR=1        LMAX LIMITED BY      FISSION BARRIER>0.1MEV       00100003
C     OUTPUT MX IS EROT+0.5 FOR MAX AM                                  00110003
C     FLGRAZ OUTPUT:GRAZING AM,                                         00120003
C     NOPRNT=1 REDUCED PRINT OUT                                        00130003
C     JCAL=0 STANDARD                                                   00140003
C     JCAL=2 SIMPLIFIED VERSION :NO AM BUT EFFECTIVE EXCIT (EXCIT-EROT) 00150003
C            WHERE EROT=YRAST FROM MODIFIED LDR                         00160003
C     JCAL=3 SAME AS '2' BUT EROT=SPHERIC NUCLEUS(RIGID)                00170003
C     FINERT,RATIOS(1-6) PARAMETERS FOR MODIFIED LDR                    00180006
C     SIGML(L) PARTIAL CROSS SECTIONS                                   00190003
C     CRSROT(EROT+0.5) SUM OF SIGML CORESPONDING TO ROT.ENERGY AROUND   00200003
C                      EROT(1MEV INTERVALS),NEEDED FOR JCAL=2,3         00210003
C-----------------------------------------------------------------      00220005
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIGML(1),CRSROT(1),RATIOS(6),ARRAY(4)                   00230010
      COMMON/FISRT/ROT0,A2MS,RKAPPA,X0,Y0                               00240007
C                                                                       00250000
C     MAYBE UNINITIALIZED: TL0,TL1
C
C     I.Giese 2007
      SAVE
C
      LMAX=200                                                          00260006
      DO1 I=1,LMAX                                                      00270006
      CRSROT(I)=0.                                                      00280000
 1    SIGML(I)=0.                                                       00290000
 2    CONTINUE                                                          00300001
C------------------------------------- POTENTIAL PARAMETERS             00310005
      F=67.                                                             00320000
      IF(AP.LE.4.)F=1100.                                               00330000
      R00=1.17                                                          00340000
      D=0.574                                                           00350000
      IF(IFUS.NE.5) GOTO 3                                              00360003
      F=70.                                                             00370000
      R00=1.25                                                          00380000
      D=0.44                                                            00390000
C3    IF(IFUS.EQ.7) READ(5,1015) F,R00,D                                00400010
 3    IF(IFUS.NE.7) GOTO 4                                              00410010
      CALL MYRD(ARRAY,3,23,5)                                           00420011
      F=ARRAY(1)                                                        00430010
      R00=ARRAY(2)                                                      00440010
      D=ARRAY(3)                                                        00450010
 4    WRITE(6,1016) F,R00,D,ELAB                                        00460010
C---------------------------------------------------------------        00470005
      ACN=AP+AT                                                         00480006
      E  =ELAB*AT/ACN                                                   00490006
      AR=AT**0.333333+AP**0.333333                                      00500005
      IF(AP.LE.4.) AR=AT**0.333333                                      00510005
      COUL= 1.4393 *ZT*ZP                                               00520006
      RA=R00*AR                                                         00530000
      U=AT*AP/(AT+AP)                                                   00540000
      H2=41.814                                                         00550000
      CROTL=H2/(1.16*ACN**1.666666)                                     00560000
      CR2=D*COUL/F                                                      00570000
      CR3P=H2*D/(F*U)                                                   00580000
      PWAVE=SQRT(2.*931.16*U*E)                                         00590005
      WAVEL=19.732/PWAVE                                                00600005
      ARWAVE=3141.59*(WAVEL**2)                                         00610005
C-------------------------------------------------------------          00620005
 8    IF(NOPRNT.EQ.1) GOTO 80                                           00630001
      WRITE(6,1002)                                                     00640001
 1002 FORMAT(3H0  )                                                     00650000
      WRITE(6,1003)                                                     00660000
      WRITE(6,1004)                                                     00670000
 1003 FORMAT(61H0 PROJECTILE  TARGET   PROJECTILE  TARGET    ENERGY    E00680000
     1NERGY  )                                                          00690000
 1004 FORMAT(61H    MASS       MASS      CHARGE    CHARGE     LAB       00700000
     1 CM    )                                                          00710000
      WRITE(6,1005)                                                     00720000
 1005 FORMAT(3H   )                                                     00730000
      WRITE(6,1006) AP,AT,ZP,ZT,ELAB,E                                  00740003
 1006 FORMAT(4X,F7.3,5(3X,F7.3))                                        00750000
      WRITE(6,1002)                                                     00760000
      WRITE(6,1007)                                                     00770000
 1007 FORMAT(124H ANGULAR TRANSMISSION BARRIER BAR-RADIUS  BARRIER    BA00780000
     1RRIER   RELATIVE CURVATURE  PROBABLE SIGMA(L) ROTOR-E EFFROT EFFSI00790000
     2GMA    )                                                          00800000
      WRITE(6,1008)                                                     00810000
 1008 FORMAT(124H MOMENTA COEFFICIENTS  RADII   CONSTANT    (MEV)    CUR00820000
     1VATURE  SIGMA(L)   SIGN     L-ENERGY  (MB)     (MEV)   (MEV)  (MB)00830000
     2       )                                                          00840000
      WRITE(6,1005)                                                     00850000
      WRITE(6,1009)                                                     00860000
 1009 FORMAT(122H     L         TL        RB        R0       VLB        00870000
     1WL       RSIGL    SD2VLB     SIGXL2   SIGL     EROTL  EROTLI SUMDX00880000
     2L    )                                                            00890000
      WRITE(6,1005)                                                     00900001
C---------------------------------------------------------------------  00910001
C                                                                       00920001
 80   SUMAL2=0.0                                                        00930005
      SUMTL=0.0                                                         00940000
      SUMDXL=0.0                                                        00950000
      EROTLI=0.0                                                        00960000
      IREP=0                                                            00970000
      TL2=0.                                                            00980000
C                                      PARTIAL WAVE LOOP                00990000
      DO 29 K=1,LMAX                                                    01000000
      L=K-1                                                             01010000
      R=RA                                                              01020000
      AL=FLOAT(L)                                                       01030000
      CR3=CR3P*AL*(AL+1)                                                01040000
 9    RB=RA-D*DLOG(CR2/R**2+CR3/R**3)                                   01050000
      DELTR=ABS(RB-R)                                                   01060000
      IF(DELTR-0.001) 11,10,10                                          01070000
 10   R=RB                                                              01080000
      GO TO 9                                                           01090000
 11   CONTINUE                                                          01100000
      R0=RB/AR                                                          01110000
      VLB=COUL/RB+(H2/2./U)*(AL*(AL+1)/RB**2)-F*EXP(-(RB-RA)/D)         01120000
      D2VLB=2.*COUL/RB**3+(H2/U)*3.*AL*(AL+1)/RB**4                     01130000
     1 -(F/D**2)*EXP(-(RB-RA)/D)                                        01140000
      IF(D2VLB) 12,13,14                                                01150000
 12   SD2VLB=-1.0                                                       01160000
      GO TO 15                                                          01170000
 13   SD2VLB=0.0                                                        01180000
      GO TO 15                                                          01190000
 14   SD2VLB=1.0                                                        01200000
 15   CONTINUE                                                          01210000
      WL2=ABS((H2/U)*D2VLB)                                             01220000
      WL=SQRT(WL2)                                                      01230000
      DUM=2.*3.1459*(VLB-E)/WL                                          01240000
      IF(DUM.LT.-20.) GOTO 152                                          01250000
      TL=1./(1.+EXP(DUM))                                               01260000
      GOTO 155                                                          01270000
 152  TL=1.                                                             01280000
 155  CONTINUE                                                          01290000
      RSIGL=TL*(2.*AL+1.)                                               01300000
      SIGL=ARWAVE*RSIGL                                                 01310000
      SIGML(K)=SIGL                                                     01320000
      IF(TL.LT.0.5) IREP=IREP+1                                         01330000
      IF(IREP.GE.2) GOTO 160                                            01340000
      FLGRAZ=L-1                                                        01350000
      TL1=TL2                                                           01360000
      TL2=TL                                                            01370000
      IF(L.EQ.0) TL0=TL                                                 01380000
 160  CONTINUE                                                          01390000
C-----------------------------------------------                        01400004
      A=ACN                                                             01410003
      Z=ZP+ZT                                                           01420000
      AN=A-Z                                                            01430000
      ONLYGS=0.                                                         01440007
      CALL FISROT(A,Z,AN,AL,DELR,DELSP,ERO,ONLYGS)                      01450007
      BAR=DELSP-DELR                                                    01460000
C     WRITE(6,161)BAR,DELSP,DELR
C 161 FORMAT('Bar=',F8.3,' DELSP=',F8.3,' DELR=',F8.3)
      IF(LIMBAR.EQ.0) GOTO 170                                          01470003
      IF(BAR.GT.0.1) GOTO 170                                           01480000
      FLGRAZ=AL                                                         01490000
      GOTO 30                                                           01500000
 170  CONTINUE                                                          01510000
C----------------------------------------------EROTL                    01520005
      EROTL=ERO                                                         01530003
      IF(JCAL.EQ.3) GOTO 19                                             01540003
      RIGITY=FINERT                                                     01550004
C     IF(AL.GE.RATIOS(6))GOTO 18                                        01560008
C     CALL IP2(RATIOS(2),RATIOS(4),RATIOS(6),AL,RATIOS(1),RATIOS(3),    01570008
C    1         RATIOS(5),RIGITY)                                        01580008
 18   EROTL=DELR*(AL+1.)/RIGITY                                         01590004
      IF(AL.GT.0.)EROTL=EROTL/AL                                        01600004
C----------------------------------------------                         01610005
 19   SIGXL2=AL*(AL+1.)*RSIGL                                           01620000
      IF(EROTL-0.5-EROTLI) 20,20,22                                     01630000
 20   SUMDXL=SUMDXL+SIGL                                                01640000
      IF(TL-0.0001    ) 22,21,21                                        01650000
 21   CONTINUE                                                          01660000
      GO TO 25                                                          01670000
 22   IF(L) 24,24,23                                                    01680000
 23   CONTINUE                                                          01690000
C23   WRITE(6,1010) EROTLI,SUMDXL                                       01700000
      LX=EROTL+.5                                                       01710000
      CRSROT(LX)=SUMDXL                                                 01720000
      MX=LX                                                             01730003
 1010 FORMAT(3H+  ,107X,F5.1,1X,F7.2)                                   01740000
 24   CONTINUE                                                          01750000
      EROTLI=EROTLI+1.0                                                 01760000
      SUMDXL=SIGL                                                       01770000
 25   CONTINUE                                                          01780000
      IF(TL.GT.0.99) GOTO 28                                            01790000
C----------------------------------------                               01800005
      IF(NOPRNT.EQ.1) GOTO 28                                           01810001
      WRITE(6,1011) L,TL,RB,R0,VLB,WL,RSIGL,SD2VLB,SIGXL2,SIGL,EROTL    01820000
 1011 FORMAT(3X,I3,4X,F10.8,2X,F6.3,4X,F6.3,4X,F6.1,4X,F6.3,4X,F6.2,4X,F01830000
     16.3,1X,F10.2,2X,F7.3,3X,F7.3)                                     01840000
      IF(L) 26,26,28                                                    01850000
 26   IF(TL-.0001)27,28,28                                              01860000
 27   WRITE(6,1010) EROTLI,SUMDXL                                       01870000
C---------------------------------------                                01880005
 28   CONTINUE                                                          01890000
      SUMTL=SUMTL+RSIGL                                                 01900000
      SUMAL2=SUMAL2+SIGXL2                                              01910000
      IF(TL-.0001)30,29,29                                              01920000
 29   CONTINUE                                                          01930000
C                                                                       01940004
C                                      END PARTIAL WAVE LOOP            01950000
C---------------------------------------------------------------------  01960004
 30   SIGMA=ARWAVE*SUMTL                                                01970005
      IF(LIMBAR.EQ.1) GOTO 40                                           01980003
      IF(TL0.LT.0.5) GOTO 32                                            01990000
      FLGRAZ=FLGRAZ+(TL1-0.5)/(TL1-TL2)                                 02000000
 32   IF(FLGRAZ.LT.0.)FLGRAZ=0.                                         02010000
      WRITE(6,1012) SIGMA,FLGRAZ                                        02020000
 1012 FORMAT(' PARAP X-SECTION(MBARN) =',  F8.1,4X,'LGRAZE',F8.1 /)     02030002
      AVEL2=SUMAL2/SUMTL                                                02040000
      AVROE=AVEL2*CROTL                                                 02050000
      RMSL2=SQRT(AVEL2)                                                 02060000
      IF(NOPRNT.EQ.0)WRITE(6,1013) RMSL2,AVEL2,AVROE                    02070001
      GOTO 50                                                           02080000
 40   WRITE(6,1017) SIGMA,FLGRAZ                                        02090000
 50   CONTINUE                                                          02100000
 1013 FORMAT(23H0ROOT MEAN SQUARE L =  ,F6.2,10X,16H RMSL SQUARED = ,F1002110000
     1.3,10X,33H AVERAGE ROTATION ENERGY (MEV) = ,F6.3)                 02120000
 1014 FORMAT(6E10.3)                                                    02130000
 1015 FORMAT(10F8.3)                                                    02140000
 1016 FORMAT(/' PARAP',6X,'V0=',F8.3,6X,'R0=',F8.3,6X,'D=',F8.3,6X,     02150001
     1 'ELAB=',F8.1)                                                    02160001
 1017 FORMAT(' ZERO FISSION BARRIER LIMITED  SIGFUS,LCRIT',F10.3,F10.0) 02170000
      RETURN                                                            02180000
      END                                                               02190000
