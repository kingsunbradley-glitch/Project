      SUBROUTINE OUT2(POP,M,N,M0,N0,ISUPP,DIG1,MXDD,INV,LIN)            00030014
C     
C  Druck-Routine fuer zweidimensionales Array, modifiziert am 18.11.88
C  Stand: 3.5.1994
C  
C     MODIFIED TO ALLOW CALL BY VALUE OF DIG AND MXD 18 NOV 86  
C     IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),AM(15)                                           00040011
      CHARACTER*1 ARR(128)                                              00050011
      COMMON/S/BLANK,NU,SYMB                                            00060000
C     MAYBE UNINITIALIZED: J
      CHARACTER*1 BLANK,NU(15),SYMB(7)                                  00070011
      CHARACTER*4 DASH(2)                                               00080011
      DATA DASH/' M> ',' N> '/                                          00090004
C---------------------------------------------------------------------  00100000
C
C     I.Giese 2007
      SAVE
C
      AMX=0.                                                            00110000
      MN=M*N                                                            00120000
      DO 10 I=1,MN                                                      00130000
 10   AMX=AMAX1(AMX,POP(I))                                             00140000
      IF(AMX.LE.0.) RETURN                                              00150000
C                                                                       00160000
      DIG=DIG1                                                          00170012
      IF(DIG.LE.1.)DIG=2.                                               00180002
      MXD=MXDD                                                          00190014
      IF(MXD.LE.0)MXD=5                                                 00200000
      MXD=MIN0(15,MXD)                                                  00210001
C                                                                       00220000
      AM(1)=AMX                                                         00230002
      IF(LIN.EQ.1) GOTO 30                                              00240002
      DO 20 I=2,15                                                      00250002
 20   AM(I)=AM(I-1)/DIG                                                 00260002
      GOTO 35                                                           00270002
 30   DIG=MXD                                                           00280002
      DUM=AMX/DIG                                                       00290002
      DO 32 I=2,15                                                      00300002
 32   AM(I)=AM(I-1)-DUM                                                 00310002
C                                                                       00320002
 35   IH=MAX0(M,N)                                                      00330008
      IV=MIN0(M,N)                                                      00340008
      IDH=1                                                             00350000
      IHH=IH                                                            00360001
      IF(IH.LT.128) GOTO 50                                             00370000
 40   IDH=IDH+1                                                         00380000
      IHH=IH /IDH                                                       00390009
      IF(IHH.GT.128)GOTO 40                                             00400000
 50   IDV=1                                                             00410000
      IVV=IV                                                            00420001
      IF(IV.LT.60) GOTO 60                                              00430000
 55   IDV=IDV+1                                                         00440000
      IVV=IV /IDV                                                       00450009
      IF(IVV.GT.60) GOTO 55                                             00460007
C                                                                       00470000
 60   IF(IH.EQ.IV) GOTO 65                                              00480001
      IF(IH.EQ.N) GOTO 70                                               00490001
C---------------------------------------- M IS HORIZONTAL               00500008
 65   IPARH=1                                                           00510001
      IF(INV.EQ.1 .OR. INV.EQ.3) IPARH=2                                00520001
      IPARV=3                                                           00530001
      IF(INV.GT.2) IPARV=4                                              00540001
      IDASH=1                                                           00550004
      J00=N0-1                                                          00560003
      GOTO 80                                                           00570001
C---------------------------------------- N IS HORIZONTAL               00580008
 70   IPARH=3                                                           00590001
      IF(INV.GE.2) IPARH=4                                              00600001
      IPARV=2                                                           00610008
      IF(INV.EQ.1 .OR. INV.EQ.3) IPARV=1                                00620008
      IDASH=2                                                           00630004
      J00=M0-1                                                          00640003
C-----------------------------------------------------------------      00650008
 80   DO 90 I=1,128                                                     00660001
 90   ARR(I)=BLANK                                                      00670001
      IF(IV.LT.3) GOTO 110                                              00680001
      DO 100 I=1,IHH,5                                                  00690001
 100  ARR(I)=SYMB(6)                                                    00700001
      WRITE(6,103) LIN,IDH,IDV,DIG,ISUPP                                00710009
      IF(IHH.EQ.128) WRITE(6,101) DASH(IDASH),(ARR(I),I=1,IHH)          00720004
      IF(IHH.LT.128) WRITE(6,102) DASH(IDASH),(ARR(I),I=1,IHH)          00730004
 101  FORMAT(A4,128A1)                                                  00740004
 102  FORMAT(A4,1X,127A1)                                               00750004
 103  FORMAT(/' OUT2:  LINFLAG',I2,4X,'STEPS H,V',2I2,4X,'CUT-PAR',     00760009
     1       F8.2,4X,'ZERO-LINE SUPPRESS FLAG',I2)                      00770010
C-----------------------------------------------------------------      00780008
 110  MXD1=MXD+1                                                        00790001
      DO 200 J1=1,IV,IDV                                                00800001
      IF(IPARV.EQ.1) J=J1                                               00810001
      IF(IPARV.EQ.2) J=IV+1-J1                                          00820001
      IF(IPARV.EQ.3) J=(J1-1)*IH                                        00830001
      IF(IPARV.EQ.4) J=(IV-J1)*IH                                       00840001
      L=0                                                               00850001
      I0=0                                                              00860001
      DO 180 I1=1,IH,IDH                                                00870001
      I0=I0+1                                                           00880001
      IF(IPARH.EQ.1) I=I1                                               00890001
      IF(IPARH.EQ.2) I=IH+1-I1                                          00900006
      IF(IPARH.EQ.3) I=(I1-1)*IV                                        00910001
      IF(IPARH.EQ.4) I=(IH-I1)*IV                                       00920001
      JI=J+I                                                            00930001
      ARR(I0)=BLANK                                                     00940002
      DUM=POP(JI)                                                       00950002
      IF(DUM.LE.0.)GOTO 180                                             00960002
      DO 170 MD=1,MXD                                                   00970002
      IF(DUM.GE.AM(MD)) GOTO 175                                        00980002
 170  CONTINUE                                                          00990002
      GOTO 180                                                          01000002
 175  L=L+1                                                             01010002
      ARR(I0)=NU(MXD1-MD)                                               01020001
 180  CONTINUE                                                          01030001
      IF(L.EQ.0 .AND. ISUPP.EQ.1) GO TO 200                             01040001
      J2=IABS(J00+J1)                                                   01050003
      IF(IHH.LT.128) GOTO 192                                           01060003
      WRITE(6,190) J2,(ARR(I),I=1,I0)                                   01070003
 190  FORMAT(I4,128A1)                                                  01080001
      GOTO 200                                                          01090003
 192  WRITE(6,193) J2,(ARR(I),I=1,I0)                                   01100003
 193  FORMAT(I4,1X,127A1)                                               01110003
 200  CONTINUE                                                          01120001
      RETURN                                                            01130001
      END                                                               01140001
