      SUBROUTINE DROPB(Z,A,INDEX,DPBAR,Y,LDM,IPRNT) 
C
C  Droplet-Modell (frei nach Myers)
C  dient zur Berechnung des Oberflaechenverhaeltnisses von Sattelpunkt-
C  und Gleichgewichtskonfiguration, was a_f/a_n beeinflusst;
C  ruft DROPL auf.
C  Stand: 3.5.1994
C 
C       DROPLET BARRIER                                                 00050000
C    /* INDEX = 0 NUCLEUS TOO LIGHT FOR Y-FAMILY                   */   00060000
C    /*       = 1 SOLN WITH Y > .44                                */   00070000
C    /*       = 2 SOLN WITH Y = .44                                */   00080000
C    /*       = 3 SOLN WITH Y < .44                                */   00090000
C    /*       = 4 SOLN FOR HEAVY NUCLEI USING FULL Y_FAMILY SEARCH */   00100000
C    /*       = 5 FISSION UNSTABLE                                 */   00110000
C    /*       = 6 SOMETHING WRONG                                  */   00120000
C                                                                       00130000
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(36),BCA(36),BSA(36),BKA(36),BRA(36),BVA(36),BWA(36)  00140000
      DIMENSION YA(3 ), AA(3 )                                          00150000
      COMMON/BEES/BS,BC,BK,BR,BV,BW                                     00160000
      DATA                                                              00170000
     1  XA                                                              00180000
     2   /.30, .32, .34, .36, .38,                                      00190000
     3    .40, .42, .44, .46, .48,                                      00200000
     4    .50, .52, .54, .56, .58,                                      00210000
     5    .60, .62, .64, .66, .68,                                      00220000
     6    .70, .72, .74, .76, .78,                                      00230000
     7    .80, .82, .84, .86, .88,                                      00240000
     8    .90, .92, .94, .96, .98,                                      00250000
     9   1.00/                                                          00260000
       DATA                                                             00270000
     1  BCA                                                             00280000
     2   /.82500, .82332, .82173, .82023, .81881,                       00290000
     3    .81747, .81621, .81503, .81393, .81293,                       00300000
     4    .81206, .81133, .81081, .81056, .81073,                       00310000
     5    .81155, .81347, .81749, .82584, .84250,                       00320000
     6    .86664, .89017, .91008, .92667, .94060,                       00330000
     7    .95238, .96239, .97088, .97807, .98409,                       00340000
     8    .98905, .99303, .99609, .99827, .99957,                       00350000
     9   1.00000/                                                       00360000
      DATA                                                              00370000
     1  BSA                                                             00380000
     2  /1.27314, 1.27418, 1.27522, 1.27627, 1.27732,                   00390000
     3   1.27837, 1.27941, 1.28042, 1.28141, 1.28235,                   00400000
     4   1.28320, 1.28394, 1.28450, 1.28477, 1.28458,                   00410000
     5   1.28352, 1.28126, 1.27619, 1.26532, 1.24296,                   00420000
     6   1.20963, 1.17623, 1.14717, 1.12229, 1.10085,                   00430000
     7   1.08224, 1.06604, 1.05195, 1.03974, 1.02927,                   00440000
     8   1.02044, 1.01319, 1.00750, 1.00338, 1.00086,                   00450000
     9   1.00000/                                                       00460000
      DATA                                                              00470000
     1BKA                                                               00480000
     2  /1.5968, 1.5917, 1.5867, 1.5817, 1.5768,                        00490000
     3   1.5718, 1.5668, 1.5617, 1.5564, 1.5509,                        00500000
     4   1.5452, 1.5389, 1.5321, 1.5243, 1.5151,                        00510000
     5   1.5036, 1.4886, 1.4668, 1.4216, 1.3727,                        00520000
     6   1.3014, 1.2418, 1.1951, 1.1576, 1.1267,                        00530000
     7   1.1009, 1.0790, 1.0609, 1.0458, 1.0328,                        00540000
     8   1.0224, 1.0141, 1.0079, 1.0035, 1.0009,                        00550000
     9   1.0000/                                                        00560000
      DATA                                                              00570000
     1BRA                                                               00580000
     2  / .4548,  .4538,  .4527,  .4516,  .4506,                        00590000
     3    .4495,  .4486,  .4476,  .4468,  .4462,                        00600000
     4    .4458,  .4458,  .4463,  .4477,  .4506,                        00610000
     5    .4560,  .4661,  .4860,  .5295,  .6247,                        00620000
     6    .7578,  .8660,  .9386,  .9845, 1.0119,                        00630000
     7   1.0264, 1.0323, 1.0326, 1.0292, 1.0241,                        00640000
     8   1.0181, 1.0123, 1.0071, 1.0032, 1.0007,                        00650000
     9   1.0000/                                                        00660000
      DATA                                                              00670000
     1BVA                                                               00680000
     2  / .7832,  .7810,  .7787,  .7765,  .7741,                        00690000
     3    .7718,  .7694,  .7669,  .7643,  .7617,                        00700000
     4    .7589,  .7561,  .7532,  .7502,  .7472,                        00710000
     5    .7444,  .7422,  .7419,  .7473,  .7682,                        00720000
     6    .8071,  .8462,  .8786,  .9048,  .9257,                        00730000
     7    .9426,  .9564,  .9673,  .9761,  .9832,                        00740000
     8    .9887,  .9930,  .9961,  .9983,  .9996,                        00750000
     9   1.0000/                                                        00760000
      DATA                                                              00770000
     1BWA                                                               00780000
     2  / .5782,  .5732,  .5685,  .5641,  .5601,                        00790000
     3    .5564,  .5531,  .5500,  .5474,  .5451,                        00800000
     4    .5434,  .5422,  .5410,  .5425,  .5447,                        00810000
     5    .5493,  .5583,  .5752,  .6096,  .6781,                        00820000
     6    .7673,  .8392,  .8899,  .9255,  .9503,                        00830000
     7    .9677,  .9797,  .9877,  .9930,  .9963,                        00840000
     8    .9983,  .9993,  .9998, 1.0000, 1.0000,                        00850000
     9   1.0000/                                                        00860000
C---------------------------------------------------------              00870000
C
C     I.Giese 2007
      SAVE
C
      DPBAR = 3.0E+30                                                   00880000
      Y = 3.0E+30                                                       00890000
                                                                        00900000
C        /* CHECK TO SEE IF THE SPHERE IS STABLE */                     00910000
      BS = BSA(35)                                                      00920000
      BC = BCA(35)                                                      00930000
      BK = BKA(35)                                                      00940000
      BR = BRA(35)                                                      00950000
      BV = BVA(35)                                                      00960000
      BW = BWA(35)                                                      00970000
      IDUM=0                                                            00980000
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )                           00990002
      CALL DROPL (Z,A,1,DROP,DEL,EPS,LDM,IPRNT)                         01000002
      IF((F2-DROP).GT.0.)  GO TO      2                                 01010000
      INDEX = 5                                                         01020000
C     INDEX = 5                       /* SPHERE IS UNSTABLE     */      01030000
      RETURN                                                            01040000
C        /* IF ICODE = 0 LOOK FOR MAX., IF ICODE = 1 LOOK FOR MIN. */   01050000
 2    ICODE = 0                                                         01060000
C        /* FOR HEAVY NUCLEI USE ALL BEES IPERT = 0             */      01070000
C        /* WHEN IPERT = 1 A PERTURBATION APPROACH IS USED      */      01080000
C        /* IF Z > 68 A NON-PERT CALCULATION IS ATTEMPTED       */      01090000
      IPERT = 0                                                         01100000
      IF(Z.GT.68)    GO TO      31                                      01110000
 30   IPERT = 1                                                         01120000
      BS = BSA(35)                                                      01130000
      BC = BCA(35)                                                      01140000
      BK = 1                                                            01150000
      BR = 1                                                            01160000
      BV = 1                                                            01170000
      BW = 1                                                            01180000
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )                           01190002
C                                                                       01200000
 31   DO 50 I=3,36                                                      01210000
         J = 37 - I                                                     01220000
         F1 = F2                                                        01230000
         BS = BSA(J)                                                    01240000
         BC = BCA(J)                                                    01250000
         IF(IPERT.EQ.1)    GO TO      3                                 01260000
         BK = BKA(J)                                                    01270000
         BR = BRA(J)                                                    01280000
         BV = BVA(J)                                                    01290000
         BW = BWA(J)                                                    01300000
 3       CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )                        01310002
         IF(ICODE.EQ.1)    GO TO      11                                01320000
C        /* LOOKING FOR A MAXIMUM */                                    01330000
         IF(I.EQ.25 .AND. IPERT.EQ.0) GOTO     30                       01340000
         IF((F2-F1).GT.0.)   GO TO      10                              01350000
         IF((I-24).LT.0)    GO TO      12                               01360000
         IF((I-24).EQ.0)    GO TO      13                               01370000
         IF((I-24).GT.0)    GO TO      14                               01380000
C        /* GOOD MAXIMUM */                                             01390000
 12      INDEX = 3                                                      01400000
         IF(IPERT.EQ.0)    INDEX = 4                                    01410000
         GO TO      15                                                  01420000
C        /* CONTINUE SEARCH LOOKING FOR MINIMUM */                      01430000
 13      ICODE = 1                                                      01440000
         GO TO      10                                                  01450000
C        /* THIS SHOULD NOT HAPPEN */                                   01460000
 14      INDEX = 6                                                      01470000
         RETURN                                                         01480000
C        /* LOOKING FOR MINIMUM */                                      01490000
 11      IF((F2-F1).LT.0.)   GO TO      10                              01500000
         IF(I.LT.14)    GO TO      14                                   01510000
         INDEX = 1                                                      01520000
         GO TO      15                                                  01530000
 10   CONTINUE                                                          01540000
 50   CONTINUE                                                          01550000
 17   IF(ICODE.EQ.0)    GO TO      16                                   01560000
C     /* TOO LIGHT, NO MINIMUM */                                       01570000
      INDEX = 0                                                         01580000
      RETURN                                                            01590000
C        /* NO MAX FOUND SO IT MUS BE AT Y = .44 */                     01600000
 16   INDEX = 2                                                         01610000
      J = 13                                                            01620000
 15   J = J + 1                                                         01630000
      BS = BSA(J)                                                       01640000
      BC = BCA(J)                                                       01650000
      BR = BRA(J)                                                       01660000
      BK = BKA(J)                                                       01670000
      BV = BVA(J)                                                       01680000
      BW = BWA(J)                                                       01690000
      RI2=(A-2.*Z)/A                                                    01700000
      RI2=RI2*RI2                                                       01710000
      XLDM= .7053 *Z*Z/(2.*A*17.9439*(1.-1.78  *RI2))                   01720000
      IF(IPRNT.EQ.1)WRITE(6,60) XA(J),XLDM                              01730000
 60   FORMAT(/' SADDLEPOINT X  =',F10.3,4X,'XLDM=',F10.3)               01740000
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IPRNT)                           01750002
      DPBAR = F2 - DROP                                                 01760000
      Y = 1.E0 - XA(J)                                                  01770000
      IF(IPRNT.EQ.1)WRITE(6,70)XA(J),DPBAR                              01780000
      YA(2)=DPBAR                                                       01790000
      IREP=0                                                            01800000
      JJ=J-1                                                            01810000
      IF(JJ.LT.1)RETURN                                                 01820000
 65   BS = BSA(JJ)                                                      01830000
      BC = BCA(JJ)                                                      01840000
      BR = BRA(JJ)                                                      01850000
      BK = BKA(JJ)                                                      01860000
      BV = BVA(JJ)                                                      01870000
      BW = BWA(JJ)                                                      01880000
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )                           01890002
      DPBAR=F2-DROP                                                     01900000
      IF(IREP.EQ.0)YA(1)=DPBAR                                          01910000
      IF(IPRNT.EQ.1)WRITE(6,70)XA(JJ),DPBAR                             01920000
 70   FORMAT(' X=',F10.3,4X,'BARRIER=',F10.2)                           01930000
      IREP=IREP+1                                                       01940000
      IF(IREP.GT.1) GOTO 80                                             01950000
      JJ=J+1                                                            01960000
      IF(JJ.GT.36)RETURN                                                01970000
      GOTO 65                                                           01980000
 80   YA(3)=DPBAR                                                       01990000
      J1=J-1                                                            02000000
      CALL PARABO(XA(J1 ),YA,AA,X0,BAR)                                 02010000
      IF(IPRNT.EQ.1) WRITE(6,90)X0,BAR                                  02020000
 90   FORMAT(' PARABOLA X=',F10.3,4X,'BARRIER=',F10.2)                  02030000
      J2=J+1                                                            02040000
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BSA(J1),BSA(J),BSA(J2),BS)        02050000
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BCA(J1),BCA(J),BCA(J2),BC)        02060000
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BRA(J1),BRA(J),BRA(J2),BR)        02070000
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BKA(J1),BKA(J),BKA(J2),BK)        02080000
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BVA(J1),BVA(J),BVA(J2),BV)        02090000
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BWA(J1),BWA(J),BWA(J2),BW)        02100000
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IPRNT)                           02110002
      DPBAR = F2 - DROP                                                 02120000
      Y = 1.E0 - X0                                                     02130000
      RETURN                                                            02140000
C                                                                       02150000
      END                                                               02160000
