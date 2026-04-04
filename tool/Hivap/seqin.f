      SUBROUTINE SEQIN(NROW,NSTEP,LDBM, CST,CLD)                        00020044
C
C  Unterprogramm zur Initialisierung spezifischer Variabler eines 
C  gegebenen, aktuellen Mutterkerns fuer Verdampfung (Separations-
C  energien, Zustandsdichteparameter etc.)
C  Stand: 3.5.1994
C 
      INCLUDE 'common.f'                                                00030051
      data r0wr,avwr,aswr,akwr/1.153,0.04543,0.1246,0.15231/                    
      data r0ts,avts,asts,akts/1.160,0.04385,0.1584,0.33200/                    
C       
C     I.Giese 2007
      SAVE
C
C                          INITIALISATIONS FOR CURRENT STEP             00040008
      IA=NSTEP                                                          00050043
      IZ=NROW                                                           00060043
      AVGF=0.                                                           00070011
      DO    J=1,JDIM                                                    00080052
        YPOPJ(J)=0.                                                     00090052
        YERJ(J)=0.                                                      00100052
        YFISJ(J)=0.                                                     00110052
      ENDDO                                                             00120052
      IF(IOPT.EQ.1) THEN                                                00130052
        DO    J=1,JDIM                                                  00140052
         FLAND(J)=0.                                                    00150052
        ENDDO                                                           00160052
      ENDIF   ! IOPT=1                                                  00170052
      DO 45 K=1,5                                                       00180052
      AVESP(K)=0.                                                       00190022
 45   SSP(K)=0.                                                         00200022
      DO 50 I=1,IPOPS                                                   00210010
 50   TABLE(I)=BLANK                                                    00220010
      IDUM=2*LN+1                                                       00230021
      DO 55 I=1,IDUM                                                    00240021
 55   DJN(I)=0.                                                         00250021
      IDUM=2*LP+1                                                       00260021
      DO 60 I=1,IDUM                                                    00270021
 60   DJP(I)=0.                                                         00280021
      IDUM=2*LA+1                                                       00290021
      DO 65 I=1,IDUM                                                    00300021
 65   DJA(I)=0.                                                         00310021
      IF(LPRINT.lt.3) then                                              00320052
        DO 70 K=1,KEN                                                   00330052
 70     SPECN(K)=0.                                                     00340052
        DO 72 K=1,KEP1                                                  00350052
 72     SPECP(K)=0.                                                     00360052
        DO 74 K=1,KEA1                                                  00370052
 74     SPECA(K)=0.                                                     00380052
        DO 76 K=1,KEF                                                   00390052
 76     SPECF(K)=0.                                                     00400052
        DO 78 K=1,KEG                                                   00410052
 78     SPECG(K)=0.                                                     00420052
        DO 79 K=1,KEG                                                   00430052
 79     SPECQ(K)=0.                                                     00440052
      endif  ! lprint<3                                                 00450052
C-----------------------------------------------------------------------00460001
C                                 Q-VALUES FOR CURRENT STEP             00470001
      Q(1)=0.                                                           00480053
      Q(5)=0.                                                           00490002
      DO 100 K=2,4                                                      00500002
 100  Q(K)=-BE(NROW,NSTEP,K-1)                                          00510002
C-----------------------------------------------------------------------00520001
C                                  LEVEL DENSITY PARAMETERS             00530001
      IF(BARFAC.LE.0) BARFAC=1.                                         00550002
      IF(DELT.LE.0.) DELT=12.                                           00560034
      IF(IPAIR.EQ.0)DELT=0.                                             00570016
      IF(TZERO.LE.0.) TZERO=0.2                                         00580002
C ------------------------                                              00600002
      IAF=KZ(NSTEP)+IN(NSTEP)                                           00610002
      AF=IAF                                                            00620000
      ZF=KZ(NSTEP)                                                      00630034
C ---------------------------------------------------                   00600002
      IF(LDBM.EQ.(-1)) then                                             00590042
        if(cld.le.0.) cld=1.                                                    
        AL(1)=ALS(IA,IZ)                                                01000016
        AL(2)=ALS(IA+1,IZ)                                              01010016
        AL(3)=ALS(IA,IZ+1)                                              01020016
        AL(4)=ALS(IA+2,IZ+2)                                            01030016
        AL(5)=CLD*AL(1)                                                 01040016
      endif   ! Ldbm=-1                                                         
c     -----------------------                                                   
      IF(LDBM.eq.0) then                                                00640042
        if(cld.le.0.) cld=1.                                                    
        if(cst.le.0.) cst=10.                                                   
        AL(1)=AF/CST                                                    00650000
        AL(2)=(AF-1.)/CST                                               00660000
        AL(3)=AL(2)                                                     00670000
        AL(4)=(AF-4.)/CST                                               00680000
        AL(5)=CLD*AL(1)                                                 00690002
      endif   ! ldbm=0                                                  00700034
c     -----------------------                                                   
      if(ldbm.gt.0) then                                                        
        if(abs(CLD-1.).lt.0.01) then   ! wr formula                             
          r0ld=r0wr                                                             
          avld=avwr                                                             
          asld=aswr                                                             
          akld=akwr                                                             
        else        ! ts formula                                                
          r0ld=r0ts                                                             
          avld=avts                                                             
          asld=asts                                                             
          akld=akts                                                             
        endif                                                                   
        if(abs(CST-1.15).lt.0.25) r0ld=cst                                      
c                                                                               
        AL(1)=avld*r0ld**3*AF+asld*r0ld**2*(AF**0.666667)+              00710034
     1        akld*r0ld*(AF**0.33333)                                   00720034
        AF1=AF-1.                                                       00730034
        AL(2)=avld*r0ld**3*AF1+asld*r0ld**2*(AF1**0.666667)+            00740034
     1        akld*r0ld*(AF1**0.33333)                                  00750034
        AL(3)=AL(2)                                                     00760034
        AF1=AF-4.                                                       00770034
        AL(4)=avld*r0ld**3*AF1+asld*r0ld**2*(AF1**0.666667)+            00780034
     1        akld*r0ld*(AF1**0.33333)                                  00790034
        CALL DROPB (ZF,AF,INDEX,DPBAR,Y,1,0)                            00800037
        AL(5)=avld*r0ld**3*AF+asld*r0ld**2*(AF**0.666667)*BS+           00810034
     1        akld*r0ld*(AF**0.33333)*BK                                00820034
      endif   ! ldbm>0                                                          
c     -----------------------                                                   
 144  IF(IPAIR.ne.4) then                                               00840043
        ODD1=MOD(KZ(NSTEP),2)                                           00850043
        ODD2=MOD(IN(NSTEP),2)                                           00860000
        DE=ABS(DELT    )/SQRT(AF)                                       00870000
        DELTA(1)=(1.-ODD1-ODD2)*DE                                      00880000
        DELTA(3)=(ODD1-ODD2)*DE                                         00890000
        DELTA(2)=(ODD2-ODD1)*DE                                         00900000
        DELTA(4)=DELTA(1)                                               00910000
        DELTA(5)=DELTA(1)                                               00920002
        IF(IPAIR   .eq.2 ) then                                         00930016
          DO     I=1,5                                                  00940053
          DELTA(I)=DELTA(I)-DE                                          00950053
          enddo                                                         00960053
        endif ! ipair=2                                                         
      endif  ! ipair not 4                                                      
C -------------------------------                                       00970001
      IF(ipair.eq.4 .or. LDBM.eq.(-1)) then                             00990043
        DELTA(1)=DELTAS(IA,IZ)                                          01050043
        DELTA(2)=DELTAS(IA+1,IZ)                                        01060016
        DELTA(3)=DELTAS(IA,IZ+1)                                        01070016
        DELTA(4)=DELTAS(IA+2,IZ+2)                                      01080016
        DELTA(5)=DELTAS(IA,IZ)                                          01090019
      endif  ! ipair=4                                                          
c     ---------------------------------------------------------                 
 240  IF(LPRINT.LT.2) WRITE(6,245) IPAIR,DELT,CST,AL,DELTA              01100034
 245  FORMAT(' SEQIN:PAIR,DELT,CST',I4,F7.2,F8.3,                       01110047
     1' AL=',5F7.2,' DELTA=',5F6.2)                                     01120047
      IF(LPRINT.LT.2 .AND. LDBM.GE.1)WRITE(6,250)CLD,BS,BK,DPBAR,AF,ZF,Y01130038
 250  FORMAT(' AF/AN=',F8.3,2X,' BS,BK=',2F9.4,' DPBAR=',F7.2,          01140047
     1   ' A=',F5.0,' Z=',F5.0,' Y=',F8.3)                              01150047
c     ---------------------------------------------------------------           
c                                     shell corrections                         
      IF(ISHELL.gt.0) then                                              01160023
        SHELK(1)=SHELLS(IA,IZ)                                          01170016
        SHELK(2)=SHELLS(IA+1,IZ)                                        01180015
        SHELK(3)=SHELLS(IA,IZ+1)                                        01190015
        SHELK(4)=SHELLS(IA+2,IZ+2)                                      01200015
        SHELK(5)=0.                                                     01210015
        IF(ISHELF.EQ.1) SHELK(5)=SHELK(1)                               01220016
      IF(ISHELF.EQ.2) SHELK(5)=BFS(NSTEP,NROW)-BFLDM(NSTEP,NROW)*BARFAC+01230031
     1                         SHELLS(NSTEP,NROW)                       01240017
        IF(LPRINT.LT.3)                                                 01250023
     1    WRITE(6,265) SHELK                                            01260023
 265    FORMAT(' SHELL CORRECTIONS',4X,5F8.1)                           01270013
      endif ! shell corections                                                  
c     --------------------------------------------------------------            
 300  RETURN                                                            01280012
      END                                                               01290001
