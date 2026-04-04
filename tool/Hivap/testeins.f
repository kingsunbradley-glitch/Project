      PROGRAM TESTEINS
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(3)
      EINS=1.
      WRITE(*,*) EINS, 1., 1.0E0, 1.0D0
      A(1)=1.0
      A(2)=2.0
      A(3)=3.0
      N=5
      CALL TESTSUB(A,N,EINS)
      CALL TESTSUB(A,N,1.)
      CALL TESTSUB(A,N,1.0)
      CALL TESTSUB(A,N,1.0E0)
      CALL TESTSUB(A,N,1.0D0)
      END
