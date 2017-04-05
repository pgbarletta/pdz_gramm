      SUBROUTINE apc_wrapper(m, n, mtx_a, mtx_b, indices)
      implicit none

      integer i, j, k, m, n, Z, Nmax
      parameter (Nmax=2700)
      integer indices(n), pjtion (3*Nmax, 3*Nmax), tmp_indices(3*Nmax)
      real*8 mtx_a(m,n), mtx_b(m,n), tmp (m,n)


      do i=1, n
         do j=1, n
            tmp(i,j)=0.0d0
            do k=1, m
               tmp(i, j) = tmp(i,j) + mtx_a(k, i) * mtx_b(k, j)
            enddo
         enddo
      enddo

      do i=1, n
        do j=1, n
           pjtion(i,j)=int(tmp(i,j)**2*1.d5)
        enddo
      enddo

      do i=1, n
         do j=1, n
            if((j.lt.(i-4)).or.(j.gt.(i+4))) then
              pjtion(i,j)=-1*ifix(sngl(1.d5))
            endif
         enddo
      enddo

      do i=1, n
         do j=1, n
            pjtion(i,j)= -1 * pjtion(i,j)
         enddo
      enddo


      call apc(n, pjtion, tmp_indices, z)

      do i=1, n
        indices(i) = tmp_indices(i)
      enddo

      RETURN 
      END


C     ******* SAMPLE CALLING PROGRAM FOR SUBROUTINE APC	   *******
C     ***     (MIN-COST	ASSIGNMENT PROBLEM)		       ***
C     ***						       ***
C     ***     THE PROGRAM IS BASED ON THE PAPER		       ***
C     ***     G. CARPANETO, S. MARTELLO, P. TOTH "ALGORITHMS   ***
C     ***       AND CODES FOR THE ASSIGNMENT PROBLEM",	       ***
C     ***	ANNALS OF OPERATIONS RESEARCH 7, 1988.	       ***
C     ***						       ***
C     ***     ALL THE SUBROUTINES ARE WRITTEN IN AMERICAN      ***
C     ***	STANDARD FORTRAN AND ARE ACCEPTED BY THE       ***
C     ***	PFORT VERIFIER.				       ***
C     ***						       ***
C     ***     QUESTIONS	AND COMMENTS SHOULD BE DIRECTED	TO     ***
C     ***     SILVANO MARTELLO AND PAOLO TOTH		       ***
C     ***     D.E.I.S.,	UNIVERSITA' DI BOLOGNA,                ***
C     ***     VIALE RISORGIMENTO 2,                            ***
C     ***     40136, BOLOGNA, ITALY.                           ***
C     ************************************************************
C
c      INTEGER A(260,260),F(260),Z
c
c      double precision b(260,260)
c
c      READ(5,10) N
c   10 FORMAT(20I4)
c
cc      DO 20 I=1,N
c        READ(5,10) (A(I,J),J=1,N)
cc   20 CONTINUE
c
c      DO 20 I=1,N
c        READ(5,*) (b(I,J),J=1,N)
c   20 CONTINUE
c
c      do i=1,n
c         do j=1,n
cc            a(i,j)=ifix(sngl(b(i,j)**2*1.d4))
c            a(i,j)=int((b(i,j)**2*1.d5))
c         enddo
c      enddo
c
c      DO I=1,N
cc        print*, (a(I,J),J=1,N)
c      enddo
c
cc     ----------------------
c
c      do i=1,n
c         do j=1,n
c            a(i,j)=-1*a(i,j)
c         enddo
c      enddo


c      CALL APC(N,A,F,Z)
c      WRITE(6,30) Z
c   30 FORMAT(26H  COST OF THE ASSIGNMENT =,I10/)
c      WRITE(6,40) (F(I),I=1,N)
c   40 FORMAT(11H ASSIGNMENT,20I4)
c
c      do i=1,n
c         print*,i,b(i,f(i))
c      enddo
c

c
c      STOP
c      END
      SUBROUTINE APC(N,A,F,Z)
C
C SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
C
C HUNGARIAN METHOD. COMPLEXITY O(N**3).
C
C
C MEANING OF THE INPUT PARAMETERS:
C N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
C A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
C ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
C
C MEANING OF THE OUTPUT PARAMETERS:
C F(I) = COLUMN ASSIGNED TO ROW  I .
C Z    = COST OF THE OPTIMAL ASSIGNMENT =
C      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
C
C ALL THE PARAMETERS ARE INTEGERS.
C VECTOR  F  MUST BE DIMENSIONED AT LEAST AT  N , MATRIX  A
C AT LEAST AT  (N,N) . AS CURRENTLY DIMENSIONED, THE SIZE
C LIMITATION IS  N .LE. 260 . IN ALL THE SUBROUTINES, THE
C INTERNAL VARIABLES WHICH ARE PRESENTLY DIMENSIONED AT
C 260 MUST BE DIMENSIONED AT LEAST AT  N .
C
C THE ONLY MACHINE-DEPENDENT CONSTANT USED IS  INF (DEFINED
C BY THE FIRST EXECUTABLE STATEMENT OF THIS SUBROUTINE). INF
C MUST BE SET TO A VERY LARGE INTEGER VALUE.
C
C THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
C LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
C MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
C THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
C GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
C CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
C RESEARCH 7, 1988.
C
C SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
C PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
C UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
C
C MEANING OF THE MAIN INTERNAL VARIABLES:
C FB(J) = ROW ASSIGNED TO COLUMN  J .
C M     = NUMBER OF INITIAL ASSIGNMENTS.
C U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
C V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
C
C APC NEEDS THE FOLLOWING SUBROUTINES: INCR
C                                      INIT
C                                      PATH
C
C ALL THE SUBROUTINES ARE WRITTEN IN AMERICAN NATIONAL
C STANDARD FORTRAN AND ARE ACCEPTED BY THE PFORT VERIFIER.
C
C
C THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
C
c      parameter (Nmax=1000)
      parameter (Nmax=2700) 
      INTEGER A(3*Nmax,3*Nmax),F(3*Nmax),Z
      INTEGER U,V,FB
      COMMON U(3*Nmax),V(3*Nmax),FB(3*Nmax)
      INF = 10**9
C SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
      CALL INIT(N,A,F,M,INF)
      IF ( M .EQ. N ) GO TO 20
C SOLUTION OF THE REDUCED PROBLEM.
      DO 10 I=1,N
        IF ( F(I) .GT. 0 ) GO TO 10
C DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
        CALL PATH(N,A,I,F,INF,J)
C ASSIGNMENT OF ROW  I  AND COLUMN  J .
        CALL INCR(F,J)
   10 CONTINUE
C COMPUTATION OF THE SOLUTION COST  Z .
   20 Z = 0
      DO 30 K=1,N
        Z = Z + U(K) + V(K)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE INCR(F,J)
C
C ASSIGNMENT OF COLUMN  J .
C
c      parameter (Nmax=1000)
      parameter (Nmax=2700)
      INTEGER F(3*Nmax)
      INTEGER U,V,FB,RC
      COMMON U(3*Nmax),V(3*Nmax),FB(3*Nmax),RC(3*Nmax)
   10 I = RC(J)
      FB(J) = I
      JJ = F(I)
      F(I) = J
      J = JJ
      IF ( J .GT. 0 ) GO TO 10
      RETURN
      END
      SUBROUTINE INIT(N,A,F,M,INF)
C
C SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
C
C P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
C
c      parameter (Nmax=1000)
      parameter (Nmax=2700)
      INTEGER A(3*Nmax,3*Nmax),F(3*Nmax)
      INTEGER U,V,FB,P,R
      COMMON U(3*Nmax),V(3*Nmax),FB(3*Nmax),P(3*Nmax)
C PHASE 1 .
      M = 0
      DO 10 K=1,N
        F(K) = 0
        FB(K) = 0
   10 CONTINUE
C SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
      DO 40 J=1,N
        MIN = INF
        DO 30 I=1,N
          IA = A(I,J)
          IF ( IA .GT. MIN ) GO TO 30
          IF ( IA .LT. MIN ) GO TO 20
          IF ( F(I) .NE. 0 ) GO TO 30
   20     MIN = IA
          R = I
   30   CONTINUE
        V(J) = MIN
        IF ( F(R) .NE. 0 ) GO TO 40
C ASSIGNMENT OF COLUMN  J  TO ROW  R .
        M = M + 1
        FB(J) = R
        F(R) = J
        U(R) = 0
        P(R) = J + 1
   40 CONTINUE
C PHASE 2 .
C SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
      DO 110 I=1,N
        IF ( F(I) .NE. 0 ) GO TO 110
        MIN = INF
        DO 60 K=1,N
          IA = A(I,K) - V(K)
          IF ( IA .GT. MIN ) GO TO 60
          IF ( IA .LT. MIN ) GO TO 50
          IF ( FB(K) .NE. 0 ) GO TO 60
          IF ( FB(J) .EQ. 0 ) GO TO 60
   50     MIN = IA
          J = K
   60   CONTINUE
        U(I) = MIN
        JMIN = J
        IF ( FB(J) .EQ. 0 ) GO TO 100
        DO 80 J=JMIN,N
          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
          R = FB(J)
          KK = P(R)
          IF ( KK .GT. N ) GO TO 80
          DO 70 K=KK,N
            IF ( FB(K) .GT. 0 ) GO TO 70
            IF ( A(R,K) - U(R) - V(K) .EQ. 0 ) GO TO 90
   70     CONTINUE
          P(R) = N + 1
   80   CONTINUE
        GO TO 110
C REASSIGNMENT OF ROW  R  AND COLUMN  K .
   90   F(R) = K
        FB(K) = R
        P(R) = K + 1
C ASSIGNMENT OF COLUMN  J  TO ROW  I .
  100   M = M + 1
        F(I) = J
        FB(J)= I
        P(I) = J + 1
  110 CONTINUE
      RETURN
      END
      SUBROUTINE PATH(N,A,II,F,INF,JJ)
C
C DETERMINATION OF AN AUGMENTING PATH STARTING FROM
C UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
C JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
C
C MEANING OF THE MAIN INTERNAL VARIABLES:
C LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
C PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
C         LABELLED AND NOT EQUAL TO  FB(J) ).
C RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
C         ALTERNATING PATH.
C UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
C
c      parameter (Nmax=1000)
      parameter (Nmax=2700)
      INTEGER A(3*Nmax,3*Nmax),F(3*Nmax),Z
      INTEGER PI(3*Nmax),LR(3*Nmax),UC(3*Nmax)
      INTEGER U,V,FB,RC,R
      COMMON U(3*Nmax),V(3*Nmax),FB(3*Nmax),RC(3*Nmax)
C INITIALIZATION.
      LR(1) = II
      DO 10 K=1,N
        PI(K) = A(II,K) - U(II) - V(K)
        RC(K) = II
        UC(K) = K
   10 CONTINUE
      NUC = N
      NLR = 1
      GO TO 40
C SCANNING OF THE LABELLED ROWS.
   20 R = LR(NLR)
      DO 30 L=1,NUC
        J = UC(L)
        IA = A(R,J) - U(R) - V(J)
        IF ( IA .GE. PI(J) ) GO TO 30
        PI(J) = IA
        RC(J) = R
   30 CONTINUE
C SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
   40 DO 50 L=1,NUC
        J = UC(L)
        IF ( PI(J) .EQ. 0 ) GO TO 100
   50 CONTINUE
C UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
      MIN = INF
      DO 60 L=1,NUC
        J = UC(L)
        IF ( MIN .GT. PI(J) ) MIN = PI(J)
   60 CONTINUE
      DO 70 L=1,NLR
        R = LR(L)
        U(R) = U(R) + MIN
   70 CONTINUE
      DO 90 J=1,N
        IF ( PI(J) .EQ. 0 ) GO TO 80
        PI(J) = PI(J) - MIN
        GO TO 90
   80   V(J) = V(J) - MIN
   90 CONTINUE
      GO TO 40
  100 IF ( FB(J) .EQ. 0 ) GO TO 110
C LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF
C COLUMN  J .
      NLR = NLR + 1
      LR(NLR) = FB(J)
      UC(L) = UC(NUC)
      NUC = NUC - 1
      GO TO 20
C DETERMINATION OF THE UNASSIGNED COLUMN  J .
  110 JJ = J
      RETURN
      END
 
