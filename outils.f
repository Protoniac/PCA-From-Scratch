      SUBROUTINE CALCUL_A(Q,LDQ,C,LDC,Y,LDY,A,LDA,p)
      IMPLICIT NONE
      INTEGER LDQ,LDC,LDY,LDA,I
      DOUBLE PRECISION Q(LDQ,*),C(LDC),Y(LDY),A(LDA),p
      CALL DGEMV('N',LDQ,LDC,(1/p)*1D0,Q,LDQ,C,1,0D0,A,1)
      DO I=1,LDA
        A(I) = Y(I) - A(I)
      END DO
      END SUBROUTINE
      
      SUBROUTINE CALCUL_D(C,LDC,H,LDH,D,LDD)
      IMPLICIT NONE 
      INTEGER I,LDC,LDH,LDD 
      DOUBLE PRECISION C(LDC), H(LDH), D(LDD)
      DO I=1,LDD
        D(I) = (C(I+1) - C(I))/(3*H(I))
      END DO
      END SUBROUTINE
      
      SUBROUTINE CALCUL_B(A,LDA,C,LDC,D,LDD,H,LDH,B,LDB)
      IMPLICIT NONE
      INTEGER I,LDA,LDC,LDD,LDH,LDB
      DOUBLE PRECISION A(LDA),C(LDC),D(LDD),H(LDH),B(LDB)
      DO I=1,LDB
        B(I) = (A(I+1) - A(I))/H(I) - C(I)*H(I) - D(I)*H(I)*H(I)
      END DO
      END SUBROUTINE
      
      SUBROUTINE ADDITION_MATRICE_CARRE(A,LDA,B,LDB,M,LDM,alpha)
      IMPLICIT NONE
      INTEGER LDA,LDB,LDM,I,J
      DOUBLE PRECISION alpha, A(LDA,*), B(LDB,*), M(LDM,*)
      DO I=1,LDA
        DO J=1,LDA
            M(I,J) = A(I,J) + alpha * B(I,J)
        END DO
      END DO
      END SUBROUTINE

      SUBROUTINE DESCENTE(U,LDU,B,LDB)
      IMPLICIT NONE
      DOUBLE PRECISION U(LDU,*), B(LDB), S
      INTEGER I, J,LDU,LDB
*f2py intent(inplace) B
      DO I=1,LDU
        S=0
        DO J=1,I-1
            S=S+U(I,J)*B(J)
        END DO
        B(I)=(B(I)-S)/U(I,I)
      END DO
      RETURN
      END SUBROUTINE


      SUBROUTINE REMONTEE(U,LDU,B,LDB)
      IMPLICIT NONE
      DOUBLE PRECISION U(LDU,*), B(LDB), S
      INTEGER I, J,LDU,LDB
*f2py intent(inplace) B
      DO I=LDU,1,-1
        S=0
        DO J=I+1,LDU
            S=S+U(J,I)*B(J)
        END DO
        B(I)=(B(I)-S)/U(I,I)
      END DO
      RETURN
      END SUBROUTINE
      
      SUBROUTINE IDENTITE(A,LDA)
      IMPLICIT NONE
      INTEGER LDA,I,J
      DOUBLE PRECISION A(LDA,LDA)
      DO I = 1,LDA
        DO J = 1,LDA
            IF (I==J) THEN
                A(I,I) = 1D0
            ELSE 
                A(I,J) = 0D0
            END IF
        END DO
      END DO
      END SUBROUTINE
      
      SUBROUTINE INVERSION(A,LDA,INVERSE,LDINVERSE,LDCINVERSE)
      IMPLICIT NONE
      INTEGER LDA,LDINVERSE,I,LDCINVERSE
      DOUBLE PRECISION A(LDA,*), LOWER(LDA,LDA), 
     $                 INVERSE(LDINVERSE,LDCINVERSE)
      CALL CHOLESKY(A,LDA,LOWER,LDA)
      CALL IDENTITE(INVERSE,LDINVERSE)
      DO I = 1,LDA
        CALL DESCENTE(LOWER,LDA,INVERSE(:,I),LDINVERSE)
        CALL REMONTEE(LOWER,LDA,INVERSE(:,I),LDINVERSE)
      END DO
      END SUBROUTINE
      
      DOUBLE PRECISION FUNCTION NORME2(U,LDU,V,LDV)
      INTEGER LDU,I
      DOUBLE PRECISION U(LDU),V(LDV)
      NORME2 = 0
      DO I=1,LDU
        NORME2 = NORME2 + U(I) * V(I)
      END DO
      RETURN
      END FUNCTION
