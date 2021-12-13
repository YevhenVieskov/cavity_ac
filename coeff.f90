      SUBROUTINE COEFF
	  use gridar
	  use phys
	  use orga
	  use times
      use coefp
	  use coef2
	  use bndcon
	  use orga2
	  use condif
      IMPLICIT NONE

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX,RDY,
!     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!c
!      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
!      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
!     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
!     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

!      COMMON /ORGA/ NF
!      INTEGER NF(IMAX,JMAX)

!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

!      COMMON /COEFP/ CN, CS, CE, CW, CC, DIV
!      REAL CN(IMAX,JMAX), CS(IMAX,JMAX), CE(IMAX,JMAX),
!     1     CW(IMAX,JMAX), CC(IMAX,JMAX), DIV(IMAX,JMAX)

!      COMMON /coef2/ CWI,CEI,CNJ,CSJ
!      REAL CWI(IMAX),CEI(IMAX),CNJ(JMAX),CSJ(JMAX)

!      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
!     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
!      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
!     1     , POutN, POutE, POutS, POutW
!      INTEGER WB, WL, WR, WT

!      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
!     >               StrtP
!      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
!      REAL  Epsi,OmStrt,Alpha,NrmRhs

!      COMMON /CONDIF/ RU, RV,divuv,divmax
!      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
!      INTEGER divuv

!c local variables

      INTEGER I,J
      REAL schaal
!C
!C***  only first time step: set coefficients for Poisson equation
!C     scaling by SCHAAL makes coefficient matrix symmetric 
!C
      IF (t_cycle.EQ.1) THEN

         DO 100 J=2,JM1Us
         DO 100 I=2,IM1Us
            schaal=-DelX(I)*DelY(J)
            IF (NF(I,J).EQ.F) THEN
!c--- full cell
               CS(I,J)=CSJ(J)*schaal
               CN(I,J)=CNJ(J)*schaal
               CE(I,J)=CEI(I)*schaal
               CW(I,J)=CWI(I)*schaal
               CC(I,J)=-CW(I,J)-CE(I,J)-CN(I,J)-CS(I,J)
!c
!c--- eliminate solid walls, inflow (u=prescribed) and outflow
!c
               IF (NF(I+1,J).GE.IN) CC(I,J)=CC(I,J)+CE(I,J)
               IF (NF(I+1,J).GE.OU) CE(I,J)=0.0
 
               IF (NF(I-1,J).GE.IN) CC(I,J)=CC(I,J)+CW(I,J)
               IF (NF(I-1,J).GE.OU) CW(I,J)=0.0

               IF (NF(I,J+1).GE.IN) CC(I,J)=CC(I,J)+CN(I,J)
               IF (NF(I,J+1).GE.OU) CN(I,J)=0.0
  
               IF (NF(I,J-1).GE.IN) CC(I,J)=CC(I,J)+CS(I,J)
               IF (NF(I,J-1).GE.OU) CS(I,J)=0.0

            ELSE
!c--- obstacle cell
               CC(I,J) = -schaal
               CN(I,J) = 0.0
               CS(I,J) = 0.0
               CE(I,J) = 0.0
               CW(I,J) = 0.0
               DIV(I,J) = 0.0
            ENDIF

100      CONTINUE

         IF (IMilu.EQ.0) THEN
!c--- scale diagonal at 1.0 for SOR
            DO 500 J=2,JM1Us
            DO 500 I=2,IM1Us
               CN(I,J)=CN(I,J)/CC(I,J)
               CS(I,J)=CS(I,J)/CC(I,J)
               CE(I,J)=CE(I,J)/CC(I,J)
               CW(I,J)=CW(I,J)/CC(I,J)
 500        CONTINUE
         ENDIF

      ENDIF
!c
      if(T_cycle==1) then
        !open(150,file="coeff.dat")
        !call PRNT(cw,"cw",150,IMaxUs,JMaxUs)
	  !call PRNT(ce,"ce",150,IMaxUs,JMaxUs)
	  !call PRNT(cs,"cs",150,IMaxUs,JMaxUs)
	  !call PRNT(cn,"cn",150,IMaxUs,JMaxUs)
	end if
!C*** each time step compute right-hand side of Poisson equation
!C
      DO 200 J=2,JM1Us
      DO 200 I=2,IM1Us
         IF (NF(I,J).EQ.F) THEN

!C.... only full cells are considered
!C

          IF (divuv.eq.1) THEN
             DIV(I,J)=( RDX(I)*(U(I,J)-U(I-1,J)) + &
                                  RDY(J)*(V(I,J)-V(I,J-1)) )/DelT  
          ELSE
             DIV(I,J)=( RDX(I)*(RU(I,J)-RU(I-1,J)) +& 
                                  RDY(J)*(RV(I,J)-RV(I,J-1)) )   
          ENDIF

!c--- outflow boundary condition p=POut (set in SETPAR)
!c
            IF (nf(i-1,j).EQ.OU) div(i,j)=div(i,j)-cwi(i)*POutW
            IF (nf(i+1,j).EQ.OU) div(i,j)=div(i,j)-cei(i)*POutE
            IF (nf(i,j-1).EQ.OU) div(i,j)=div(i,j)-csj(j)*POutS
            IF (nf(i,j+1).EQ.OU) div(i,j)=div(i,j)-cnj(j)*POutN
!c 
            schaal=-SIGN(1.0, cc(i,j)) *DelX(I)*DelY(J) 
            DIV(I,J)=DIV(I,J)*schaal

            IF (IMilu.EQ.0) DIV(I,J)=DIV(I,J)/CC(I,J)

         ENDIF
  200 CONTINUE 

!C*** Compute norm of right-hand side;
!C    to be used in SOR convergence criterion for Poisson equation.
!C
      NrmRhs=0.0
      DO J=2,JM1Us
      DO I=2,IM1Us
         NrmRhs=NrmRhs + DIV(I,J)*DIV(I,J)
      ENDDO
      ENDDO
      NrmRhs = SQRT(NrmRhs)

      RETURN
      END