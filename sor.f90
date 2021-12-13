SUBROUTINE SOR
use gridar
use phys
use orga2
use coefp
use times
      IMPLICIT NONE

      !INTEGER IMAX,JMAX
      !PARAMETER (IMAX=130, JMAX=130)      

      !COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
     !1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
     ! REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     !1   RDX(IMAX),RDY(JMAX),CX,CY
     ! INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

     ! COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
     ! REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     !1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     !2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

     ! COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     !>               StrtP
     ! INTEGER IMilu,Iter,ItMax,ItSum,StrtP
     ! REAL  Epsi,OmStrt,Alpha,NrmRhs      

      !COMMON /COEFP/ CN, CS, CE, CW, CC, DIV
     ! REAL CN(IMAX,JMAX), CS(IMAX,JMAX), CE(IMAX,JMAX),
     !1     CW(IMAX,JMAX), CC(IMAX,JMAX), DIV(IMAX,JMAX)

      !COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      !REAL  T, DelT, DelTMx, TFin, TStart
      !INTEGER Cycle

!c local variables
      REAL OMEG,DELTA,DIFF,Epsi2,DEL0
      INTEGER i,j

!C     EUCLIDEAN NORM IS USED and SCALING WITH RHS
!C october 2003: only relative criterion 
!C      Epsi2=Epsi*SQRT(FLOAT(IM2Us*JM2Us))
!C      Epsi2=AMAX1(Epsi2, Epsi*NrmRhs)
       Epsi2=Epsi*NrmRhs

      OMEG=OmStrt
      Iter=0

      IF (StrtP.eq.0) then
        DO 10 j=1,JMaxUs
        DO 10 i=1,IMaxUs
	   p(i,j)=0.0
 10	CONTINUE
      ENDIF

 140  Iter=Iter+1
      DELTA=0.0

      DO 150 j=2,Jm1Us
      DO 150 i=2,Im1Us        
         DIFF = - P(i,j) + DIV(i,j)&
                - CW(i,j)*P(i-1,j) - CE(i,j)*P(i+1,j )&
                - CS(i,j)*P(i,j-1) - CN(i,j)*P(i,j+1 )
         P(i,j) = P(i,j) + OMEG*DIFF
         DELTA = DELTA + DIFF**2
         IF (Iter.EQ.1) DEL0=DELTA
 150  CONTINUE

      DELTA=SQRT(DELTA)
      
      IF (DELTA.GT.Epsi2.and.Iter.LT.ItMax) goto 140

      ItSum = ItSum + Iter

!c      IF (DELTA.GT.0.1*NrmRhs .AND. DELTA.GT.0.1*DEL0) THEN  
!c           WRITE(6,'(A17,2E13.4,1F8.4,3E13.4)')
!c     >       'slow convergence ', NrmRhs, DEL0, OMEG,
!c     >         DELTA, DELTA/NrmRhs, DELTA/DEL0
!c      ENDIF

      RETURN
      END

