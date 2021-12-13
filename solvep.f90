SUBROUTINE SOLVEP
use gridar
use phys
use orga
use times
use orga2
      IMPLICIT NONE

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX,Y,YJ,DelY,RDX,RDY,
!     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
 
!      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
!      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
!     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
!     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

!      COMMON /ORGA/ NF
!      INTEGER NF(IMAX,JMAX)

!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

!c local variables
      INTEGER I,J
	  REAL(8) BETA
!C
!      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
!     >               StrtP
!      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
!      REAL  Epsi,OmStrt,Alpha,NrmRhs

!      CALL COEFF

!      IF (IMilu.EQ.0) THEN
!c--- solve Poisson equation with SOR
!          CALL SOR
!      ELSE
!c--- solve with MILU
 !        CALL MILU
!      ENDIF

!C
!C**** VELOCITY-UPDATE FOR PRESSURE GRADIENT
!C
!      DO 180 J=1,JM1Us
!      DO 180 I=1,IM1Us
!        IF (NF(I,J).LE.OU .AND. NF(I+1,J).LE.OU)&
!         U(I,J)=U(I,J)-2.0*DelT*(P(I+1,J)-P(I,J))/(DelX(I+1)+DelX(I))
!        IF (NF(I,J).LE.OU .AND. NF(I,J+1).LE.OU)&
!         V(I,J)=V(I,J)-2.0*DelT*(P(I,J+1)-P(I,J))/(DelY(J+1)+DelY(J))
! 180  CONTINUE

      !CALL BCBND
	  BETA=1
      DO I=1,IM1Us
        DO J=1,JM1Us
		  P(I,J)=PN(I,J)-DelT2*BETA*((U(I,J)-U(I-1,J))*DelY(J)+(V(I,J)-V(I,J-1))*DelX(J))*RDX(I)*RDY(J)
		END DO
	  END DO
	  RETURN
	  END
