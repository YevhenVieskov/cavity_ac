SUBROUTINE DIVERR
use gridar
use phys
use condif
!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX, RDY,
!     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

!      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
!      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
!     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
!     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

!      COMMON /CONDIF/ RU,RV,divuv,divmax
!      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
!      INTEGER divuv

      INTEGER i,j
      REAL divu

      divmax=0.0
      DO 15 J=2,JM1Us
      DO 15 I=2,IM1Us
	 divu=RDX(I)*(U(I,J)-U(I-1,J)) +& 
            RDY(J)*(V(I,J)-V(I,J-1))
	 divmax=MAX(ABS(divu),divmax)
   15 CONTINUE

       RETURN
       END