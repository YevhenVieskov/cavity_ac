SUBROUTINE INIT
use gridar
use phys
use rel
use condif
      IMPLICIT NONE

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

!      common /rel/ relatief
!      real relatief      

!      COMMON /CONDIF/ RU,RV,divuv,divmax
!      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
!      INTEGER divuv

!c local variables
      REAL POld,verschil,absoluut
      INTEGER i,j

      verschil=0.0
      absoluut=0.0
      DO 5 J=2,JM1Us
      DO 5 I=2,IM1Us
         verschil = verschil +&
                 (UN(I,J)-U(I,J))**2 + (VN(I,J)-V(I,J))**2
         absoluut = absoluut + U(I,J)**2 + V(I,J)**2
 5    CONTINUE   
      verschil=sqrt(verschil/(IM2Us*JM2Us))
      absoluut=sqrt(absoluut/(IM2Us*JM2Us))
      relatief =  verschil / (absoluut + 1e-10)
      relatief2 =  verschil / (absoluut + 1e-10)
      DO 10 J=1,JMaxUs
         UN(0,J)=U(0,J)
         U(0,J)=0.0
      DO 10 I=1,IMaxUs
         UN(I,J) = U(I,J)
         VN(I,J) = V(I,J)
         U(I,J) = 0.0
         V(I,J) = 0.0
         POLD = P(I,J)
!c         P(I,J)=2.0*P(I,J)-PN(I,J)
         PN(I,J)=POLD
 10   CONTINUE
      DO 11 I=1,IMaxUs
        VN(I,0)=V(I,0)
        V(I,0)=0.0
 11   CONTINUE

!c-- restore non-zero boundary conditions
      !CALL BCBND

      RETURN
      END
