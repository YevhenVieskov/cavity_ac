SUBROUTINE BCBND
use gridar
use phys
use times
use orga
use bndcon
      IMPLICIT NONE

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX, RDY,
!     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

!      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
!      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
!     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
!     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

!      COMMON /ORGA/ NF
!      INTEGER NF(IMAX,JMAX)

!      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
!     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
!      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
!     1     , POutN, POutE, POutS, POutW
!      INTEGER WB, WL, WR, WT

!c  local variables
      INTEGER i,j

!c... velocity conditions: SgnL=1 slip wall; SgnL=-1 no-slip wall

      DO 100 J=1,JMaxUs
!c-- western boundary
         V(1,J)=SgnL*V(2,J)
         IF (nf(1,j).EQ.9) THEN
            U(1,J)=0.0
            P(1,J)=P(2,J)
         ELSEIF (nf(1,j).EQ.8) THEN
            u(1,j)=UIn*COS(FreqIn*t)
            u(0,j)=u(1,j)
            v(1,j)=2.*VIn*COS(FreqIn*t) - v(2,j)
            P(1,J)=P(2,J)
         ELSEIF (nf(1,j).EQ.7) THEN
            IF(nf(2,j).NE.0) u(1,j)=u(2,j)
            u(0,j)=2.0*u(1,j)-u(2,j)
!c... negative outflow
            IF (u(1,j).GT.0.0) u(0,j)=u(1,j)
            P(1,J)=POutW 
         ENDIF
!C
!c-- eastern boundary
         V(IMaxUs,J)=SgnR*V(IM1Us,J)
         IF (nf(IMaxUs,j).EQ.9) THEN
            U(IM1Us,J)=0.0
            P(IMaxUs,J)=P(IM1Us,J)
         ELSEIF (nf(IMaxUs,j).EQ.8) THEN
            u(IM1Us,j)=UIn*COS(FreqIn*t)
            u(IMaxUs,j)=u(IM1Us,j)
            v(IMaxUs,j)=2.0*VIn*COS(FreqIn*t)-v(IM1Us,j)
            P(IMaxUs,J)=P(IM1Us,J)
         ELSEIF (nf(IMaxUs,j).EQ.7) THEN
            IF (nf(IM1Us,j).NE.0) u(IM1Us,j)=u(IM2Us,j)
            u(IMaxUs,j)=2.0*u(IM1Us,j)-u(IM2Us,j)
!c... negative outflow
            IF (u(IM1Us,j).LT.0.0) u(IMaxUs,j)=u(IM1Us,j)
            P(IMaxUS,J)=POutE 
!C            P(IMaxUS,J)=P(IM1Us,J)
         ENDIF
!C
 100  CONTINUE
!C
      DO 200 I=1,IMaxUs
!c-- southern boundary
         U(I,1)=SgnB*U(I,2)
         IF (nf(i,1).EQ.9) THEN
            V(I,1)=0.0
            P(I,1)=P(I,2)
         ELSEIF (nf(i,1).EQ.8) THEN
            U(I,1)=2.0*UIn*COS(FreqIn*t)-U(I,2)
            v(i,1)=VIn*COS(FreqIn*t)
            v(i,0)=v(i,1)
            P(I,1)=P(I,2)
         ELSEIF (nf(i,1).EQ.7) THEN
            IF (nf(i,2).NE.0) v(i,1)=v(i,2)
            v(i,0)=2.0*v(i,1)-v(i,2)
!c... negative outflow
            IF (v(i,1).GT.0.0) v(i,0)=v(i,1)
            P(I,1)=POutS
         ENDIF
!C
!c-- northern boundary
         U(I,JMaxUs)=SgnT*U(I,JM1Us)
         IF (nf(i,JMaxUs).EQ.9) THEN
            V(I,JM1Us)=0.0
            P(I,JMaxUs)=P(I,JM1Us)
         ELSEIF (nf(i,JMaxUs).EQ.8) THEN
            U(I,JMaxUs)=2.0*UIn*COS(FreqIn*t)-U(I,JM1Us)
            v(i,JM1Us)=VIn*COS(FreqIn*t)
            v(i,JMaxUs)=v(i,JM1Us)
            P(I,JMaxUs)=P(I,JM1Us)
         ELSEIF (nf(i,JMaxUs).EQ.7) THEN
            IF (nf(i,JM1Us).NE.0) v(i,JM1Us)=v(i,JM2Us)
            v(i,JMaxUs)=2.0*v(i,JM1Us)-v(i,JM2Us)
!c... negative outflow
            IF (v(i,JM1Us).LT.0.0) v(i,JMaxUs)=v(i,JM1Us)
            P(I,JMaxUs)=POutW
         ENDIF
!C
  200 CONTINUE

      RETURN
      END