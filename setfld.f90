SUBROUTINE SETFLD
use gridar
use phys
use orga
use bndcon
use tank
use apert  
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

!      COMMON /ORGA/ NF
!      INTEGER NF(IMAX,JMAX)

!      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
!     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
!      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
!     1     , POutN, POutE, POutS, POutW
!      INTEGER WB, WL, WR, WT

!      COMMON /TANK/ XMax, XMin, YMax, YMin
!      REAL  XMax, XMin, YMax, YMin

!c local variables
      INTEGER i,j
!C
!C**** PRESSURE and VELOCITIES
!C
      DO 10 J=1,JMaxUs
      DO 10 I=1,IMaxUs
         PN(I,J)=0.0
         P(I,J)=0.0
         U(I,J)=0.0
         V(I,J)=0.0
         UN(I,J)=U(I,J)
         VN(I,J)=V(I,J)
         FS(I,J)=1.0
         NF(I,J)=0
 10   CONTINUE

      DO J=1,JMaxUs
         U(0,J)=U(1,J)
         UN(0,J)=U(0,J)
      ENDDO
      DO I=1,IMaxUs
         V(I,0)=V(I,1)
         VN(I,0)=V(I,0)
      ENDDO
!C
!c... define boundary conditions for F
!c... solid: NF=9, F=0;  inflow: NF=8, F=1;  outflow: NF=7, F=1
!c
      DO 70 I=2,IM1Us
         NF(I,1)=B
         NF(I,JMaxUs)=B
   70 CONTINUE
!C
       DO 71 J=1,JMaxUs
         NF(1,J)=B
         NF(IMaxUs,J)=B
  71   CONTINUE

!c...  INFLOW / OUTFLOW boundaries
!c       
      IF (wb.EQ.OU) THEN
         DO i=1,IMaxUs
            nf(i,1)=wb
         ENDDO
      ELSEIF (wb.EQ.IN) THEN
         DO i=2,IM1Us
            nf(i,1)=wb
         ENDDO
      ENDIF

      IF (wt.EQ.OU) THEN
         DO i=1,IMaxUs
            nf(i,JMaxUs)=wt
         ENDDO
      ELSEIF (wt.EQ.IN) THEN
         DO i=2,IM1Us
            nf(i,JMaxUs)=wt
         ENDDO
      ENDIF

      IF (wl.EQ.OU) THEN
         DO j=1,JMaxUs
            nf(1,j)=wl
         ENDDO
      ELSEIF (wl.EQ.IN) THEN
         DO j=2,JM1Us
            nf(1,j)=wl
         ENDDO
      ENDIF

      IF (wr.EQ.OU) THEN
         DO j=1,JMaxUs
            nf(IMaxUs,j)=wr
         ENDDO
      ELSEIF (wr.EQ.IN) THEN
         DO j=2,JM1Us
            nf(IMaxUs,j)=wr
         ENDDO
      ENDIF
!c
!c... adaption of vof-function to cell labels
!c
      DO j=1,JMaxUs
      DO i=1,IMaxUs
         IF (nf(i,j).EQ.B) THEN
            fs(i,j)=0.0
         ELSEIF (nf(i,j).EQ.OU .OR. nf(i,j).EQ.IN) THEN
            fs(i,j)=1.0
         ENDIF

      ENDDO
      ENDDO

!c      call bcbnd
!C
      RETURN
      END