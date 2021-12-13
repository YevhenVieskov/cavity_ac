SUBROUTINE GRID
use gridar
use dxrl
use coef2
use const
use tank
use rek

      IMPLICIT NONE

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX,RDY,
!     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

!      COMMON /DXRL/ DXRI,DXLI,DYTJ,DYBJ
!      REAL DXRI(IMAX),DXLI(IMAX),DYTJ(JMAX),DYBJ(JMAX)

!      COMMON /coef2/ CWI,CEI,CNJ,CSJ
!      REAL CWI(IMAX),CEI(IMAX),CNJ(JMAX),CSJ(JMAX)

!      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI
!      REAL EMF,EMF1,EM6,EM10,EP10,PI

!      COMMON /TANK/ XMax, XMin, YMax, YMin
!      REAL  XMax, XMin, YMax, YMin

!      COMMON /REK/ xpos,ypos
!      REAL xpos,ypos

!c  local variables
      REAL stepx,stepy,scalex,scaley,ymid,xmid,ccx,ccy&
          ,reli,kladl,kladr,AbsCX,AbsCY
      INTEGER i,j,jmid,imid
!C
!C**********     GRID IN X-DIRECTION     ***************************
!C
!C... definition of stretched grid
!C    when cx = 0 uniform grid
!C         cx > 0 refinement near position x=Xpos
!C         cx < 0 refinement near both endpoints

         xpos=MAX(xpos,xmin)
         xpos=MIN(xpos,xmax)
         AbsCX=ABS(CX)

      IF (AbsCX.GE.EM6) THEN
!c.. determine number of gridpoints per interval (logarithmic split)
         reli=LOG(AbsCX*(xpos-xmin)+1.)/&
          (LOG(AbsCX*(xpos-xmin)+1.)+LOG(AbsCX*(xmax-xpos)+1.))
         IMID=NINT(reli*float(IM2Us))+1
         XMID=xpos
!c.. stretching factors for both intervals
         kladl=1.0+2.0*(reli-0.5)/(0.5*AbsCX**2+1.0)
         kladr=1.0+2.0*(0.5-reli)/(0.5*AbsCX**2+1.0)
      ENDIF

      IF (AbsCX.LT.EM6) THEN
         STEPX=(XMax-XMin)/FLOAT(IM2Us)
         DO 2 I=1,IM1Us
       2 X(I)=XMin+(I-1)*STEPX

      ELSEIF (CX.GT.0.0) THEN

!c.. refinement near X=Xpos
!c.. left interval
         CCX=AbsCX*kladl
         SCALEX=1.0
         IF (CCX.GT.em6) SCALEX=(XMID-XMin)/TANH(CCX)
         DO 4 I=1,IMID-1
       4 X(I)=XMin+SCALEX*TANH(CCX*FLOAT(I-1)/FLOAT(IMID-1))
!c.. `mid'point
         X(IMID)=XMID
!c.. right interval
         CCX=AbsCX*kladr
         SCALEX=1.0
         IF (CCX.GT.em6) SCALEX=(XMax-XMID)/TANH(CCX)
         DO 5 I=IMID+1,IM1Us
       5 X(I)=XMax-SCALEX*TANH(CCX*FLOAT(IM1Us-I)/FLOAT(IM1Us-IMID))

      ELSEIF (CX.LT.0.0) THEN

!c.. refinement near both endpoints, away from X=Xpos
!c         IMID=INT(0.5*float(IM2Us))+1
!c         XMID=(xmax+xmin)*0.5
!c.. left interval
         CCX=AbsCX
         SCALEX=(XMID-XMin)/TANH(CCX)
         DO 7 I=1,IMID-1
       7 X(I)=XMID-SCALEX*TANH(CCX*FLOAT(IMID-I)/FLOAT(IMID-1))
!c.. midpoint
         X(IMID)=XMID
!c.. right interval
         CCX=AbsCX
         SCALEX=(XMax-XMID)/TANH(CCX)
         DO 8 I=IMID+1,IM1Us
       8 X(I)=XMID+SCALEX*TANH(CCX*FLOAT(I-IMID)/FLOAT(IM1Us-IMID))

      ENDIF
!C
!C... auxiliary quantities
!C
      DO 10 I=2,IM1Us
         DelX(I)=X(I)-X(I-1)
         RDX(I)=1./DelX(I)
         XI(I)=0.5*(X(I-1)+X(I))
   10 CONTINUE
!C
      DelX(1)=DelX(2)
      RDX(1)=1.0/DelX(1)
      XI(1)=XI(2)-DelX(2)
      DelX(IMaxUs)=DelX(IM1Us)
      RDX(IMaxUs)=1.0/DelX(IMaxUs)
      XI(IMaxUs)=XI(IM1Us)+DelX(IM1Us)
      X(IMaxUs)=X(IM1Us)+DelX(IM1Us)
!C

      DO 11 I=2,IM1Us
         DXLI(I)=0.5*(DelX(I)+DelX(I-1))
         DXRI(I)=0.5*(DelX(I)+DelX(I+1))
         CWI(I)=RDX(I)/DXLI(I)
         CEI(I)=RDX(I)/DXRI(I)
  11 CONTINUE
      DXLI(1)=DXLI(2)
      DXRI(1)=DXRI(2)
!C
!C***********     GRID IN Y-DIRECTION      ***************************
!C
!C... definition of stretched grid
!C    when cy = 0 uniform grid
!C         cy > 0 refinement near Ypos
!C         cy < 0 refinement away from Ypos
!C
         ypos=MAX(ypos,ymin)
         ypos=MIN(ypos,ymax)
         AbsCY=ABS(CY)

      IF (AbsCY .GE. EM6) THEN
!c.. determine number of gridpoints per interval (logarithmic split)
         reli=LOG(AbsCY*(ypos-ymin)+1.)/&
          (LOG(AbsCY*(ypos-ymin)+1.)+LOG(AbsCY*(ymax-ypos)+1.))
         JMID=NINT(reli*float(JM2Us))+1
         YMID=ypos
!c.. stretching factors for both intervals
         kladl=1.0+2.0*(reli-0.5)/(0.5*AbsCY**2+1.0)
         kladr=1.0+2.0*(0.5-reli)/(0.5*AbsCY**2+1.0)
      ENDIF
!C
      IF (AbsCY.LT.EM6) THEN

         STEPY=(YMax-YMin)/FLOAT(JM2Us)
         DO 12 J=1,JM1Us
     12  Y(J)=YMin + (J-1)*STEPY

      ELSEIF (CY.GT.0.0) THEN

!c.. refinement near Y=Ypos
!c.. lower interval
         CCY=AbsCY*kladl
         SCALEY=1.0
         IF (CCY.GT.em6) SCALEY=(YMID-YMin)/TANH(CCY)
         DO 14 J=1,JMID-1
      14 Y(J)=YMin+SCALEY*TANH(CCY*FLOAT(J-1)/FLOAT(JMID-1))
!c.. `mid'point
         Y(JMID)=YMID
!c.. upper interval
         CCY=AbsCY*kladr
         SCALEY=1.0
         IF (CCY.GT.em6) SCALEY=(YMax-YMID)/TANH(CCY)
         DO 15 J=JMID+1,JM1Us
      15 Y(J)=YMax-SCALEY*TANH(CCY*FLOAT(JM1Us-J)/FLOAT(JM1Us-JMID))

      ELSEIF (CY.LT.0.0) THEN

!c.. refinement near both endpoints, away from Y=Ypos
!c         JMID=INT(0.5*float(JM2Us))+1
!c         YMID=(ymax+ymin)*0.5
!c.. lower interval
         SCALEY=(YMID-YMin)/TANH(CY)
         DO 17 J=1,JMID-1
      17 Y(J)=YMID-SCALEY*TANH(CY*FLOAT(JMID-J)/FLOAT(JMID-1))
!c.. midpoint
         Y(JMID)=YMID
!c.. upper interval
         SCALEY=(YMax-YMID)/TANH(CY)
         DO 18 J=JMID+1,JM1Us
     18  Y(J)=YMID+SCALEY*TANH(CY*FLOAT(J-JMID)/FLOAT(JM1Us-JMID))

      ENDIF

!C
!C... auxiliary quantities
!C
      DO 20 J=2,JM1Us
         DelY(J)=Y(J)-Y(J-1)
         RDY(J)=1.0/DelY(J)
         YJ(J)=0.5*(Y(J-1)+Y(J))
   20 CONTINUE
!C
      DelY(1)=DelY(2)
      RDY(1)=RDY(2)
      YJ(1)=YJ(2)-DelY(2)
      DelY(JMaxUs)=DelY(JM1Us)
      RDY(JMaxUs)=1.0/DelY(JMaxUs)
      YJ(JMaxUs)=YJ(JM1Us)+DelY(JM1Us)
      Y(JMaxUs)=Y(JM1Us)+DelY(JM1Us)
!C
      DO 25 J=2,JM1Us
         DYBJ(J)=0.5*(DelY(J)+DelY(J-1))
         DYTJ(J)=0.5*(DelY(J)+DelY(J+1))
         CSJ(J)=RDY(J)/DYBJ(J)
         CNJ(J)=RDY(J)/DYTJ(J)
   25 CONTINUE
      DYBJ(1)=DYBJ(2)
      DYTJ(1)=DYTJ(2)
!C

!MAIN CODE
     ! DX=DelX(3)-DelX(2)
      !DY=DelY(3)-DelY(2)
      RETURN
      END