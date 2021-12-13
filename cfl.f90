SUBROUTINE CFLCHK2(FlgAdj)
use const
use gridar
use phys
use times
use adjust
      IMPLICIT NONE

      INTEGER FlgAdj

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
!     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
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

!      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
!      REAL CFLMin, CFLMax
!      INTEGER NDTAd, MxDTAd, SlwFlw

!c-- local variables
      INTEGER i,j
      REAL kladx,klady,CFL2,CFLx,CFLy,CFLtmp

      FlgAdj = 0
      CFL2=0.0
      DO 20 J=2,JM1Us
         klady=delt2/dely(j)
      DO 20 I=2,IM1Us
         kladx=delt2/delx(i)
         CFLx=ABS(u(i,j))*kladx
         CFLy=ABS(v(i,j))*klady
         CFLtmp=AMAX1(CFLx,CFLy)
         CFL2=AMAX1(CFLtmp,CFL2)
 20   CONTINUE

!c--- check for large time step

      IF (CFL2.GT.5.0) THEN
         WRITE(6,*) 'CFL2 is larger than 5'
         STOP 'CFL2 is larger than 5'
      ENDIF

      IF (CFL2.GT.CFLMax2) THEN
!c--- reduce time step
         FlgAdj=1
         IF (NDTAd2.GE.ABS(MxDTAd2)) THEN
            FlgAdj=0
         ENDIF
      ENDIF

!c--- check for small time step;
!c    if 25 times in a row  CFL < CFLMin  then enlarge time step

      IF (CFL2.LT.CFLMin2) THEN
         SlwFlw2=SlwFlw2+1
      ELSE
         SlwFlw2=0
      ENDIF

      IF (SlwFlw2.EQ.25 .AND. DelT2.LT.0.5*DelTMx2) THEN
         SlwFlw2=0
         FlgAdj=-1
      ENDIF


!c--- adjust time step if required

      IF (FlgAdj .NE. 0) CALL DTADJ2(FlgAdj)

      RETURN
      END


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE DTADJ2(FlgAdj)
use const
use gridar
use phys
use times
use adjust

      IMPLICIT NONE

      INTEGER FlgAdj

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

!      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
!      REAL CFLMin, CFLMax
!      INTEGER NDTAd, MxDTAd, SlwFlw

!c local variables
      INTEGER i,j

!c    when MxDTAd > 0 then DelT can be reduced and enlarged
!c    when MxDTAd < 0 then DelT may only be reduced

!c--- NDTAd is effective number of time-step reductions
!c    FlgAdj = 1 reduction of DelT; FlgAdj = -1 enlargement

      NDTAd2 = NDTAd2 + FlgAdj

!c-- reset time level

      t2 = t2 - delt2

!c-- time step is adjusted

      IF (FlgAdj.EQ.1) THEN
         delt2 = 0.5*delt2
         WRITE(6,4000) t2, delt2
      ELSEIF (FlgAdj.EQ.-1 .AND. MxDTAd2.GT.0) THEN
         delt2 = 2.0*delt2
         WRITE(6,4001) t2, delt2
      ENDIF

!c-- reset all variables

      DO 10 J=1,JMaxUs
         U(0,J)=UN(0,J)
      DO 10 I=1,IMaxUs
         P(I,J)=PN(I,J)
         U(I,J)=UN(I,J)
         V(I,J)=VN(I,J)
 10   CONTINUE
      DO 11 I=1,IMaxUs
         V(I,0)=VN(I,0)
 11   CONTINUE

 4000 FORMAT (/,2X,' T =',E12.5,'  TIME STEP REDUCED: DelT2 =',E11.4,/)
 4001 FORMAT (/,2X,' T =',E12.5,'  TIME STEP DOUBLED: DelT2 =',E11.4,/)

      RETURN
      END

SUBROUTINE CFLCHK(FlgAdj)
use const
use gridar
use phys
use times
use adjust
      IMPLICIT NONE

      INTEGER FlgAdj

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
!     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
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

!      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
!      REAL CFLMin, CFLMax
!      INTEGER NDTAd, MxDTAd, SlwFlw

!c-- local variables
      INTEGER i,j
      REAL kladx,klady,CFL,CFLx,CFLy,CFLtmp

      FlgAdj = 0
      CFL=0.0
      DO 20 J=2,JM1Us
         klady=delt/dely(j)
      DO 20 I=2,IM1Us
         kladx=delt/delx(i)
         CFLx=ABS(u(i,j))*kladx
         CFLy=ABS(v(i,j))*klady
         CFLtmp=AMAX1(CFLx,CFLy)
         CFL=AMAX1(CFLtmp,CFL)
 20   CONTINUE

!c--- check for large time step

      IF (CFL.GT.5.0) THEN
         WRITE(6,*) 'CFL is larger than 5'
         STOP 'CFL is larger than 5'
      ENDIF

      IF (CFL.GT.CFLMax) THEN
!c--- reduce time step
         FlgAdj=1
         IF (NDTAd.GE.ABS(MxDTAd)) THEN
            FlgAdj=0
         ENDIF
      ENDIF

!c--- check for small time step;
!c    if 25 times in a row  CFL < CFLMin  then enlarge time step

      IF (CFL.LT.CFLMin) THEN
         SlwFlw=SlwFlw+1
      ELSE
         SlwFlw=0
      ENDIF

      IF (SlwFlw.EQ.25 .AND. DelT.LT.0.5*DelTMx) THEN
         SlwFlw=0
         FlgAdj=-1
      ENDIF


!c--- adjust time step if required

      IF (FlgAdj .NE. 0) CALL DTADJ(FlgAdj)

      RETURN
      END


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE DTADJ(FlgAdj)
use const
use gridar
use phys
use times
use adjust

      IMPLICIT NONE

      INTEGER FlgAdj

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

!      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
!      REAL CFLMin, CFLMax
!      INTEGER NDTAd, MxDTAd, SlwFlw

!c local variables
      INTEGER i,j

!c    when MxDTAd > 0 then DelT can be reduced and enlarged
!c    when MxDTAd < 0 then DelT may only be reduced

!c--- NDTAd is effective number of time-step reductions
!c    FlgAdj = 1 reduction of DelT; FlgAdj = -1 enlargement

      NDTAd = NDTAd + FlgAdj

!c-- reset time level

      t = t - delt

!c-- time step is adjusted

      IF (FlgAdj.EQ.1) THEN
         delt = 0.5*delt
         WRITE(6,4000) t, delt
      ELSEIF (FlgAdj.EQ.-1 .AND. MxDTAd.GT.0) THEN
         delt = 2.0*delt
         WRITE(6,4001) t, delt
      ENDIF

!c-- reset all variables

      DO 10 J=1,JMaxUs
         U(0,J)=UN(0,J)
      DO 10 I=1,IMaxUs
         P(I,J)=PN(I,J)
         U(I,J)=UN(I,J)
         V(I,J)=VN(I,J)
 10   CONTINUE
      DO 11 I=1,IMaxUs
         V(I,0)=VN(I,0)
 11   CONTINUE

 4000 FORMAT (/,2X,' T =',E12.5,'  TIME STEP REDUCED: DelT =',E11.4,/)
 4001 FORMAT (/,2X,' T =',E12.5,'  TIME STEP DOUBLED: DelT =',E11.4,/)

      RETURN
      END
