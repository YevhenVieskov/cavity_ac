subroutine setpar
use gridar
use phys
use bndcon
use orga2
use condif
use const
use times
use adjust
use tank
use fluid
use prints
use rek
 IMPLICIT NONE

!      INTEGER IMAX,JMAX,NSMAX,SDTMAX,NRFMAX,NRHMAX,mmpnts
!      PARAMETER (IMAX=130, JMAX=130, NSMAX=20, SDTMAX=10, 
!     >           NRFMAX=10, NRHMAX=10, mmpnts=25)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
!     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

!      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
!      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
!     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
!     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

!      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
!     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
!      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
!     1     , POutN, POutE, POutS, POutW
!      INTEGER WB, WL, WR, WT

!      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
!     >               StrtP
!      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
!      REAL  Epsi,OmStrt,Alpha,NrmRhs

!      COMMON /CONDIF/ RU,RV,divuv,divmax
!      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
!      INTEGER divuv

!      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI
!      REAL EMF,EMF1,EM6,EM10,EP10,PI

!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

!      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
!      REAL CFLMin, CFLMax
!      INTEGER NDTAd, MxDTAd, SlwFlw

!      COMMON /PRINTS/ TPrt,PrtDT,TAdPrt,TCnf,NP
!     >   ,PGnu,PAvs,Puvpf,PVelop,PForce
!      REAL TPrt, PrtDT, TAdPrt, TCnf(2001)
!      INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce

!      COMMON /TANK/ XMax, XMin, YMax, YMin
!      REAL  XMax, XMin, YMax, YMin

!      COMMON /FLUID/ Nu
!      REAL  Nu

!      COMMON /REK/ posx,posy
!      REAL posx,posy

! local variables
      REAL dxx,dyy,h1,PIn
      INTEGER i,j,IPIn
      open(5,file='cfd.in')
!C**** CONSTANT PARAMETERS
 !     PI = 3.14159265359
 !     EP10 = 1.0e10
!      EM6  = 1.0e-6
!      EM10 = 1.0e-10
!      EMF  = 1.0e-5
!      EMF1 = 1.0-EMF

!C**** TANK GEOMETRY
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) XMIN,XMAX,YMIN,YMAX
      READ(5,*)

!C**** GRID DEFINITION
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) IMaxUs,JMaxUs,CX,CY,xpos,ypos
      READ(5,*)
      IF (IMaxUs.GT.IMAX .OR. JMaxUs.GT.JMAX) STOP  'gridsize too large'
      IM1Us = IMaxUs -1
      JM1Us = JMaxUs -1
      IM2Us = IMaxUs -2
      JM2Us = JMaxUs -2

      CALL GRID   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!C**** MATERIAL QUANTITIES
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) Nu
      READ(5,*)

!C**** BOUNDARY CONDITIONS AND IN-/OUTFLOW
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) WL, WR, WT, WB, UIn, VIn
      READ(5,*)
      IF ((WL.NE.1).and.(WL.NE.2).and.(WL.NE.7).and.(WL.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition left'
         stop
      ENDIF
      IF ((WR.NE.1).and.(WR.NE.2).and.(WR.NE.7).and.(WR.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition right'
         stop
      ENDIF
      IF ((WT.NE.1).and.(WT.NE.2).and.(WT.NE.7).and.(WT.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition top'
         stop
      ENDIF
      IF ((WB.NE.1).and.(WB.NE.2).and.(WB.NE.7).and.(WB.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition bottom'
         stop
      ENDIF
      FreqIn=0
      IPIn=0
      PIn=-1.0
      POutN=0.0
      POutE=0.0
      POutS=0.0
      POutW=0.0
      IF (IPIn.EQ.1) POutW=PIn
      IF (IPIn.EQ.2) POutE=PIn
      IF (IPIn.EQ.3) POutN=PIn
      IF (IPIn.EQ.4) POutS=PIn
!c
      SgnL=1.0
      IF (WL.EQ.2) SgnL=-1.0
      SgnR=1.0
      IF (WR.EQ.2) SgnR=-1.0
      SgnT=1.0
      IF (wt.EQ.2) SgnT=-1.0
      SgnB=1.0
      IF (wb.EQ.2) SgnB=-1.0

!C**** NUMERICAL MODEL PARAMETERS
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) Alpha, Divuv
      READ(5,*)

      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) Epsi, ItMax, StrtP
      ItSum = 0
      OmStrt=1.0
      IMilu=1  !! VIBOR MILU
      READ(5,*)

!C**** TIME PARAMETERS
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) TFin,DelT,CFLMax
      READ(5,*)
      T=0.0
	  T2=0d0
	  TFin2=50*TFin
      DelT2=DelT
	  CFLMax2=CFLMax
!C--- timestep DelT adjustment to stability limits
!c    look for smallest mesh size
      DXX = XMax-XMin
      DO I=1,IM1Us
         IF (DelX(I).LT.DXX) DXX=DelX(I)
      ENDDO
      DYY = YMax-YMin
      DO J=1,JM1Us
         IF (DelY(J).LT.DYY) DYY=DelY(J)
      ENDDO
      H1 = 1./(DXX*DXX) + 1./(DYY*DYY)
      DelTMx = .5/(H1*(Nu+em10))
      DelT = MIN(0.9*DelTMx, DelT)

!C**** TIME-STEP ADJUSTMENT CONTROL
      CFLMin=0.4*CFLMax
      MxDTAd=5
      NDTAd=0
      SlwFlw=0

	  CFLMin2=0.4*CFLMax2
      MxDTAd2=5
      NDTAd2=0
      SlwFlw2=0

!C**** PRINT/PLOT CONTROL
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) PrtDT,Puvpf
      READ(5,*)

!C**** ORGANISATIONAL PARAMETERS
      NP = 0
      t_cycle = 0
      VMAX=0.0
      pmin=0.0
      pmax=0.0

      RETURN

 9999 STOP 'error in reading CFD.IN'

      

end subroutine setpar