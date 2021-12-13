SUBROUTINE PRT(N)
use gridar
use phys
use apert
use bndcon
use orga
use orga2
use times
use prints
use tank
use fluid
use condif
      IMPLICIT NONE

!      INTEGER IMAX,JMAX,NSMAX,SDTMAX,mmpnts,NRFMAX,NRHMAX
!      PARAMETER (IMAX=130, JMAX=130, NSMAX=20, SDTMAX=10, 
!     >           mmpnts=25, NRFMAX=10, NRHMAX=10)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX,RDY,
!     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
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

!      COMMON /ORGA/ NF
!      INTEGER NF(IMAX,JMAX)

!      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
!     >               StrtP
!      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
!      REAL  Epsi,OmStrt,Alpha,NrmRhs

!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

!      COMMON /PRINTS/ TPrt,PrtDT,TAdPrt,TCnf,NP
!     >   ,PGnu,PAvs,Puvpf,PVelop,PForce
!      REAL TPrt, PrtDT, TAdPrt, TCnf(2001)
!      INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce

!      COMMON /TANK/ XMax, XMin, YMax, YMin
!      REAL  XMax, XMin, YMax, YMin

!      COMMON /FLUID/ Nu
!      REAL  Nu

!      COMMON /CONDIF/ RU,RV,divuv,divmax
!      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
!      iNTEGER divuv

!c local variables
      INTEGER i,j,n,imid,jstap
      REAL dx,dy

!c----------------------------------------------------------------------
    1 FORMAT (/)
    2 FORMAT ('XMin=',F7.3,3X,'XMax=',F7.3,3X,'YMin=',F7.3,3X,&
         'YMax=',F7.3)    
    3 FORMAT (I3,' + 2 X-INTERVALS , MEAN MESHSIZE =',F8.4)
    4 FORMAT (I3,' + 2 Y-INTERVALS , MEAN MESHSIZE =',F8.4)
    5 FORMAT ('TIME STEP DelT           =',E12.4)
    7 FORMAT ('BOUNDARY CONDITIONS (1=SLIP, 2=NO-SLIP, 7=OUTFLOW, ',&
      '8=INFLOW)',/,&
      '  RIGHT',i3,3X,'LEFT',i3,3X,'TOP',i3,3X,'BOTTOM',i3)
   71 FORMAT ('INFLOW AMPLITUDE UIn =',E11.4, ' VIn =',E11.4 )
    8 FORMAT ('KINEMATIC VISCOSITY Nu   =',E10.3)
   12 FORMAT ('PRESSURE CONVERGENCE Epsi=',E9.2)
   13 FORMAT ('UPWIND PARAMETER Alpha   =',F4.1)
   15 FORMAT (8F10.4)
   22 FORMAT ('T=',E12.5,2X,'div=',E12.5,2X,I4,3X,14(F9.5))
   41 FORMAT (1H ,48I2)
   42 FORMAT (1H ,17F9.4)
   43 FORMAT (1H ,20E12.4)

!C--------------------------------------------------------------------
      IF (N .EQ. 1) THEN
!C
!C--- initializing print actions
!C
         WRITE(6,1)
         WRITE(6,2) XMin,XMax,YMin,YMax
         WRITE(6,*)
         DX=(XMax-XMin)/IM2Us
         WRITE(6,3) IM2Us,DX
         WRITE(6,15) (X(I), I=1,IM1Us)
         WRITE(6,*)
         DY=(YMax-YMin)/JM2Us
         WRITE(6,4) JM2Us, DY
         WRITE(6,15) (Y(J), J=1,JM1Us)
         WRITE(6,*)
         WRITE(6,5) DelT
         WRITE(6,7) WR, WL, WT, WB
         WRITE(6,71) UIn,VIn
         WRITE(6,8) Nu
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,12) Epsi
         WRITE(6,13) Alpha
         WRITE(6,*)

!C.... print initial configuration
         CALL PRTCNF(1)
         WRITE(6,1)

!C--------------------------------------------------------------------
!C        regular printout every PrtDt time steps
!C--------------------------------------------------------------------
      ELSEIF (N .EQ. 2) THEN
!C 
!c** info for screen and 'cfd.out'
!c--- horizontal velocity along x-mid plane is chosen here
!c
         jstap=MAX0(JM2Us/6,1)
         IMID=IMaxUs/2+1
!c         WRITE(6,22) T,divmax,Iter,(U(IMID,J),J=2,JM1Us,jstap)
         WRITE(6,22) T,divmax,Iter,U(IMID,2),U(IMID,JM2Us/2),&
                    U(IMID,JM2Us)
        
         IF (Puvpf.EQ.1) THEN
!c... increase counter for number of calls of postpro-routines
        	 NP = NP + 1
!C--- velocities, forces, fluxes and streamlines
!c--  output of solution in all grid points
	         CALL PRTCNF(2)
         ENDIF   

!c-- compute next small print time
         TPrt = TPrt + PrtDT
!C
!C--------------------------------------------------------------------
!C--------------------------------------------------------------------
!C
      ELSEIF (N .EQ. 3) THEN
!C
!C*****  afsluitende uitvoer
!C
!c-- averages for time step and iterations
         WRITE(6,*)
         WRITE(6,'(a21,i6)') 'No. of time steps  = ', t_cycle
         WRITE(6,'(a21,e9.2)') 'average time step  = ',&
                                           (TFin-TStart)/t_cycle 
         WRITE(6,'(a21,f6.1)') 'average iterations = ',&
                                 FLOAT(ItSum)/FLOAT(t_cycle)
         WRITE(6,*)         

!c... increase counter for number of calls of postpro-routines
         NP = NP + 1 
!c--- closing of files
!C--- velocities, forces, fluxes and streamfunctions
         CALL PRTCNF(2)

      ENDIF

      RETURN
      END
!c  End of subroutine PRT
!c
!cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c
SUBROUTINE PRTCNF(N)
use gridar
use phys
use apert
use bndcon
use orga
use orga2
use times
use prints
use tank
use fluid
use condif
use workar
      IMPLICIT NONE
      INTEGER N
!c
!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX,RDY,
!     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

!      COMMON /TANK/ XMax, XMin, YMax, YMin
!      REAL  XMax, XMin, YMax, YMin
!
!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

!      COMMON /PRINTS/ TPrt,PrtDT,TAdPrt,TCnf,NP
!     >   ,PGnu,PAvs,Puvpf,PVelop,PForce
!      REAL TPrt, PrtDT, TAdPrt, TCnf(2001)
!      INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce

!      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
!      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
!     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
!     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

!      COMMON /DXRL/ DXRI,DXLI,DYTJ,DYBJ
!      REAL DXRI(IMAX),DXLI(IMAX),DYTJ(JMAX),DYBJ(JMAX)

!      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
!     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
!      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
!     1     , POutN, POutE, POutS, POutW
!      INTEGER WB, WL, WR, WT

!      COMMON /ORGA/ NF
!      INTEGER NF(IMAX,JMAX)

!      COMMON /WORKAR/ PSI
!      REAL PSI(IMAX,JMAX)

!c  local variables
      REAL UC(IMAX,JMAX),VC(IMAX,JMAX)
      INTEGER I,J
      CHARACTER*4 nummer
!c
      IF (N.EQ.1) THEN

         OPEN (UNIT=28, FILE='config.dat')
         WRITE(28,'(I3,4F12.5)') 0, XMin, XMax, YMin, YMax
         WRITE(28,'(3I5)') IMaxUs, JMaxUs
         WRITE(28,'(10F12.5)') (XI(I),I=1,IMaxUs)
         WRITE(28,'(10F12.5)') (YJ(J),J=1,JMaxUs)
         DO J=1,JMaxUs
            WRITE(28,'(32I3)') (NF(I,J),I=1,IMaxUs)
         ENDDO         
         WRITE(28,'(6F15.5)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
         CLOSE(28)

      ELSEIF (N.EQ.2) THEN

         IF (NP .GT. 2001) THEN
            WRITE(6,*) 'too many uvpf-files created'
            RETURN
         ENDIF

         WRITE (nummer,'(I4)') NP
         DO i=1,4
            IF (nummer(i:i).EQ.' ') nummer(i:i)='0'
         ENDDO

!c--- average velocities to cell centers

         DO J=1,JMaxUs
         DO I=1,IMaxUs
            UC(I,J)=0.0
            VC(I,J)=0.0
            IF (NF(I,J).NE.9) THEN
               UC(I,J)=0.5*(U(I-1,J)+U(I,J))
               VC(I,J)=0.5*(V(I,J-1)+V(I,J))
               vmax=MAX(SQRT(UC(I,J)**2+VC(I,J)**2),vmax)
               pmin=MIN(p(i,j),pmin)
               pmax=MAX(p(i,j),pmax)
            ENDIF
            PSI(I,J)=0.0
         ENDDO
         ENDDO

!c-- slip boundary velocities
         DO I=2,IM1Us
            J=1
            IF (NF(I,J).EQ.9 .AND. WB.EQ.1) UC(I,J)=UC(I,J+1)
            J=JMaxUs
            IF (NF(I,J).EQ.9 .AND. WT.EQ.1) UC(I,J)=UC(I,J-1)
         ENDDO

         DO J=2,JM1Us
            I=1
            IF (NF(I,J).EQ.9 .AND. WL.EQ.1) VC(I,J)=VC(I+1,J)
            I=IMaxUs
            IF (NF(I,J).EQ.9 .AND. WR.EQ.1) VC(I,J)=VC(I-1,J)
         ENDDO

	 PSI(1,1)=0.0
	 DO I=2,IM1Us
	   PSI(I,1)=PSI(I-1,1)-DelX(I)*V(I,1)
	 ENDDO

	 DO I=1,IM1Us
	 DO J=2,JM1Us
	   PSI(I,J)=PSI(I,J-1)+DelY(J)*U(I,J)
	 ENDDO
	 ENDDO

         OPEN (UNIT=28, FILE='uvpf'//nummer//'.dat')

         WRITE(28,'(F12.5)') T
         DO J=1,JMaxUs
         DO I=1,IMaxUs
!c.. small numbers set to 0 to avoid format problems on Cray
            IF (ABS(UC(I,J)) .LE. 1.E-20) UC(I,J)=0.0
            IF (ABS(VC(I,J)) .LE. 1.E-20) VC(I,J)=0.0
            IF (ABS(P(I,J)) .LE. 1.E-20) P(I,J)=0.0
            IF (ABS(PSI(I,J)) .LE. 1.E-20) PSI(I,J)=0.0

            WRITE(28,'(3E12.4, F7.3, E12.4)') &
              UC(I,J),VC(I,J),P(I,J),FS(I,J),PSI(I,J)
         ENDDO
         ENDDO

         CLOSE(28)

         WRITE(6,'(A13,A4,A21)') 'plotfile UVPF', nummer, & 
                                      '.DAT has been created '

      ENDIF

      RETURN
      END



      SUBROUTINE PRNT(FI,TITLE,numfile,NNI,NNJ)
	integer,parameter::ni=130,nj=130
	real(8),intent(in)::FI(ni,nj)    
	CHARACTER*6,intent(in):: TITLE
	integer,intent(in)::NNI,NNJ,numfile
      WRITE(numfile,20) TITLE 
      IS=-11
  100 IS=IS+12
      IE=IS+11
      IE=MIN(NNI,IE)
      WRITE(numfile,21) (I,I=IS,IE) 
      WRITE(numfile,22) 
      DO J=NNJ,1,-1
        WRITE(numfile,23) J,(FI(I,J),I=IS,IE) 
      END DO
      IF(IE.LT.NNI) GO TO 100
   20 FORMAT(2X,26('*-'),7X,A6,7X,26('-*')) 
   21 FORMAT(3X,'I = ',I3,11I10)
   22 FORMAT(2X,'J')
   23 FORMAT(1X,I3,1P12E10.2) 

      RETURN
      END 


	  subroutine output(numfile)
	  use const
	  use gridar
	  use phys
	  integer numfile
	  character filename  
      !open(numfile,file=filename)
	  call PRNT(U,"U",numfile,IMaxUs,JMaxUs)
	  call PRNT(V,"V",numfile,IMaxUs,JMaxUs)
	  call PRNT(P,"P",numfile,IMaxUs,JMaxUs)
      close(numfile)

	  end subroutine output