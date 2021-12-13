
      SUBROUTINE TILDE
	  use gridar
	  use phys
	  use orga
	  use orga2
	  use fluid
	  use times
	  use condif
      IMPLICIT NONE

!C NEW version -- conservative symmetry preserving (November 1998)
!C Only plane geometry -- axisymmetry has still to be implemented!!!!!!!!

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
!c
!      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX,RDY,
!     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
!      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
!     1   RDX(IMAX),RDY(JMAX),CX,CY
!      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

!      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
!      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
 !    1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
!     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

!      COMMON /ORGA/ NF
!      INTEGER NF(IMAX,JMAX)

!      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
!     >               StrtP
!      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
!      REAL  Epsi,OmStrt,Alpha,NrmRhs

!      COMMON /FLUID/ Nu
!      REAL  Nu

!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

!      COMMON /CONDIF/ RU,RV,divuv,divmax
!      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
!      INTEGER divuv

!c local variables
      REAL dxc,dyc,rdxc,rdxe,rdxw,rdyc,rdyn,rdys
      REAL hue,huw,hvn,hvs,csom,ce,cw,cn,cs,de,dw,dn,ds,klad,temp
      INTEGER i,j

!c      double precision dsecnd,t1,t2
!C
!CALCULATION TIME STEP
temp=1d0/(1d0/DelT+1d0/DelT2)

!C**** CALCULATION OF TEMPORARY VELOCITIES
!C

!C....       TILDE IN X-DIRECTION
!C
      DO 10 J=2,JM1Us
         RDYN=1./(YJ(J+1)-YJ(J))
         RDYS=1./(YJ(J)-YJ(J-1))
      DO 10 I=1,IM1Us

!c--- zeroise all velocities but do not change IN and BND cells
         IF (nf(i,j).LE.OU .AND. nf(i+1,j).LE.OU) u(i,j)=0.0

         IF ( (NF(I,J).EQ.F .OR. nf(i,j).EQ.OU) .AND.&
                (NF(I+1,J).EQ.F .OR. nf(i+1,j).EQ.OU) ) THEN
!c.. both i,j and i+1,j FUL or OUT cell

            DXC=XI(I+1)-XI(I)
            RDXC=1./DXC
            RDXE=RDX(I+1)
            RDXW=RDX(I)

!C
!c.. nieuwe formulering
            HUE=0.25*(UN(I,J)+UN(I+1,J))
            HUW=0.25*(UN(I,J)+UN(I-1,J))
            HVN=0.25*(DelX(I)*VN(I,J)+DelX(I+1)*VN(I+1,J))*RDXC
            HVS=0.25*(DelX(I)*VN(I,J-1)+DelX(I+1)*VN(I+1,J-1))*RDXC

!c coefficients to be scaled by DXC*DELY(J) for (skew)symmetry

            DE=Nu*RDXE+ALPHA*ABS(HUE)
            CE=( HUE - DE) * RDXC
            DW=Nu*RDXW+ALPHA*ABS(HUW)
            CW=(-HUW - DW) * RDXC
            DN=Nu*RDYN+ALPHA*ABS(HVN)
            CN=( HVN - DN) * RDY(J)
            DS=Nu*RDYS+ALPHA*ABS(HVS)
            CS=(-HVS - DS) * RDY(J)
            CSom=-(DE+DW)*RDXC-(DN+DS)*RDY(J)

!c expliciet
           klad=-CE*UN(I+1,J)-CW*UN(I-1,J)-&
                 CN*UN(I,J+1)-CS*UN(I,J-1)+Csom*UN(I,J)
           RU(I,J)=KLAD
           U(I,J)=UN(I,J)+temp*KLAD-temp*(P(i+1,j)-P(i,j))*RDXC

         ENDIF
 10   CONTINUE

!C
!C....       TILDE IN Y-DIRECTION
!C
      DO 20 J=1,JM1Us
         DYC=YJ(J+1)-YJ(J)
         RDYC=1./DYC
         RDYN=RDY(J+1)
         RDYS=RDY(J)
      DO 20 I=2,IM1Us

!c--- zeroise all velocities but do not change IN and BND cells
         IF (nf(i,j).LE.OU .AND. nf(i,j+1).LE.OU) v(i,j)=0.0

         IF ( (NF(I,J).EQ.F .OR. nf(i,j).EQ.OU) .AND.&
                (NF(I,J+1).EQ.F .OR. nf(i,j+1).EQ.OU) ) THEN
!c.. both i,j and i,j+1 FUL or OUT cell

            RDXE=1./(XI(I+1)-XI(I))
            RDXW=1./(XI(I)-XI(I-1))

            HUE=0.25*(DelY(J)*UN(I,J)+DelY(J+1)*UN(I,J+1))*RDYC
            HUW=0.25*(DelY(J)*UN(I-1,J)+DelY(J+1)*UN(I-1,J+1))*RDYC
            HVN=0.25*(VN(I,J)+VN(I,J+1))
            HVS=0.25*(VN(I,J)+VN(I,J-1))

!c coefficients to be scaled by DELX(I)*DYC for (skew)symmetry
            DE=Nu*RDXE+ALPHA*ABS(HUE)
            CE=( HUE - DE) * RDX(I)
            DW=Nu*RDXW+ALPHA*ABS(HUW)
            CW=(-HUW - DW) * RDX(I)
            DN=Nu*RDYN+ALPHA*ABS(HVN)
            CN=( HVN - DN) * RDYC
            DS=Nu*RDYS+ALPHA*ABS(HVS)
            CS=(-HVS - DS) * RDYC
            CSom=-(DE+DW)*RDX(I)-(DN+DS)*RDYC

!c expliciet
            KLAD=-CE*VN(I+1,J)-CW*VN(I-1,J)-&
                 CN*VN(I,J+1)-CS*VN(I,J-1)+Csom*VN(I,J)
            RV(I,J)=KLAD
            V(I,J)=VN(I,J)+temp*KLAD-temp*(P(i,j+1)-P(i,j))*RDYC


         ENDIF
 20   CONTINUE

      RETURN
      END
