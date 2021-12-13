SUBROUTINE MILU
use gridar
use phys
use coefp
use orga2
use times
      IMPLICIT NONE

!      INTEGER IMAX,JMAX
!      PARAMETER (IMAX=130, JMAX=130)
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

!      COMMON /COEFP/ CN, CS, CE, CW, CC, DIV
!      REAL CN(IMAX,JMAX), CS(IMAX,JMAX), CE(IMAX,JMAX),
!     1     CW(IMAX,JMAX), CC(IMAX,JMAX), DIV(IMAX,JMAX)

!      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
!     >               StrtP
!      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
!      REAL  Epsi,OmStrt,Alpha,NrmRhs

!      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
!      REAL  T, DelT, DelTMx, TFin, TStart
!      INTEGER Cycle

      REAL rkr,rkrold,zaz,maxkmr,maxp,somw,somp,&
          alf,beta,delta,cee,cnn,css,cww

      REAL r(IMAX,JMAX),z(IMAX,JMAX),az(IMAX,JMAX),&
          w(IMAX,JMAX),d1(IMAX,JMAX),dsi(IMAX,JMAX)
      INTEGER i,j,M,N
      SAVE d1,dsi

!C ##################        INIT         ###################

      M=IMaxUs
      N=JMaxUs
      DO 5 I=1,IMaxUs
        r(i,1)=0.0
        r(i,N)=0.0
        z(i,1)=0.0
        z(i,N)=0.0
        az(i,1)=0.0
        az(i,N)=0.0
        w(i,1)=0.0
        w(i,N)=0.0
        dsi(i,1)=0.0
        dsi(i,N)=0.0
    5 CONTINUE
      DO 10 J=1,JMaxUs
        r(1,j)=0.0
        r(M,j)=0.0
        z(1,j)=0.0
        z(M,j)=0.0
        az(1,j)=0.0
        az(M,j)=0.0
        w(1,j)=0.0
        w(M,j)=0.0
        dsi(1,j)=0.0
        dsi(M,j)=0.0
   10 CONTINUE

      IF (StrtP.eq.0) then
        DO 15 j=1,JMaxUs
        DO 15 i=1,IMaxUs
	   p(i,j)=0.0
   15   CONTINUE
      ENDIF

!C ##################   Preconditioning   ###################
!C Onderstaande loop kan veel efficienter
!C Ook verderop is winst te behalen. Vectorisatie is mogelijk door
!C langs diagonalen (zuidoost-noordwest) te nummeren

      IF (T_CYCLE.EQ.1) THEN

         beta=(1.0+5./(IM1Us*JM1Us))
         DO 20 j=2,JM1Us
         DO 20 i=2,IM1Us 
           cee=0.0
           cww=0.0
           cnn=0.0
           css=0.0
           IF (i.GT.2) THEN
              cee=ce(i-1,j)
              IF (j.LT.JM1Us) css=cs(i-1,j+1)
           ENDIF
           IF (j.GT.2) THEN
             cnn=cn(i,j-1)
             IF (i.LT.IM1Us) cww=cw(i+1,j-1)
           ENDIF
           Dsi(i,j)=1./(beta*Cc(i,j)&
                           -Cw(i,j)*(cee+css)*Dsi(i-1,j)&
                           -Cs(i,j)*(cnn+cww)*Dsi(i,j-1))
           d1(i,j)=Cc(i,j)-2./Dsi(i,j)
   20    CONTINUE

      ENDIF

!C ##################  Conjugate gradient  ###################

      Iter = 0
      rkr=0.0
      DO 40 j=2,JM1Us
      DO 40 i=2,IM1Us
        r(i,j)=div(i,j)-(Cc(i,j)*p(i,j)+Cs(i,j)*p(i,j-1)+&
              Cw(i,j)*p(i-1,j)+Cn(i,j)*p(i,j+1)+Ce(i,j)*p(i+1,j))
        r(i,j)=(r(i,j)-Cs(i,j)*r(i,j-1)-Cw(i,j)*r(i-1,j))*Dsi(i,j)
        z(i,j)=r(i,j)/Dsi(i,j)
        rkr=rkr+r(i,j)*z(i,j)
   40 CONTINUE

      somw=0.0
      somp=0.0
      DO j=2,JM1Us
      DO i=2,IM1Us
         somw=somw+z(i,j)*z(i,j)
         somp=somp+p(i,j)*p(i,j)
      ENDDO
      ENDDO
      maxkmr=SQRT(somw)
      maxp=SQRT(somp)
      delta=maxkmr/(maxp+1.0e-10)
      IF (delta.LT.1e-6*Epsi) GOTO 2000

      DO 45 j=2,JM1Us
      DO 45 i=2,IM1Us
        p(i,j)=p(i,j)/Dsi(i,j)+Cn(i,j)*p(i,j+1)+Ce(i,j)*p(i+1,j)
 45   CONTINUE

 1000 Iter=Iter+1
      rkrold=rkr

      DO 50 j=JM1Us,2,-1
      DO 50 i=IM1Us,2,-1
        az(i,j)=(z(i,j)-Cn(i,j)*az(i,j+1)-Ce(i,j)*az(i+1,j))*Dsi(i,j)
        w(i,j)=z(i,j)+d1(i,j)*az(i,j)
 50   CONTINUE

      zaz=0.0

      DO 60 j=2,JM1Us
      DO 60 i=2,IM1Us
        w(i,j)=(w(i,j)-Cs(i,j)*w(i,j-1)-Cw(i,j)*w(i-1,j))*Dsi(i,j)
        az(i,j)=az(i,j)+w(i,j)
        zaz=zaz+az(i,j)*z(i,j)
 60   CONTINUE

      alf=rkr/zaz
      rkr=0.0
      DO 70 j=2,JM1Us
      DO 70 i=2,IM1Us
        p(i,j)=p(i,j)+alf*z(i,j)
        r(i,j)=r(i,j)-alf*az(i,j)
        w(i,j)=r(i,j)/Dsi(i,j)
        rkr=rkr+r(i,j)*w(i,j)
   70 CONTINUE

      somw=0.0
      somp=0.0
      DO j=2,JM1Us
      DO i=2,IM1Us
         somw=somw+w(i,j)*w(i,j)
         somp=somp+p(i,j)*p(i,j)
      ENDDO
      ENDDO
      maxkmr=SQRT(somw)
      maxp=SQRT(somp)
      delta=maxkmr/(maxp+1.0e-10)

      beta=rkr/rkrold
      DO 80 j=2,JM1Us
      DO 80 i=2,IM1Us
        z(i,j)=w(i,j)+beta*z(i,j)
   80 CONTINUE

      IF (delta.GT.Epsi .AND. Iter.LT.ItMax) GOTO 1000

      DO 90 j=JM1Us,2,-1
      DO 90 i=IM1Us,2,-1
        p(i,j)=(p(i,j)-Cn(i,j)*p(i,j+1)-Ce(i,j)*p(i+1,j))*Dsi(i,j)
   90 CONTINUE

 2000 ItSum = ItSum + Iter

      return
      end
