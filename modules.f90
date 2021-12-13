module const
integer,parameter::imax=130,jmax=130
real(8),parameter::ep10=1d-10,em6=1d-6,em10=1d-10,emf=1d-5,emf1=emf-1d0,pi=3.14159265359
integer,parameter::B=9,IN=8,OU=7,F=0,S=1,E=2

integer,parameter::BB=100,BF=101,BS=102,BE=103,BI=104,BO=105
integer,parameter::FBa=106,FF=107,FSa=108,FE=109,FI=110,FO=111
integer,parameter::SB=112,SF=113,SS=114,SE=115,SI=116,SO=117
integer,parameter::EB=118,EF=119,ES=120,EE=121,EI=122,EO=123
integer,parameter::IB=124,IFa=125,IS=126,IE=127,II=128,IO=129
integer,parameter::OB=130,OFa=131,OS=132,OE=133,OI=134,OO=135
end module const

module apert
use const
real(8) ax(imax,jmax),ay(imax,jmax),fb(imax,jmax),fs(imax,jmax)
real(8) axn(imax,jmax),ayn(imax,jmax),fsn(imax,jmax)
end module apert

module times
integer t_cycle,t_cycle2
real(8) T,DelT,DelTMx,TFin,TStart
real(8) T2,DelT2,DelTMx2,TFin2,TStart2
end  module times

module gridar
use const
real(8) x(imax),xi(imax),DelX(imax),y(jmax),yj(jmax),DelY(jmax)
real(8) rdx(imax),rdy(jmax),cx,cy   !,dx,dy
integer IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
end module gridar

module phys
use const
real(8) u(0:imax,jmax),un(0:imax,jmax)
real(8) v(imax,0:jmax),vn(imax,0:jmax)
real(8) p(imax,jmax),pn(imax,jmax),vmax,pmin,pmax
end module phys

module dxrl
use const
real(8) dxri(imax),dxli(imax),dytj(jmax),dybj(jmax)
end module dxrl

module bndcon
real(8) SgnB,SgnL,SgnR,SgnT,UIn,VIn,FreqIn
real(8) PoutN,PoutE,POutS,POutW
integer WB,WL,WR,WT
end module bndcon

module tank
real(8) Xmax,Xmin,Ymax,Ymin
end module tank

module orga2
integer IMilu
integer Iter,Itmax,ItSum,StrtP
real(8) epsi,OmStrt,Alpha,NrmRhs
end module orga2

module condif
use const
integer divuv
real(8) ru(imax,jmax),rv(imax,jmax),divmax
end module condif

module adjust
real(8) CFLMin,CFLMax
integer NDTAd,MxDTAd,SlwFlw
real(8) CFLMin2,CFLMax2
integer NDTAd2,MxDTAd2,SlwFlw2
end module adjust

module fluid
real(8) rho,mu,nu,sigma
end module fluid

module prints
REAL(8) TPrt, PrtDT, TAdPrt, TCnf(2001)
INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce
end module prints

module rek
REAL(8) xpos,ypos
!REAL(8) posx,posy,xpos,ypos
end module rek

module coef2
use const
REAL(8) CWI(IMAX),CEI(IMAX),CNJ(JMAX),CSJ(JMAX)
end module coef2

module orga
use const
INTEGER NF(IMAX,JMAX),NFN(IMAX,JMAX)
end module orga

module rel
integer relatief,relatief2
end module rel

module coefp
use const
REAL(8) CN(IMAX,JMAX), CS(IMAX,JMAX), CE(IMAX,JMAX)
REAL(8)     CW(IMAX,JMAX), CC(IMAX,JMAX), DIV(IMAX,JMAX)
end module coefp

module workar
use const
REAL PSI(IMAX,JMAX)
end module workar

module label
use const
real(8)::ulabel(imax,jmax),ulabfs(imax,jmax),vlabel(imax,jmax),vlabfs(imax,jmax)
real(8)::plabel(imax,jmax),plabfs(imax,jmax),plabfsn(imax,jmax)
end module label







