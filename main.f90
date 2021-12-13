program main
!use phys
!use gridar
!!!!!!!!!!!!!!!!!!!!!!!!!!
use times
use orga2
use rel
use prints
character(200) filename1,filename2,filename3
character(50) szNumber
integer IOpst,i
REAL HalfDT,HalfDT2
INTEGER FlgAdj,FlgAdj2
IOpst=1
!filename1=Trim("../vof2d16_ac3/data/")//Trim("grid")//Trim(".vtk")
call setpar
call setfld

call diverr

TStart = T
TPrt = T
TAdPrt = T

TStart2 = T2
TPrt2 = T2
TAdPrt2 = T2

CALL PRT(1)
FlgAdj=0
t_Cycle=0

FlgAdj2=0
t_Cycle2=0
100 CONTINUE 
!внешние итерации по времени
T = T + DelT
t_Cycle = t_Cycle + 1 
HalfDT = 0.5*DelT
10 CONTINUE 
!внутренние итерации по псевдовремени 
T2 = T2 + DelT2
t_Cycle2 = t_Cycle2 + 1 

HalfDT2 = 0.5*DelT2
 call init





 call bcbnd

 
 IF (FlgAdj2.NE.0) relatief2 = 1.0  
 call tilde

 call solvep

 call bcbnd

 CALL CFLCHK2(FlgAdj2) 


 IF (FlgAdj2.EQ.0) THEN  
   call diverr
   !IF (T2+HalfDT2 .GT. TPrt) CALL PRT(2)
   !if(mod(t_cycle2,100)==0) then

    !Write(szNumber,'(i6.6)') t_cycle2
    !filename1=Trim("../vof2d16_ac4/vtk/")//Trim("uvpf")//Trim(szNumber)//Trim(".vtk")
    !filename3=Trim("../vof2d16_ac4/data/")//Trim("uvpf")//Trim(szNumber)//Trim(".dat")
	
	!open(100,file=filename3)
	!call output(100)
	!close(100)
    !call vtk(IOpst,filename1,filename2)
	
  !end if
 end if 
 IF (t_cycle2.LT.4) relatief2 = 1.0
         

IF (T2+HalfDT2 .LT. TFin2 .AND.&
            relatief2.GE. -1e-6*delT2) GOTO 10         

IF (FlgAdj.NE.0) relatief = 1.0 
CALL CFLCHK(FlgAdj)
IF (FlgAdj.EQ.0) THEN  
   call diverr
   IF (T+HalfDT .GT. TPrt) CALL PRT(2)
   if(mod(t_cycle,100)==0) then

    Write(szNumber,'(i6.6)') t_cycle
    filename1=Trim("../vof2d16_ac5/vtk/")//Trim("uvpf")//Trim(szNumber)//Trim(".vtk")
    !filename3=Trim("../vof2d16_ac5/data/")//Trim("uvpf")//Trim(szNumber)//Trim(".dat")
	
	!open(100,file=filename3)
	!call output(100)
	!close(100)
    call vtk(IOpst,filename1,filename2)
	
  end if
 end if 
IF (t_cycle.LT.4) relatief = 1.0 
IF (T+HalfDT .LT. TFin .AND.&
            relatief.GE. -1e-6*delT) GOTO 100 
end



