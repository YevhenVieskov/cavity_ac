subroutine vtk(IOpst,filename1,filename2)
use gridar
use phys
use apert
use const
implicit none
integer IOpst,ios
character(200),intent(in)::filename1,filename2
real(8) utemp(imax,jmax),vtemp(imax,jmax)
integer i,j
!header vtk file
open(IOpst,file=filename1,status='new',iostat=ios)
write(IOpst,'(A)')   '# vtk DataFile Version 2.0'
write(IOpst,'(A,A,A)')   'vof result output'
write(IOpst,'(A)')      'ASCII'

!grid data

write(IOpst,'(A)')      'DATASET RECTILINEAR_GRID' 
write(IOpst,'(A,1X,I5,1X,I5,1X,I5)')      'DIMENSIONS',IM1Us,JM1Us,1
write(IOpst,'(A,1X,I5,1X,A)')      'X_COORDINATES',IM1Us,'double'  
write(IOpst,*)   X(1:IM1Us)  
write(IOpst,'(A,1X,I5,1X,A)')      'Y_COORDINATES',JM1Us,'double' 
write(IOpst,*)   Y(1:JM1Us)  
write(IOpst,'(A,1X,I5,1X,A)')      'Z_COORDINATES',1,'double'
write(IOpst,*)   0.0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
write(IOpst,'(A,1X,I5)')   'POINT_DATA', (IM1Us)*(JM1Us)
write(IOpst,'(A,1X,A,1X,A)') 'VECTORS', 'velocity_point', 'double'
utemp=0d0
vtemp=0d0

do i=1,IM1Us
  do j=1,JM1Us !!!
     utemp(i,j)=0.5d0*(u(i,j+1)+u(i,j))
  end do
end do

do i=1,IM1Us !!!
  do j=1,JM1Us
     vtemp(i,j)=0.5d0*(v(i+1,j)+v(i,j))
  end do
end do

do j=1,JM1Us
  do i=1,IM1Us 
    write(IOpst,'(f8.3,1X,f8.3,1X,f8.3)') utemp(i,j),vtemp(i,j),0.0  
  end do
end do

write(IOpst,'(A,1X,I5)')   'CELL_DATA', (IM2Us)*(JM2Us)

write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'pressure', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=2,JM1Us
  do i=2,IM1Us 
    write(IOpst,'(f20.3)') p(i,j)
  end do
end do

write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'vof_fraction', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=2,JM1Us
  do i=2,IM1Us  
    write(IOpst,'(f8.3)') fs(i,j)
  end do
end do  

end subroutine vtk



