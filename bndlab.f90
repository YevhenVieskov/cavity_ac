subroutine bndlab
!define edge apertures, cell label, velocity label
use gridar
use phys
use label
use orga
implicit none
!local variables
integer i,j
do i=1,IMaxUs
  do j=1,JMaxUs
    !geometry label
    if(nf(i,j)==B.and.nf(i+1,j)==B) ulabel=BB
	if(nf(i,j)==B.and.nf(i,j+1)==B) vlabel=BB

    if(nf(i,j)==B.and.nf(i+1,j)==F) ulabel=BF
	if(nf(i,j)==B.and.nf(i,j+1)==F) vlabel=BF

    if(nf(i,j)==B.and.nf(i+1,j)==S) ulabel=BS
	if(nf(i,j)==B.and.nf(i,j+1)==S) vlabel=BS

	if(nf(i,j)==B.and.nf(i+1,j)==E) ulabel=BE
	if(nf(i,j)==B.and.nf(i,j+1)==E) vlabel=BE

    if(nf(i,j)==F.and.nf(i+1,j)==B) ulabel=FBa
	if(nf(i,j)==F.and.nf(i,j+1)==B) vlabel=FBa

    if(nf(i,j)==S.and.nf(i+1,j)==B) ulabel=SB
	if(nf(i,j)==S.and.nf(i,j+1)==B) vlabel=SB

	if(nf(i,j)==E.and.nf(i+1,j)==B) ulabel=EB
	if(nf(i,j)==E.and.nf(i,j+1)==B) vlabel=EB
    !surface label
    if(nf(i,j)==F.and.nf(i+1,j)==S) ulabfs=FSa
	if(nf(i,j)==F.and.nf(i,j+1)==S) vlabfs=FSa

	if(nf(i,j)==I.and.nf(i+1,j)==S) ulabfs=IS
	if(nf(i,j)==I.and.nf(i,j+1)==S) vlabfs=IS

	if(nf(i,j)==OU.and.nf(i+1,j)==S) ulabfs=OS
	if(nf(i,j)==OU.and.nf(i,j+1)==S) vlabfs=OS

	if(nf(i,j)==S.and.nf(i+1,j)==F) ulabfs=SF
	if(nf(i,j)==S.and.nf(i,j+1)==F) vlabfs=SF

	if(nf(i,j)==S.and.nf(i+1,j)==I) ulabfs=SI
	if(nf(i,j)==S.and.nf(i,j+1)==I) vlabfs=SI

	if(nf(i,j)==S.and.nf(i+1,j)==OU) ulabfs=SO
	if(nf(i,j)==S.and.nf(i,j+1)==OU) vlabfs=SO

    if(nf(i,j)==F.and.nf(i+1,j)==F) ulabfs=FF
	if(nf(i,j)==F.and.nf(i,j+1)==F) vlabfs=FF

	if(nf(i,j)==I.and.nf(i+1,j)==I) ulabfs=II
	if(nf(i,j)==I.and.nf(i,j+1)==I) vlabfs=II

	if(nf(i,j)==OU.and.nf(i+1,j)==OU) ulabfs=OO
	if(nf(i,j)==OU.and.nf(i,j+1)==OU) vlabfs=OO
  end do
end do
end subroutine bndlab