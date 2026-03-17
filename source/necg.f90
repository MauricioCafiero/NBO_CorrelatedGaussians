!PROJECT TITLE: N-electron correlated gaussians
!
!     AUTHOR: Mauricio Cafiero
!     PROJECT SUMMARY: N-electron correlated gaussians, born-oppenheimer, using shifted centers.
PROGRAM  necg
!SHARED MODULE ***********************************************************************************
USE MATRICES_TO_BE_SHARED



!Declarations ************************************************************************************
IMPLICIT NONE 

include 'mpif.h'


DOUBLE PRECISION:: DP,dp2,dp3,dp4,d1,d2,d3,d4,d5,d6,d7,d8,d9,store,eps,ETA,STEPmx,ACCRCY,XTOL	
double precision, allocatable,dimension(:)::hiv	,w,hiv2,gr,v1,v2,work
INTEGER::I,J,JJ	,l,k,lw,ii,MAXIT,MSGLVL,MAXFUN,kk,stf,lwork,hostnm,stfe	
double precision, allocatable,dimension(:,:)::ea,eb,iden		

CHARACTER (LEN=10) :: TT			
CHARACTER (LEN=5)  :: TZ
CHARACTER (LEN=8)  :: TD
CHARACTER (LEN=80) :: CHARSTRING,recode	

external outprod

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!***************************************************************************************************
!Prepare output
  
!!$flag=hostnm(hostname)
!!$write(*,*)'hello from ',hostname

if(rank==0)then

   OPEN (7, FILE = 'necg.out', ACCESS = 'SEQUENTIAL')	

   CALL BIOUT('N-particle Non-Born-Oppenheimer Correlated Gaussians')
   CALL BIOUT('===================================='//'===============')
   CALL BIOUT('ADAMOWICZ RESEARCH                 '  &
  			//'[UofA chemistry]')
   CALL BIOUT('Programming: Mauricio Cafiero ') 
   CALL BIOUT('===================================='  &
  	               //'===============')


  													
   CALL BIOUT(':PROGRAM BEGINS:')				
   CALL Get_Current_Time(TD,TT,TZ)
   CHARSTRING='Date: '//TD//' | Time: '//TT//' / '//TZ//' gmt'
   CALL BIOUT(CHARSTRING)
end if

call read_input
 				
if(rank==0)then
DO I=6,7
  	WRITE(I,*)'M, NN, NE,NST,  =',M,ne, nst
  	WRITE(I,*)'Electric Field  =',efs
  	WRITE(I,*)'OUTPUT FILENAME     =','necg.out'
  
END DO
end if


lwork=3*m-1

Allocate(hiv(siv+m),hiv2(siv+m),ea(m,m),eb(m,m),work(lwork),  &
         v1(svh+skp),v2(svh+skp),gr(siv+m),iden(m,m))


!**************


!!$call matel(1,1,2,d1,d2,d3,d4)
!!$write(*,*)'hello'
!!$d5=d4
!!$eps=0.00001
!!$c(1)=c(1)+eps
!!$iv(91)=iv(91)+eps
!!$call matel(1,1,2,d1,d2,d3,d4)
!!$write(*,*)'hello'
!!$call gradel(1,1,2,v1,v2)
!!$d6=(d4-d5)/eps
!!$write(*,*)'gradient:',d6
!!$stop

!******************************


h=zero
s=zero
kinetic=zero
potential=zero
ki=zero
po=zero
hh=zero
ss=zero

if(rank==0)then 
d5=MPI_WTIME() 
end if       

do i=rank+1,m,nproc
   do j=i,m
      do stf=1, nst
      call matel(stf,j,i,dp,dp2,dp3,dp4)
      kinetic(j,i)=kinetic(j,i)+scvv(stf)*dp2
      potential(j,i)=potential(j,i)+scvv(stf)*(dp3)
      h(j,i)= h(j,i)+scvv(stf)*dp2+scvv(stf)*(dp3)-scvv(stf)*dp4 
      s(j,i)=s(j,i)+scvv(stf)*dp
   end do
end do
end do

call MPI_REDUCE(kinetic,ki,m*m, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(potential,po,m*m,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(h,hh,m*m, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(s,ss,m*m, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

if(rank==0)then
d6=MPI_WTIME()
d7=d6-d5
write(*,*)'time for matrix calc. is',d7,'seconds'
call flush(6)
end if
!!$goto 97

if(rank==0)then
s=ss
h=hh
kinetic=ki
potential=po
end if
!!$!stop

if(rank==0)then
ea=h
eb=s
!****
!!$iden=zero
!!$do i=1,m
!!$   iden(i,i)=one
!!$end do
!!$allocate(w(m))
!!$call DSYGV(1,'V','L',m,eb,m,iden,m,w,work,lwork,jj)
!!$write(*,*)'error from eigenvalue:',jj
!!$write(*,*)'Eigenvalue=',w
!!$deallocate(w,work)
!!$do i=1,m
!!$   write(*,*)eb(i,1)
!!$end do
!!$end if
!!$stop
!!$if(rank==0)then
!****
do j=1,m
   do i=1,j
      h(i,j)=h(j,i)
      s(i,j)=s(j,i)
      kinetic(i,j)=kinetic(j,i)
      potential(i,j)=potential(j,i)
   end do
end do

!!$!call tab(s,m,m,m,m)
!stop
!goto 16
allocate(w(m))
call DSYGV(1,'V','L',m,ea,m,eb,m,w,work,lwork,jj)
write(*,*)'error from eigenvalue:',jj
write(*,*)'Eigenvalue=',w(1)
do i=1,m
   c(i)=ea(i,1)
end do
deallocate(w,work)

16 j=1
do i=1,siv
   hiv(i)=iv(i)
   j=j+1
end do
do i=1,m
   hiv(j)=c(i)
   j=j+1
end do

DO I=6,7
  	WRITE(I,*)'********** The Input Vectors ****************************************'
  	WRITE(I,*)''
END DO
!write(*,*)iv
!write(*,*)''
write(*,*)c
DO I=6,7
  	WRITE(I,*)'********** End Input Vectors ****************************************'
  	WRITE(I,*)''
END DO

d1=zero
d2=zero
d4=zero
d5=zero
do i=1,m
   do j=1,m
      d1=d1+c(i)*h(i,j)*c(j)
      d2=d2+(c(i)*s(i,j)*c(j))
      d4=d4+c(i)*kinetic(i,j)*c(j)
      d5=d5+c(i)*potential(i,j)*c(j)
   end do
end do

d3=d1/d2
d4=d4/d2
d5=d5/d2
d6=-two*d4/d5

write(*,*)'*********************************************************************'
write(*,*)''
write(*,*)'Rayleigh Quotient: ', d3
write(*,*)''
write(*,*)'Average Kinetic  : ', d4
write(*,*)''
write(*,*)'Average Potential: ', d5
write(*,*)''
write(*,*)'Virial Coef.     : ', d6
write(*,*)''
write(*,*)'*********************************************************************'

end if
if(rank==0)then
d5=MPI_WTIME()
end if
call outprod(siv+m,hiv,d1,gr)
if(rank==0)then
d6=MPI_WTIME()
d7=d6-d5
write(*,*)'time for energy/gradient calc. is',d7,'seconds'
write(*,*)'energy is ',d1,'hartree'
call flush(6)
end if
call MPI_BCAST(hiv,siv+m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d3,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!****************  numerical optimization ***********
!l=siv+m
!lw=l*(l+10)
!lw=lw+10
!write(*,*)'lw:',lw
!allocate(w(lw))
!write(*,*)'hello'
!call uncmnd(l,hiv,outprod,hiv2,d3,jj,w,lw)

!j=1
!do i=1,siv
!   iv(i)=hiv2(i)
!   j=j+1
!end do
!do i=1,m
!   c(i)=hiv2(j)
!   j=j+1
!end do
!goto 23
!****************************************************

!***************** analytical optimization **********

!Set up parameters for the optimization program
      MAXIT = 20
      MSGLVL = 1
      MAXFUN = 150*(siv+m)
      ETA = 0.21D0
      STEPMX = 0.7D1
      ACCRCY = 2.220446649250313D-16
      XTOL = DSQRT(ACCRCY)
      l=(siv+m)
      lw=14*l
      allocate(w(lw))

CALL LMQN(jj,l,hiv,d3,Gr,W,LW,outprod,MSGLVL,MAXIT,MAXFUN,ETA,STEPMX,ACCRCY,XTOL,rank)
  

if(rank==0)then

wRITE(*,*)"IERROR",jj
j=1
do i=1,siv
   iv(i)=hiv(i)
   j=j+1
end do
do i=1,m
   c(i)=hiv(j)
   j=j+1
end do
!end if
!***********************************************************************

23 DO I=6,7
  	WRITE(I,*)'********** The Output Vectors ***************************************'
  	WRITE(I,*)''
END DO
!write(*,*)iv
!write(*,*)''
write(*,*)c
write(*,*)''
DO I=6,7
  	WRITE(I,*)'********** End Output Vectors ***************************************'
  	WRITE(I,*)''
END DO

!DO I=6,7
!  	WRITE(I,*)'********** The Final Gradient Vector ********************************'
!  	WRITE(I,*)''
!END DO
!write(*,*)gr
!write(*,*)''
!DO I=6,7
! 	WRITE(I,*)'*********************************************************************'
!  	WRITE(I,*)''
!END DO

!call tab (H,m,m,m,m)
h=zero
s=zero
kinetic=zero
potential=zero
do stf=1,nst
!stf=1
do i=1,m
   do j=i,m
      call matel(stf,j,i,dp,dp2,dp3,dp4)
      kinetic(j,i)=kinetic(j,i)+scvv(stf)*dp2
      potential(j,i)=potential(j,i)+scvv(stf)*(dp3)
      h(j,i)= h(j,i)+scvv(stf)*dp2+scvv(stf)*(dp3)-scvv(stf)*dp4
      s(j,i)=s(j,i)+scvv(stf)*dp
   end do
end do
end do


!!$!stop
do j=1,m
   do i=1,j
      h(i,j)=h(j,i)
      s(i,j)=s(j,i)
      kinetic(i,j)=kinetic(j,i)
      potential(i,j)=potential(j,i)
   end do
end do

d1=zero
d2=zero
d4=zero
d5=zero
do i=1,m
   do j=1,m
      d1=d1+c(i)*h(i,j)*c(j)
      d2=d2+(c(i)*s(i,j)*c(j))
      d4=d4+c(i)*kinetic(i,j)*c(j)
      d5=d5+c(i)*potential(i,j)*c(j)
   end do
end do

d3=d1/d2
d4=d4/d2
d5=d5/d2
d6=-two*d4/d5
!!$
!!$
write(*,*)'********************************'
write(*,*)''
write(*,*)'Rayleigh Quotient: ', d3
write(*,*)''
write(*,*)'Average Kinetic  : ', d4
write(*,*)''
write(*,*)'Average Potential: ', d5
write(*,*)''
write(*,*)'Virial Coef.     : ', d6
write(*,*)''
write(*,*)'********************************'

d9=zero
d7=zero
d8=zero
kk=1
do ii=1,ne
do jj=1,3


h=zero
s=zero

do stfe=1,nst
do stf=1,nst
!stf=1
do i=1,m
   do j=i,m

      call rrms(stfe,stf,j,i,kk,dp,dp2)
      h(j,i)= h(j,i)+scvv(stf)*scvv(stfe)*dp2
      s(j,i)=s(j,i)+scvv(stf)*scvv(stfe)*dp
   end do
end do
end do
end do

do j=1,m
   do i=1,j
      h(i,j)=h(j,i)
      s(i,j)=s(j,i)
   end do
end do
!!$call biout('Hamiltonian Matrix')
!!$call tab(h,m,m,m,m)
!!$call biout('Overlap Matrix')
!!$call tab(s,m,m,m,m)
!!$stop


d1=zero
d2=zero
do i=1,m
   do j=1,m
      d1=d1+c(i)*h(i,j)*c(j)
      d2=d2+(c(i)*s(i,j)*c(j))
   end do
end do
!!$write(*,*)d1,d2
d3=d1/d2

if(jj==1)then
   d7=d7+ncv(ii+1)*d3
else if(jj==2)then
   d8=d8+ncv(ii+1)*d3
else if(jj==3)then
   d9=d9+ncv(ii+1)*d3
end if


kk=kk+1
end do
end do

d4=d7*2.541765d+000
d5=d8*2.541765d+000
d6=d9*2.541765d+000

write(*,*)''
write(*,*)'Dipole x-component=', d7,'(',d4,')'
write(*,*)''

write(*,*)''
write(*,*)'Dipole y-component=', d8,'(',d5,')'
write(*,*)''

write(*,*)''
write(*,*)'Dipole z-component=', d9,'(',d6,')'
write(*,*)''

write(*,*)'*********************************************************************'

  CALL BIOUT('PROGRAM TERMINATES')
  CALL Get_Current_Time(TD,TT,TZ)
  CHARSTRING='Date: '//TD//' | Time: '//TT//' / '//TZ//' gmt'
  CALL BIOUT(CHARSTRING)

end if

call MPI_FINALIZE(ierr)

end program necg
 









