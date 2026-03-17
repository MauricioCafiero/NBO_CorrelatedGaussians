!SUBROUTINE outprod(n,x,f)
subroutine outprod(n,x,f,g)
use matrices_to_be_shared
implicit none
include 'mpif.h'

integer::n
double precision::x(n),g(n),f
integer::i,j,stf      
double precision::dp,dp2,dp3,dp4,d1,d2,d3      

!write(*,*)'hello'

if(rank==0)then
j=1
do i=1,siv
   iv(i)=x(i)
   j=j+1
end do
do i=1,m
   c(i)=x(j)
   j=j+1
end do
end if

call MPI_BCAST(iv,siv,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(c,m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

h=zero
s=zero
hh=zero
ss=zero

do i=rank+1,m,nproc
   do j=1,i
      do stf=1,nst	
      call matel(stf,i,j,dp,dp2,dp3,dp4)
      h(i,j)= h(i,j)+scvv(stf)*dp2+scvv(stf)*(dp3)-scvv(stf)*dp4
      s(i,j)=s(i,j)+scvv(stf)*dp
   end do
end do
end do

call MPI_REDUCE(h,hh,m*m, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(s,ss,m*m, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
if(rank==0)then
s=ss
h=hh

do j=1,m
   do i=1,j
      h(i,j)=h(j,i)
      s(i,j)=s(j,i)
   end do
end do


d1=zero
d2=zero
do i=1,m
   do j=1,m
      d1=d1+c(i)*h(i,j)*c(j)
      d2=d2+(c(i)*s(i,j)*c(j))
   end do
end do

d3=d1/d2


f=d3
end if
call MPI_BCAST(f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!write(*,*)f,d2, rank
!stop

d2=one/d2
!goto 25
call make_grad(f,d2)

if(rank==0)then
do i=1,(siv+m)
   g(i)=grad(i)
end do

!stop
!*************** inserting a shift and scale factor
!f=(f*1.0d+002)+1.15d+002
!g=g*1.0d+002
!**************************************************
end if
call MPI_BCAST(f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(g,siv+m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!$if(rank==0)then
!!$write(*,*)g
!!$end if
!!$stop
25 END SUBROUTINE outprod
