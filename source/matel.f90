SUBROUTINE matel(stf,ii,jj,OV,tk,ert,zht)
!**********************************
!4/28/00 -- ov, tk, na, er are great; 
!	    

USE matrices_to_be_shared
implicit none		

double precision::zht
double precision::ov,tk,er,na,i1,i2,i3,j1,j2,j3,detakl,detak,detal,dett,f0,df0,ert,nat
double precision, dimension(ne,ne)  :: ak,al,lk,lkt,ell,ellt,alt		  
double precision, dimension(ne,ne)  ::akl,akli,jij,ism1,jije,aki,ism2,ism3,jii
double precision, dimension(skp,2)  ::si
double precision, dimension(3,ne)   ::dot1,dot1e
double precision, dimension(ne,3)   ::dot2,dot2e
double precision, dimension(skp,skp)::akb,alb,aklib,im1,im2,im3,iden,jijb,aklp,aklpt,jiib,massb
double precision, dimension(skp)    ::iv1,iv2,ee,iv3,st,st1,st2,en,nv
double precision, dimension(svh)    ::vech1
integer::zq,I,J,k,ii,jj,kk,ll,pp,stf
!******** make ak, al, sk, sl **************************

lk=zero
lkt=zero
ell=zero
ellt=zero
kk=(ii-1)*(svh+skp)+1
ll=(jj-1)*(svh+skp)+1
do i=1,ne
   do j=i,ne
      lk(j,i)=iv(kk)
      ell(j,i)=iv(ll)
      kk=kk+1
      ll=ll+1
   end do
end do

do j=1,ne 
     do i=1,j
      lkt(i,j)=lk(j,i)
      ellt(i,j)=ell(j,i)
   end do
end do

ak=matmul(lk,lkt)
al=matmul(ell,ellt)
!call tab(lk,ne,ne,ne,ne)
!stop
!******* perform symmetry operation on the ket
alt=zero
do i=1,ne
   do j=1,ne
      do k=1,ne
         alt(i,j)=alt(i,j)+al(i,k)*sml(k,j,stf)
      end do
   end do
end do

al=zero
do i=1,ne
   do j=1,ne
      do k=1,ne
         al(i,j)=al(i,j)+sml(k,i,stf)*alt(k,j)
      end do
   end do
end do
!call tab(al,ne,ne,ne,ne)
!stop
!********************************************
!!$call tab(ak,ne,ne,ne,ne)
!call tab(al,ne,ne,ne,ne)

kk=(ii-1)*(svh+skp)+svh+1
ll=(jj-1)*(svh+skp)+svh+1
do i=1,skp
   si(i,1)=iv(kk)
   si(i,2)=iv(ll)
   kk=kk+1
   ll=ll+1
end do

!write(*,*)si(1,2), si(2,2),si(3,2), si(4,2),si(5,2), si(6,2),si(7,2), si(8,2),si(9,2)
!******* perform symmetry operation on the ket shift vector
iv1=zero
do i=1,skp
   do j=1,skp
      iv1(i)=iv1(i)+sms(i,j,stf)*si(j,2)
   end do
end do
do i=1,skp
   si(i,2)=iv1(i)
end do
!***********************************************************
!write(*,*)si(1,2), si(2,2),si(3,2), si(4,2),si(5,2), si(6,2),si(7,2), si(8,2),si(9,2)
!stop
!*** st's defined here 
do i=1,skp
   st1(i)=si(i,1)
   st2(i)=si(i,2)
end do



akb=zero
alb=zero
aklib=zero
massb=zero

!write(*,*)st1
!write(*,*)st2
!!$call tab(ak,2,2,2,2)
!!$call tab(al,2,2,2,2)

akl=ak+al
akli=akl
aki=ak
call DPOTRF('l',ne,akli,ne,jjj)
call DPOTRF('l',ne,aki,ne,jjj)
dett=one
do i=1,ne
	dett=dett*akli(i,i)
end do
dett=dett**two
call DPOTRI('l',ne,akli,ne,jjj)
call DPOTRI('l',ne,aki,ne,jjj)

do i=1,ne
	do j=(i+1),ne
	aki(i,j)=aki(j,i)
	akli(i,j)=akli(j,i)
	end do
end do

!****************************************
!do the kronecker product

do i=1,ne
   kk=(i-1)
   do j=1,ne
      zq=(j-1)
      do k=1,3
      akb(3*kk+k,3*zq+k)=ak(i,j)
      massb(3*kk+k,3*zq+k)=mass(i,j)
      alb(3*kk+k,3*zq+k)=al(i,j)
      aklib(3*kk+k,3*zq+k)=akli(i,j)
      end do
   end do
end do
!call tab(ak,3,3,3,3)
!call tab(akb,9,9,9,9)
!stop

!******************************************

detal=one   
detak=one
do i=1,ne
	detal=detal*ell(i,i)
	detak=detak*lk(i,i)
end do
detal=detal**two
detak=detak**two
   
detal=dsqrt(detal)
detak=dsqrt(detak)



!****************************************************
! make it agree with L&K



!****************************************************
!   Do the integrals
!
!overlap integral 
   
iv1=zero
i1=zero
do i=1,skp
   do j=1,skp
      iv1(i)=iv1(i)+akb(i,j)*st1(j)
   end do
end do
do i=1,skp
   i1=i1+st1(i)*iv1(i)
end do

iv2=zero
i2=zero
do i=1,skp
   do j=1,skp
      iv2(i)=iv2(i)+alb(i,j)*st2(j)
   end do
end do
do i=1,skp
   i2=i2+st2(i)*iv2(i)
end do

!******* ee is defined for here on out
ee=iv1+iv2

iv3=zero
i3=zero
do i=1,skp
   do j=1,skp
      iv3(i)=iv3(i)+aklib(i,j)*ee(j)
   end do
end do
do i=1,skp
   i3=i3+ee(i)*iv3(i)
end do

!write(*,*)i3
st=iv3

ov=exp(-(i1+i2)+i3)*dsqrt(eight**ne)*dsqrt((detak*detal/dett)**3)
!!$write(*,*)i1,i2,i3
!!$ov=exp(-(i1+i2)+i3)

!********************************************************
!Kinetic energy 

iv1=zero
iv2=zero
i1=zero
i2=zero
i3=zero


im1=matmul(massb,alb)
im2=matmul(akb,im1)

do i=1,skp
   do j=1,skp
      i3=i3+(st(i)-st1(i))*im2(i,j)*(st(j)-st2(j))
   end do
end do


!!$do i=1,skp
!!$   do j=1,skp
!!$      iv1(i)=iv1(i)+alb(i,j)*(st(j)-st2(j))
!!$   end do
!!$end do
!!$
!!$
!!$do i=1,skp
!!$   do j=1,skp
!!$      iv2(i)=iv2(i)+massb(i,j)*iv1(j)
!!$   end do
!!$end do
!!$
!!$do i=1,skp
!!$   do j=1,skp
!!$      iv1(i)=iv1(i)+akb(i,j)*iv2(j)
!!$   end do
!!$end do
!!$
!!$do i=1,skp
!!$   i3=i3+(st(i)-st1(i))*iv1(i)
!!$end do

im1=matmul(aklib,alb)
im2=matmul(akb,im1)
im1=matmul(massb,im2)
im2=im1

do i=1,skp
   i2=i2+im2(i,i)
end do


tk=two*ov*(two*i3+i2)


!********************************************
! electron replsion -- associated gradients

er=zero
ert=zero

do pp=0,ne
   do ll=pp,ne
      if (pp==ll)then
         er=er+zero
      else   
!******* fill in Jij ***********
jij=zero
if(pp==0)then
   jij(ll,ll)=one
else
jij(pp,pp)=one
jij(ll,ll)=one
jij(pp,ll)=-one
jij(ll,pp)=-one
jije=jij
end if
jijb=zero

!*** form product akli*jij*akli and store in jij ********
ism1=matmul(jij,akli)
jij=matmul(akli,ism1)

!****************************************
!do the kronecker product

do i=1,ne
   kk=(i-1)
   do j=1,ne
      zq=(j-1)
      do k=1,3
      jijb(3*kk+k,3*zq+k)=jij(i,j)
      end do
   end do
end do


!*************************************
!electron repulsion integral proper

iv3=zero
i3=zero
do i=1,skp
   do j=1,skp
      iv3(i)=iv3(i)+jijb(i,j)*ee(j)
   end do
end do
do i=1,skp
   i3=i3+ee(i)*iv3(i)
end do

i2=zero
do i=1,ne
   i2=i2+ism1(i,i)
end do

i1=i3/i2

er=ncv(pp+1)*ncv(ll+1)*two*ov*f0(i1)/dsqrt(i2*pi)
ert=ert+er
!write(*,*)'er,ert',er,ert


end if
end do
end do
!stop
zht=zero
do i=1,ne
	zht=zht+ncv(i+1)*st(3*i)
end do
zht=ov*zht*efs

END SUBROUTINE matel
