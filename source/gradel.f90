SUBROUTINE gradel(stf,ii,jj,tgvl,tgvls)
!***********************************
!4/28/00 -- ov, tk, na, er are great; 
!	    5/6/00 -- error is in electron repulsion;


USE matrices_to_be_shared
implicit none		

double precision::zht
double precision::ov,tk,er,na,i1,i2,i3,j1,j2,j3,detakl,detak,detal,dett, f0,df0,ert,nat
double precision, dimension(ne,ne)  :: ak,al,lk,lkt,ell,ellt,aklps,aklpst,sds,sdsk,slds,sldsk,skdsk		  
double precision, dimension(ne,ne)  ::akl,akli,jij,ism1,jije,aki,ism2,ism3,jii,kcm,alt
double precision, dimension(skp,2)  ::si
double precision, dimension(3,ne)   ::dot1,dot1e
double precision, dimension(ne,3)   ::dot2,dot2e
double precision, dimension(skp,skp)::akb,alb,aklib,im1,im2,im3,iden,jijb,aklp,aklpt,jiib,massb
double precision, dimension(skp)    ::iv1,iv2,ee,iv3,st,st1,st2,ogs,kgs,ergs,en,nv,nags
double precision, dimension(skp)    ::ergst,nagst,igvs,zgs,voo
double precision,dimension(svh)::vech1,ogl,kgl,ergl,nagl,naglt,erglt,igvl,zgl
double precision,dimension(svh)::vech2
double precision, dimension(svh+skp)    ::tgvl,tgvls
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
!********************************************


!!$call tab(ak,ne,ne,ne,ne)
!!$call tab(al,ne,ne,ne,ne)

kk=(ii-1)*(svh+skp)+svh+1
ll=(jj-1)*(svh+skp)+svh+1
do i=1,skp
   si(i,1)=iv(kk)
   si(i,2)=iv(ll)
   kk=kk+1
   ll=ll+1
end do

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

!*** st's defined here 
do i=1,skp
   st1(i)=si(i,1)
   st2(i)=si(i,2)
end do


akb=zero
alb=zero
aklib=zero
massb=zero

!!$write(*,*)st1
!!$write(*,*)st2
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

!***********************************************************************
!derivative of ovrlap with respect to sk --good

ogs=zero
kgs=zero

do i=1,skp
   do j=1,skp
      ogs(i)=ogs(i)+akb(i,j)*(st(j)-st1(j))
   end do
end do
ogs=ogs*two*ov

!*************************************
!derivative of Skl wrt vech(Lk') -- numerical testing good!!!

ism1=matmul(aki,lk)
call vech(ism1,ne,vech1)
ogl=three*vech1/two

ism1=matmul(akli,lk)
call vech(ism1,ne,vech1)
ogl=ogl-three*vech1


!**** redefine st's as just the shifts
do i=1,skp
   st1(i)=si(i,1)
   st2(i)=si(i,2)
end do
k=1
do i=1,ne
   do j=1,3
      dot1(j,i)=st1(k)
      k=k+1
   end do
end do

dot2=transpose(dot1)
ism1=matmul(dot2,dot1)
ism2=matmul(ism1,lk)
call vech(ism2,ne,vech1)
ogl=ogl-two*vech1

k=1
do i=1,ne
   do j=1,3
      dot1(j,i)=ee(k)
      k=k+1
   end do
end do

dot2=transpose(dot1)
ism1=matmul(dot2,dot1)
ism2=matmul(ism1,akli)
ism1=matmul(akli,ism2)
ism2=matmul(ism1,lk)
call vech(ism2,ne,vech1)
ogl=ogl-two*vech1

k=1
do i=1,ne
do j=1,3
dot1(j,i)=st1(k)
k=k+1
end do
end do
ism1=matmul(dot2,dot1)
ism2=matmul(akli,ism1)
ism3=matmul(ism2,lk)
call vech(ism3,ne,vech1)
ogl=ogl+two*vech1


ism3=transpose(ism2)
ism1=matmul(ism3,lk)
call vech(ism1,ne,vech1)
ogl=ogl+two*vech1

ogl=ov*ogl


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


!do i=1,skp
!   do j=1,skp
!      iv1(i)=iv1(i)+alb(i,j)*(st(j)-st2(j))
!   end do
!end do

!do i=1,skp
!   do j=1,skp
!      iv2(i)=iv2(i)+akb(i,j)*iv1(j)
!   end do
!end do

!do i=1,skp
!   i3=i3+(st(i)-st1(i))*iv2(i)
!end do

im1=matmul(aklib,alb)
im2=matmul(akb,im1)
im1=matmul(massb,im2)
im2=im1

do i=1,skp
   i2=i2+im2(i,i)
end do


tk=two*ov*(two*i3+i2)

!*******************************
!derivative of T with respect to sk --numerical testing very good

al=matmul(mass,al)
aklp=matmul(massb,alb)
aklp=matmul(akb,aklp)
!call tab(aklp,6,6,6,6)
aklpt=transpose(aklp)
aklps=matmul(ak,al)
!call tab(aklps,2,2,2,2)
aklpst=transpose(aklps)

iv1=zero
iv2=zero

im1=matmul(aklib,akb)
im2=matmul(aklp,im1)
im2=transpose(im2)

do i=1,skp
   do j=1,skp
      iv2(i)=iv2(i)+im2(i,j)*(st(j)-st1(j))
   end do
end do

im1=matmul(aklib,akb)
im2=matmul(aklpt,im1)
im3=im2-aklpt
im3=transpose(im3)

do i=1,skp
   do j=1,skp
      iv1(i)=iv1(i)+im3(i,j)*(st(j)-st2(j))
   end do
end do

kgs=two*(two*i3+i2)*ogs + two*two*ov*iv2 + two*two*ov*iv1 

!***************************************
!derivative of K wrt to vech(Lk')

kgl=(two*i3+i2)*ogl

ism1=matmul(akli,lk)
ism2=matmul(al,ism1)
ism3=matmul(ak,ism2)
ism1=matmul(akli,ism3)
call vech(ism1,ne,vech1)

kgl=kgl-ov*six*vech1

ism1=matmul(akli,lk)
ism2=matmul(al,ism1)
call vech(ism2,ne,vech1)

kgl=kgl+ov*six*vech1

!****** up to here is trace part--tests perfect!!! ; other part below

k=1
do i=1,ne
   do j=1,3
      dot1(j,i)=st(k)        !for s
      dot1e(j,i)=st1(k)      !for sk
      k=k+1
   end do
end do

dot2=transpose(dot1)
dot2e=dot2                   !for s
sds=matmul(dot2,dot1)
sdsk=matmul(dot2,dot1e)
dot2=transpose(dot1e)
skdsk=matmul(dot2,dot1e)

k=1
do i=1,ne
   do j=1,3
      dot1(j,i)=st2(k)        !for sl
      k=k+1
   end do
end do

dot2=transpose(dot1)          !for sl
sldsk=matmul(dot2,dot1e)
dot1=transpose(dot2e)
slds=matmul(dot2,dot1)

ism1=aklps+aklpst
ism1=matmul(akli,ism1)
ism1=al-ism1
ism1=matmul(ism1,sds)
kcm=ism1

ism1=aklps+aklpst
ism1=matmul(akli,ism1)
ism1=ism1-al
ism1=matmul(ism1,sdsk)
kcm=kcm+ism1


sdsk=transpose(sdsk)
ism1=matmul(akli,aklpst)
ism1=matmul(ism1,sdsk)
kcm=kcm+ism1

ism1=matmul(akli,aklps)
ism1=ism1-al
ism1=matmul(ism1,slds)
kcm=kcm+ism1

ism1=matmul(akli,aklps)
ism1=al-ism1
ism1=matmul(ism1,sldsk)
kcm=kcm+ism1


ism1=matmul(akli,aklpst)
ism1=matmul(ism1,skdsk)
kcm=kcm-ism1

ism1=matmul(kcm,lk)
call vech(ism1,ne,vech1)
kgl=kgl+two*ov*vech1
kcm=transpose(kcm)
ism1=matmul(kcm,lk)
call vech(ism1,ne,vech1)
kgl=kgl+two*ov*vech1
kgl=two*kgl

!********************************************
! electron replsion -- associated gradients

er=zero
ert=zero
ergst=zero
erglt=zero

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

!!$call tab(jij,2,2,2,2)
!!$stop

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

!***************************************
! derivative of er wrt sk -- check numerically good --4/21/00
!  needs i2 and i3 from er above

iv3=zero

!*** form product aklib*jijb*aklib*akb and store in jijb ********
jijb=matmul(akb,jijb)
do i=1,skp
   do j=1,skp
      iv3(i)=iv3(i)+jijb(i,j)*ee(j)
   end do
end do
ergs=(ncv(pp+1)*ncv(ll+1)*two*f0(i1)/dsqrt(i2*pi))*ogs + ncv(pp+1)*ncv(ll+1)*four*ov*df0(i1)*iv3/dsqrt(pi*i2*i2*i2)
ergst=ergst+ergs

!*****************************************
!derivative of er wrt vech(Lk)
! used i2 and i3 from er above 

ism1=matmul(jij,lk)
call vech(ism1,ne,vech1)
ergl=(two*f0(i1)/dsqrt(i2*pi))*ogl + (two*ov*f0(i1)/dsqrt(pi*i2*i2*i2))*vech1

!up to here tests out numerically, but does not include f0 part

ergl=ergl+four*ov*df0(i1)*i3*vech1/dsqrt(pi*i2*i2*i2*i2*i2)

!make m2
k=1
do i=1,ne
   do j=1,3
      dot1(j,i)=ee(k)
      k=k+1
   end do
end do

dot2=transpose(dot1)
ism1=matmul(dot2,dot1)
ism2=matmul(akli,ism1)
ism1=matmul(ism2,jij)

ism2=matmul(ism1,lk)
call vech(ism2,ne,vech1)
ergl=ergl-four*ov*df0(i1)*vech1/dsqrt(pi*i2*i2*i2)

ism2=transpose(ism1)
ism1=matmul(ism2,lk)
call vech(ism1,ne,vech1)
ergl=ergl-four*ov*df0(i1)*vech1/dsqrt(pi*i2*i2*i2)



!!$!make m1
k=1
do i=1,ne
   do j=1,3
      dot1e(j,i)=st1(k)
      k=k+1
   end do
end do

dot2=transpose(dot1e)
ism1=matmul(dot2,dot1)
ism2=matmul(ism1,jij)

ism1=matmul(ism2,lk)
call vech(ism1,ne,vech1)
ergl=ergl+four*ov*df0(i1)*vech1/dsqrt(pi*i2*i2*i2)
!!$
ism1=transpose(ism2)
ism2=matmul(ism1,lk)
call vech(ism2,ne,vech1)
ergl=ergl+four*ov*df0(i1)*vech1/dsqrt(pi*i2*i2*i2)
erglt=erglt+ncv(pp+1)*ncv(ll+1)*ergl

end if
end do
end do

!**************************************
!the stark type term

zht=zero
do i=1,ne
	zht=zht+ncv(i+1)*st(3*i)
end do
zht=zht*ov*efs

zgs=zero
im1=matmul(akb,aklib)
do i=1,ne
	do j=1,skp
	zgs(j)=zgs(j)+ncv(i+1)*(ogs(j)*st(3*i)+ov*im1(j,3*i))*efs
	end do			
end do

zgl=zero
do kk=1,ne
voo=zero
voo(3*kk)=one
k=1
do i=1,ne
   do j=1,3
      dot1(j,i)=st1(k)-st(k)        !for s
      dot1e(j,i)=voo(k)      !for sk
      k=k+1
   end do
end do

dot2=transpose(dot1)
ism1=matmul(dot2,dot1e)
ism2=matmul(ism1,akli)
ism2=matmul(ism2,lk)
call vech(ism2,ne,vech1)

ism1=transpose(ism1)
ism2=matmul(akli,ism1)
ism2=matmul(ism2,lk)
call vech(ism2,ne,vech2)

zgl=zgl+ncv(kk+1)*efs*(st(3*kk)*ogl+ov*(vech1+vech2))


end do

!**************************************
!       gradient stuff
!                      works??????
igvl=kgl+erglt-zgl
igvs=kgs+ergst-zgs

!!$igvl=kgl+erglt-naglt+(dis)*ogl
!!$igvs=kgs+ergst-nagst+(dis)*ogs

do i=1,svh
   tgvl(i)=igvl(i)
   tgvls(i)=ogl(i)
end do
kk=1
do i=(svh+1),(svh+skp)
   tgvl(i)=igvs(kk)
   tgvls(i)=ogs(kk)
   kk=kk+1
end do

if(ii==jj)then
   tgvl=two*tgvl
   tgvls=two*tgvls
end if
   
!!$!*********** test printing of gradients
!!$!            can be removed after testing
!!$write(*,*)'ogs'
!!$write(*,*)ogs
!!$write(*,*)''
!!$write(*,*)'ergst'
!!$write(*,*)ergst
!!$write(*,*)''
!!$write(*,*)'nagst'
!!$write(*,*)nagst
!!$write(*,*)''
!!$write(*,*)'kgs'
!!$write(*,*)kgs
!!$write(*,*)''
!!$write(*,*)'ogl'
!!$write(*,*)ogl
!!$write(*,*)''
!!$write(*,*)'erglt'
!!$write(*,*)erglt
!!$write(*,*)''
!!$write(*,*)'naglt'
!!$write(*,*)naglt
!!$write(*,*)''
!!$write(*,*)'kgl'
!!$write(*,*)kgl
!!$write(*,*)''

END SUBROUTINE gradel
