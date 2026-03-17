subroutine make_grad(eng,norm)

use matrices_to_be_shared
implicit none
include 'mpif.h'

integer::i,j,k,l,ii,jj,ic,nx,stf
double precision::eng,norm
double precision, dimension(svh+skp)::tgv,tgvs,tgvt,tgvst

dss=zero
dhh=zero
dh=zero
ds=zero

do l=rank+1,m,nproc
   j=1
   do k=1,m

      tgvt=zero
      tgvst=zero
      
      do stf=1,nst
      call gradel(stf,l,k,tgv,tgvs)
      tgvt=tgvt+scvv(stf)*tgv
      tgvst=tgvst+scvv(stf)*tgvs
      end do

      do i=0,(svh+skp-1)
         dh((l-k)*(svh+skp)+j+i,k)=tgvt(i+1)
         ds((l-k)*(svh+skp)+j+i,k)=tgvst(i+1)
      end do
      j=j+(svh+skp)
   end do
end do

!call tab(dh,18,2,18,2)
call MPI_REDUCE(dh,dhh,siv*m,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(ds,dss,siv*m,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,MPI_COMM_WORLD,ierr)
!call tab(dhh,18,2,18,2)
!stop
deallocate(grad)
allocate(gradm(m,m*(svh+skp)),grada(m*(svh+skp)),gradc(m),grad(m*(svh+skp+1)))

if(rank==0)then
ds=dss
dh=dhh

do i=1,m
   do j=1,(m*(svh+skp))
      gradm(i,j)=dh(j,i)-eng*ds(j,i)
      gradm(i,j)=norm*gradm(i,j)
   end do
end do

ii=(svh+skp)
do j=1,m
   do k=1,ii
      ic=(j-1)*ii+k
      grada(ic)=zero
      do i=1,m
         if(i .ne. j)then
            grada(ic)=grada(ic)+two*c(i)*c(j)*gradm(i,ic)
         else
            grada(ic)=grada(ic)+c(i)*c(j)*gradm(i,ic)
         end if
      end do
   end do
end do

gradc=zero

do j=1,m
   do i=1,m
      gradc(i)=gradc(i)+two*norm*(h(i,j)-eng*s(i,j))*c(j)
   end do
end do

nx=m*ii
do i=1,nx
   grad(i)=grada(i)
end do

do i=1,m
   grad(i+nx)=gradc(i)
end do
end if
deallocate(gradm,grada,gradc)


call MPI_BCAST(grad,siv+m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end subroutine make_grad
