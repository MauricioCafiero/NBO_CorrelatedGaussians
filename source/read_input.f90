SUBROUTINE read_input

USE matrices_to_be_shared
implicit none
include 'mpif.h'

CHARACTER (LEn=80) :: CHARSTRING,cs,cs2
integer::i,j,k,jj,qq
double precision::disi


if(rank==0)then

   open(unit=12, file="oldquit.dat", status="old", action="readwrite",   &
        position="rewind", iostat=jj)
   if(jj .ne. 0)then
      call biout('Could not open restart file!')
   end if
   read(12,*)M
   read(12,*)ne
   read(12,*)nst
   read(12,*)efs
   
end if

call MPI_BCAST(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(efs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


skp=3*ne
svh=ne*(ne+1)/2
siv=m*((ne*(ne+1)/2)+3*ne)




!***********************************************************************************************

ALLOCATE(ncv(ne+1),scvv(nst),mass(ne,ne),ki(M,M),po(M,M),   &
         hh(M,M),ss(M,M),   &
         C(M),S(M,M),H(M,M),grad(m*(svh+skp)),   &
         dss(m*(svh+skp),m),dhh(m*(svh+skp),m),  &
  			CO(M),iv(siv),sml(ne,ne,nst),sms(skp,skp,nst),	&
  			DS(M*(svh+skp),M),DH(M*(svh+skp),M),Kinetic(M,M),Potential(M,M))



if(rank==0)then
   read(12,*)( ( mass(i,j), j=1,ne), i=1,ne)
end if

call MPI_BCAST(mass,ne*ne,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


if(rank==0)then
   do i=1,ne+1
      read(12,*)ncv(i)
   end do
end if

call MPI_BCAST(ncv,ne+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)




if(rank==0)then
do i=1,siv
   read(12,*)iv(i)
end do
end if

call MPI_BCAST(iv,siv,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


if(rank==0)then
   do i=1,m
      read(12,*)C(i)
   end do
end if

call MPI_BCAST(c,m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

if(rank==0)then
do i=1,nst
   read(12,*)scvv(i)
end do
end if

call MPI_BCAST(scvv,nst,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!write(*,*)scvv
scvv=scvv/nst
!write(*,*)scvv
!stop

if(rank==0)then
do k=1,nst
   !write(*,*)k
   read(12,*)( ( sml(i,j,k), j=1,ne), i=1,ne)
end do
end if

call MPI_BCAST(sml,ne*ne*nst,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


if(rank==0)then
do k=1,nst
   read(12,*)( ( sms(i,j,k), j=1,skp), i=1,skp)
end do
end if

call MPI_BCAST(sms,skp*skp*nst,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


!write(*,*)spm
if(rank==0)then
close(12)
end if



END SUBROUTINE read_input
