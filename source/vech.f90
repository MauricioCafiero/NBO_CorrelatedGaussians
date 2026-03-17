SUBROUTINE vech(Ma,dim,vec)

implicit none
integer::dim,dim2,i,j,k
double precision, dimension(dim,dim)::Ma      
double precision, dimension(dim*(dim+1)/2)::vec      

!call tab(Ma,dim,dim,dim,dim)

k=1
do j=1,dim
do i=j,dim
   !write(*,*)k,i,j
      vec(k)=Ma(i,j)
   k=k+1
end do
end do
 
!write(*,*)vec
!stop
!call tab(Ma,dim,dim,dim,dim)
   


END SUBROUTINE vech
