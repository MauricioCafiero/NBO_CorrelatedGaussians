SUBROUTINE inv(Ma,mi)


USE matrices_to_be_shared
implicit none			  

double precision, dimension(2,2)::Ma,MI,test


mi(1,1)=one/(ma(1,1)-(ma(1,2)*ma(1,2)/ma(2,2)))
mi(1,2)=-ma(1,2)/ma(1,1)/(ma(2,2)-(ma(1,2)*ma(1,2)/ma(1,1)))
mi(2,1)=-ma(1,2)/ma(2,2)/(ma(1,1)-(ma(1,2)*ma(1,2)/ma(2,2)))
mi(2,2)=one/(ma(2,2)-(ma(1,2)*ma(1,2)/ma(1,1)))


test=matmul(ma,mi)
!call tab(test,2,2,2,2)

end subroutine inv


