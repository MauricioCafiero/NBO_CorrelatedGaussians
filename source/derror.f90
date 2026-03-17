DOUBLE PRECISION FUNCTION dF0(Z)

  USE matrices_to_be_shared

  implicit none
  double precision::dF0, Z, f0,anum,aden,danum,daden,r,r3

   

!write(*,*)'hello from f0',z




      IF(Z .LE. 16.3578D+00) THEN





      F0=DSQRT((((((C1 *Z  	&
                  + C2)*Z  	&
                  + C3)*Z  	&
                  + C4)*Z  	&
                  + C5)*Z  	&
                  + ONE)/  	&
             ((((((C6 *Z  	&
                  + C7)*Z  	&
                  + C8)*Z  	&
                  + C9)*Z  	&
                  + C10)*Z  	&
                  + C11)*Z  	&
                  + ONE))

ANUM=C1*(Z**5)+C2*(Z**4)+C3*(Z**3)+C4*(Z**2)+C5*Z+ONE
      ADEN=C6*(Z**6)+C7*(Z**5)+C8*(Z**4)+C9*(Z**3)+C10*(Z**2)+C11*Z+ONE
      DANUM=5.00*C1*(Z**4)+FOUR*C2*(Z**3)+THREE*C3*(Z**2)+TWO*C4*Z+C5
      DADEN=6.00*C6*(Z**5)+5.00*C7*(Z**4)+FOUR*C8*(Z**3)+THREE*C9*(Z**2)+TWO*C10*(Z)+C11
      dF0=(ONE/(TWO*f0))*((ADEN*DANUM-ANUM*DADEN)/ADEN**2)


  	ELSE
  	    R=TWO*Z
  		R3=PI/TWO
  		F0=DSQRT(R3/R)

                dF0= -HALF*(1/DSQRT(PI/(FOUR*Z)))*PI/(FOUR*Z*Z)

     END IF

!write(*,*)'f0',f0


END FUNCTION dF0
