MODULE  matrices_to_be_shared

  !double precision, dimension(6)::ogs,kgs

character(len=32)::hostname
integer::flag, ierr, rank, nproc
  

double precision, dimension(:,:), allocatable :: ngv,mass

double precision, dimension(:,:,:), allocatable :: sml,sms	

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::	&
  						S,H,CSTORE,TSTORE1,TSTORE2,	&
  						TSTORE3,RZ1STORE,RZ2STORE,	&
  						Kinetic,Potential,ki,po,hh,ss,			&
  						OVV1,OVV2,OVV3,OVV4
  													

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::	&
  						DS,DH,DSS,GRADM,dhh
  
  DOUBLE PRECISION, DIMENSION(2,2) :: Q,Y

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::		&
  						iv,ncv,scvv,			&	!alpha1,...
  						C,		&	!alpha2,...,correlation,linear coeff
  						GRAD,CO,			&	!gradient vector,another lin coeff
  						GRADA,GRADC	!storage for cholesky params & grad vec
  								


DOUBLE PRECISION, PARAMETER :: ZERO    = 0.00D+00		!REAL NUMERICAL CONSTANTS...
  DOUBLE PRECISION, PARAMETER :: HALF    = 0.50D+00
  DOUBLE PRECISION, PARAMETER :: QUARTER = 0.25D+00
  DOUBLE PRECISION, PARAMETER :: EIGHTH  = QUARTER*HALF
  DOUBLE PRECISION, PARAMETER :: ONE     = 1.00D+00
  DOUBLE PRECISION, PARAMETER :: TWO     = 2.00D+00
  DOUBLE PRECISION, PARAMETER :: THREE   = 3.00D+00
  DOUBLE PRECISION, PARAMETER :: FOUR    = 4.00D+00
  DOUBLE PRECISION, PARAMETER :: FIVE    = 5.00D+00
  DOUBLE PRECISION, PARAMETER :: SIX     = 6.00D+00
  DOUBLE PRECISION, PARAMETER :: SEVEN   = 7.00D+00
  DOUBLE PRECISION, PARAMETER :: EIGHT   = 8.00D+00
  DOUBLE PRECISION, PARAMETER :: SIXTEEN = 16.00D+00

  DOUBLE PRECISION, PARAMETER :: C1  = 0.564691197D-04	!REAL PADE' APPROXIMATION CONSTANTS...
  DOUBLE PRECISION, PARAMETER :: C2  = 0.758433197D-03	!	...for the error fnc approximation...
  DOUBLE PRECISION, PARAMETER :: C3  = 0.769838037D-02
  DOUBLE PRECISION, PARAMETER :: C4  = 0.629344460D-01
  DOUBLE PRECISION, PARAMETER :: C5  = 0.213271302D+00
  DOUBLE PRECISION, PARAMETER :: C6  = 0.720266520D-04
  DOUBLE PRECISION, PARAMETER :: C7  = 0.955528842D-03
  DOUBLE PRECISION, PARAMETER :: C8  = 0.101431553D-01
  DOUBLE PRECISION, PARAMETER :: C9  = 0.738522953D-01
  DOUBLE PRECISION, PARAMETER :: C10 = 0.338450368D+00
  DOUBLE PRECISION, PARAMETER :: C11 = 0.879937801D+00


  DOUBLE PRECISION, PARAMETER :: PI=3.1415926535898D+00	!f(circle) -> (circumference/diameter)

  INTEGER, PARAMETER :: NMAX = 5				!NUMBER OF NONLINEAR PARAMETERS

  DOUBLE PRECISION :: dis,R_ENER,EMIN,SQNORMGRAD,efs			!BOND DISTANCE/ENERGIES/|GRADIENT|**2
  

      DOUBLE PRECISION, PARAMETER :: EZ1=  0.0d+00
  DOUBLE PRECISION, PARAMETER :: Ex1=  0.0d+00
  DOUBLE PRECISION, PARAMETER :: EY1=  0.952627944162921D+0

DOUBLE PRECISION, PARAMETER ::      EZ2=  0.0d+00
   DOUBLE PRECISION, PARAMETER ::   EZ3=  0.0d+00
     DOUBLE PRECISION, PARAMETER :: EX2=  0.825000000000027D+0 
     DOUBLE PRECISION, PARAMETER :: EX3=  -0.825000000000027D+0 
  
     DOUBLE PRECISION, PARAMETER :: EY2=  -0.476313972081468D+0 
     DOUBLE PRECISION, PARAMETER :: EY3=  -0.476313972081468D+0


  DOUBLE PRECISION :: ELAPSED_TIME


!***** Initialize some integers:
  INTEGER ::	M, jjj,NN,NE,SVH,skp,nst,siv				
!
!               m -- nyumer of basis functions, size of hamiltonian matrix
!               nn -- number of nuclei
!               ne -- number of electrons, dimension of Ak, 
!               skp -- 3*ne = size of sk vector, dimension of Kroneker products
!               svh -- 1/2*[n*(n+1)] = size of vech
!               nst -- number of symmetry terms
!               siv -- size of input vector


 CHARACTER(LEN=66) ::	FILENAME,RESTARTNAME,		&
  						RESTARTOUTPUT,OPTIMIZEDWAVE,&
  						METHODNAME,GENNAME

character(len=1)::symcode



END MODULE matrices_to_be_shared
