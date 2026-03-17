core1 : necg.o biout.o matrices_to_be_shared.o \
	make_grad.o gradel.o tn.o \
	vech.o outprod.o derror.o \
	matel.o inv.o rrms.o \
	error_function_approximation.o \
	get_current_date_and_time.o matrix_print.o read_input.o 
	mpif77 -o necg necg.o biout.o matrices_to_be_shared.o \
	make_grad.o gradel.o tn.o \
	vech.o outprod.o derror.o \
	matel.o inv.o rrms.o \
	error_function_approximation.o \
	get_current_date_and_time.o matrix_print.o read_input.o   -llapack -lblas

necg.o : necg.f90
	mpif77 -c -O3 necg.f90

rrms.o : rrms.f90
	mpif77 -c -O3 rrms.f90

tn.o : tn.f
	mpif77 -c -O3 tn.f


make_grad.o : make_grad.f90
	mpif77 -c -O3 make_grad.f90

gradel.o : gradel.f90
	mpif77 -c gradel.f90

#uncmnd.o : uncmnd.f
#	mpif77 -c -O3 uncmnd.f

derror.o : derror.f90
	mpif77 -c -O3 derror.f90

vech.o : vech.f90
	mpif77 -c -O3 vech.f90

outprod.o : outprod.f90
	mpif77 -c -O3 outprod.f90

inv.o : inv.f90
	mpif77 -c -O3 inv.f90

matel.o : matel.f90
	mpif77 -c -O3 matel.f90

read_input.o : read_input.f90
	mpif77 -c -O3 read_input.f90

matrices_to_be_shared.o :  matrices_to_be_shared.f90
	mpif77 -c -O3  matrices_to_be_shared.f90


biout.o : biout.f90
	mpif77 -c -O3 biout.f90

get_current_date_and_time.o : get_current_date_and_time.f90
	mpif77 -c -O3 get_current_date_and_time.f90


error_function_approximation.o : error_function_approximation.f90
	mpif77 -c -O3 error_function_approximation.f90 

matrix_print.o : matrix_print.f90
	mpif77 -c -O3 matrix_print.f90

