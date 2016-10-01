debug=false
mpi=false

ifeq ($(mpi),true)
CC := mpif90
#mpi_mod := set_mpi.o create_mpi_types.o
#main := veff_self_mpi_prova
main := veff_self_new_mpi
mpi_mod := set_mpi.o  utilities_mpi.o
vegas_mod := vegas_module_mpi.o
else
CC := gfortran
main := veff_self_new
mpi_mod :=  utilities.o
vegas_mod := vegas_module.o
endif

ifeq ($(debug),true)
OPT   := -g
FLAGS := -pg -fbounds-check -ffree-line-length-172  -Wall -Wextra  -fno-backtrace -finit-real=nan\
	--pedantic   -ffpe-trap=invalid,zero,overflow
else
OPT   := -pg -O3 #-Ofast
FLAGS := -fbounds-check -ffree-line-length-172 
endif



LIB = -lcuba

.SUFFIXES := # delete default known suff. for implicit rules
.SUFFIXES := .x .o .f90 .F90 .mod

HF := veff_self_new_HF





mod := precision.o time.o typeveffdata_new.o algebra.o polarisation.o representation.o omffst6.o der_er_new.o \
	integration.o  selfenergy.o Boundary.o selfenergy_first.o $(vegas_mod) selfenergy_second.o $(mpi_mod)

#mod := precision.o time.o typeveffdata.o algebra.o polarisation.o representation.o omffst6.o utilities.o \
	der_er_new.o  integration.o  selfenergy_save.o Boundary.o $(mpi_mod)	
#vpath %.c
#vpath %.o exec
#vpath %.mod exec
#vpath %  exec

all: clean $(main) $(HF)

$(main): $(mod) $(main).o
	$(CC) $(OPT) $(FLAGS) $(LIB) $(mod)  $(main).o -o $(main).x


$(HF): $(mod) $(HF).o
	$(CC) $(OPT) $(FLAGS) $(LIB) $(mod)  $(HF).o -o $(HF).x


%.o: %.f90
	$(CC) $(OPT) $(FLAGS)  $(LIB) -c $<

%.o: %.f
	$(CC) $(OPT) -c $<

%.o: %.F90
	$(CC) $(OPT)  $(FLAGS) $(LIB) -c $<


#prova:  matrices.o typeveffdata.o selfenergy.f90
#	$(CC) $(OPT)  matrices.o typeveffdata.o -c  selfenergy.f90 


.PHONY : clean deepclean
clean:
	rm -f *.x *.o *~ *# *.mod $(main) *.out*
deepclean:
	rm -r *.dSYM
