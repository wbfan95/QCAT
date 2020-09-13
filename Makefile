FC = ifort
EXE = QCAT
FLAGS = -mkl
# -O2 -static-intel
pot = ohch4.o

objects = main.o readInput.o initMol.o initTraj.o

default : $(objects) $(modules)
	$(FC) $(FLAGS) $(objects) $(modules) $(pot) \ 
	${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \ 
	-o $(EXE)
	time ./$(EXE) ohch4.in | tee stdout.log

clean :
	rm -f *.o *.mod stdout.log

# compile module
module.o : module.f90
	$(FC) $(FLAGS) -c module.f90

log.o : log.f90
	$(FC) $(FLAGS) -c log.f90

atomProp.o : atomProp.f90
	$(FC) $(FLAGS) -c atomProp.f90

surface.o : surface.f90
	$(FC) $(FLAGS) -c surface.f90 -llapack

modules = module.o log.o atomProp.o surface.o

# compile main program
main.o : main.f90 $(modules)
	$(FC) $(FLAGS) -c main.f90

readInput.o : readInput.f90 $(modules)
	$(FC) $(FLAGS) -c readInput.f90
    
initMol.o : initMol.f90 $(modules)
	$(FC) $(FLAGS) -c initMol.f90

initTraj.o : initTraj.f90 $(modules)
	$(FC) $(FLAGS) -c initTraj.f90
