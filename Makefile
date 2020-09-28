include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
include ${SLEPC_DIR}/lib/slepc/conf/slepc_variables

CFLAGS	         = 
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       = 
FCC=mpicc
FFC=ifort

flags_gsl = -I${GSL_DIR}/include -L${GSL_DIR}/lib -lgsl -lgslcblas 

DEPENDENCIES_BOUND = num_mod.f90 quads.f90 gq.f90 dvr_mod.f90 mol_mod.f90 setgrid.f90 spherical.f90 wigner.f90

DEPENDENCIES_BOUND.o = num_mod.o quads.o gq.o dvr_mod.o mol_mod.o setgrid.o spherical.o wigner.o

DEP.o:
	-${FFC} -c ${DEPENDENCIES_BOUND}

Hatom.o: DEP.o
	-${FFC} -c Hatom.f90 -o Hatom.o -mkl
Hatom: Hatom.o
	-${FFC} -o Hatom Hatom.o ${DEPENDENCIES_BOUND} -mkl


DEPENDENCIES_grid = overlord_gsl.cpp utils.cpp class_grid.cpp orb.cpp sphere_lebedev_rule.cpp sasha_mpi.cpp bieleclist.cpp 
DEPENDENCIES_grid.o = overlord_gsl.o utils.o class_grid.o orb.o sphere_lebedev_rule.o sasha_mpi.o bieleclist.o

overlord_gsl: 
	-${CLINKER} -O2 -I${PETSC_DIR}/include${PETSC_SUBDIR} -c ${DEPENDENCIES_grid} ${SLEPC_INCLUDE} ${flags_gsl} -mkl  
	-${CLINKER} -o overlord_gsl ${DEPENDENCIES_grid.o} ${SLEPC_LIB} ${flags_gsl} -mkl 

SOLWRITEROBJS = solution_writer.o
solution_writer: $(SOLWRITEROBJS)
	-${CLINKER} -o solution_writer ${SOLWRITEROBJS} ${PETSC_LIB}
	${RM} ${SOLWRITEROBJS}

all:
	make Hatom
	make overlord_gsl
	make solution_writer
	$(RM) *.o
