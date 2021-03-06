CFLAGS	         = 
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       = 
FCC=mpicxx
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


DEPENDENCIES_grid = overlord_gsl.cpp utils.cpp class_grid.cpp orb.cpp sphere_lebedev_rule.cpp sasha.cpp
DEPENDENCIES_grid.o = overlord_gsl.o utils.o class_grid.o orb.o sphere_lebedev_rule.o sasha.o

overlord_gsl: 
	-${FCC} -O2 -c ${DEPENDENCIES_grid} ${flags_gsl}
	-${FCC} -o overlord_gsl ${DEPENDENCIES_grid.o} ${flags_gsl}


all:
	make Hatom
	make overlord_gsl
	$(RM) *.o
