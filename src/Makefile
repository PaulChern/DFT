# Build file for the atomic project

#MPIF90=mpif90
CPP = gfortran
RM = rm
FLIBS = -llapack -lblas
FFLAGS = -m64 -fopenmp
FFLAGSOPT = -O3 -fopenmp

all: atomic.f90 base.f90 grid.f90 linear.f90 quad.f90 consts.f90 main.f90
	$(CPP) $(FFLAGS) $(FLIBS) atomic.f90 base.f90 grid.f90 linear.f90 quad.f90 consts.f90 main.f90 -o atomic  

opt: atomic.f90 base.f90 grid.f90 main.f90
	$(CPP) $(FFLAGSOPT) $(FFLAGS) $(FLIBS) atomic.f90 base.f90 grid.f90 linear.f90 quad.f90 consts.f90 main.f90 -o atomic

clean:
	$(RM) atomic.mod base.mod grid.mod linear.mod quad.mod consts.mod atomic
