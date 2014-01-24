
FOR = gfortran
RM = rm
LN = ln
FLIBS = -llapack -lblas
FFLAGS = -c -g -Wall -m64 -fopenmp -fPIC 
FFLAGSOPT = -O3 -fopenmp


all: main.f90
	$(FOR)  main.f90 -o atomic -L../src -ldft

clean:
