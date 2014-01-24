# Build the DFT library

FOR = gfortran
RM = rm
LN = ln
FLIBS = -llapack -lblas
FFLAGS = -c -g -Wall -m64 -fopenmp -fPIC 
FFLAGSOPT = -O3 -fopenmp


all: atomic.f90 base.f90 grid.f90 linear.f90 quad.f90 consts.f90 coupling.f90
	$(FOR) $(FFLAGS) consts.f90 -o consts.o
	$(FOR) $(FFLAGS) base.f90 -o base.o
	$(FOR) $(FFLAGS) grid.f90 -o grid.o
	$(FOR) $(FFLAGS) linear.f90 -o linear.o
	$(FOR) $(FFLAGS) quad.f90 -o quad.o
	$(FOR) $(FFLAGS) atomic.f90 -o atomic.o
	$(FOR) $(FFLAGS) coupling.f90 -o coupling.o
#	$(FOR) -fopenmp -shared -o libdft.so.1 atomic.o base.o grid.o linear.o quad.o consts.o
#	$(FOR) -fopenmp -shared -fpic atomic.f90 base.f90 grid.f90 linear.f90 quad.f90 consts.f90 \
#		-o libdft.so.1
	$(FOR) -fopenmp -shared -fpic -o libdft.so.1 \
		atomic.o base.o grid.o linear.o quad.o consts.o coupling.o
	$(LN) -sf libdft.so.1 libdft.so

clean:
	$(RM) atomic.mod base.mod grid.mod linear.mod quad.mod consts.mod coupling.mod
	$(RM) atomic.o base.o consts.o linear.o quad.o grid.o coupling.o
	$(RM) libdft.so libdft.so.1
