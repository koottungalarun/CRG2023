#FC = ifort
#FFLAGS = -O3 -diag-disable vec -r8 -free -Tf

FC = gfortran
FFLAGS = -O3 -fdefault-real-8  -g -Wall -fcheck=all
LDLIBS  = -L/usr/lib -lumfpack -lamd -lm -lblas -llapack -lcholmod -lcolamd
LDLIBS  += -L/usr/local/lib

TARGETS = euler

ALL: $(TARGETS)

SRC = $(wildcard *.f95)
OBJ = $(patsubst %.f95,%.o,$(SRC))

euler: $(OBJ)
	$(FC) -o euler $(OBJ)	$(LDLIBS)

comvar.mod: comvar.f95
	$(FC) -c $(FFLAGS) comvar.f95

mUMFPACK.mod: umfpack.f95
	$(FC) -c $(FFLAGS) umfpack.f95

%.o: %.f95 comvar.mod mUMFPACK.mod
	$(FC) -c $(FFLAGS) $*.f95 	

clean:
	rm -f $(OBJ) $(TARGETS) *.mod *.plt
