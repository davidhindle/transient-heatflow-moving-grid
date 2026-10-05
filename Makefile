#! /bin/csh -f

# fortran pour compiler.

MYFOR = gfortran

# linker also fortran.
LD = ${MYFOR}

# flags for compilation
CFLAGS = -O3 -C -fdollar-ok -fbounds-check -fsignaling-nans -ffpe-trap=invalid,zero,overflow

# debug flag
DFLAGS = -g

# Linker flags 
LDFLAGS =


PROG = pdt.exe

OBJECTS = module_tdma.o pdtransient.o

$(PROG): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS)

module_tdma.o: src/module_tdma.f90
	$(MYFOR) -c $(CFLAGS) $(DFLAGS) $< -o $@

pdtransient.o: src/pdtransient.f90 module_tdma.o
	$(MYFOR) -c $(CFLAGS) $(DFLAGS) $< -o $@

clean:
	rm -f *.o *.mod $(PROG)
