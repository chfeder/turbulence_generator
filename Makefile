#
#  Makefile for turbulence forcing generator.
#
#  Written by
#
#          Christoph Federrath
#

# C++ compiler (e.g., g++), flags, and HDF5 library path
CCOMP = mpicxx
CFLAGS = -O3

# binary target
BIN = TurbGenDemo

$(BIN) : $(BIN).o
	$(CCOMP) $(CFLAGS) -o $@ $(BIN).o

.SUFFIXES: .cpp .h

.cpp.o:
	$(CCOMP) $(CFLAGS) -c $*.cpp

clean :
	rm -f *.o *~ $(BIN)

# dependencies
$(BIN).o : TurbGen.h
