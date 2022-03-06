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

BINS = TurbGen TurbGenDemo

all: $(BINS)

clean:
	rm -f *.o *~ $(BINS)

$(BINS): %: %.cpp
	$(CCOMP) $(CFLAGS) -o $@ $<

# dependencies
$(BINS).o: TurbGen.h
