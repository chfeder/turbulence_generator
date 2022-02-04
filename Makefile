#
#  Makefile for forcing generator.
#
#  Written by
#
#          Christoph Federrath
#
# The easiest way to compile this is to './setup unitTest/StirFromFile -auto'
# in the main FLASH directory, so the Makefile.flash is linked to the correct site file
# for the machine that you are compiling on. After that, simply 'make' in this directory
# should compile the forcing generator (key here is that it needs to be compiled with the real -> real*8 flag).
# Finally, run the executable 'forcing_generator' to produce the forcing file. The forcing file then
# needs to be linked or copied into the flash4 run directory.

# Setting the compiler 'by hand':
# For gnu compiler:
# FCOMP = mpif90 -fdefault-real-8 -O3
# FFLAGS_OPT = -c -fdefault-real-8 -O3
# LFLAGS_OPT = -fdefault-real-8 -o
# For Intel compiler:
# FCOMP = mpif90 -r8 -i4 -O3
# FFLAGS_OPT = -c -r8 -O3
# LFLAGS_OPT = -r8 -o

# Automatically define the compiler settings based on the FLASH ./setup selected site Makefile
# Note that 'Makefile.flash' is automatically generated as a link when the ./setup script is run

FLASHMAKEFILE = Makefile.flash

ifeq ($(wildcard $(FLASHMAKEFILE)),)
FLASHMAKEFILE_MISSING=true
else
$(eval $(shell grep 'FCOMP ' $(FLASHMAKEFILE)))
$(eval $(shell grep 'FFLAGS_OPT ' $(FLASHMAKEFILE)))
$(eval $(shell grep 'LFLAGS_OPT ' $(FLASHMAKEFILE)))
endif

ifeq ($(FLASHMAKEFILE_MISSING),true)
$(info *** $(FLASHMAKEFILE) not present; run ./setup ... for your FLASH setup including the forcing module, and try again ***)
all:
else
all: forcing_generator
endif

forcing_generator : forcing_generator.o
	$(FCOMP) $(LFLAGS_OPT) $@ forcing_generator.o

.SUFFIXES: .F90

.F90.o:
	$(FCOMP) $(FFLAGS_OPT) $*.F90

clean :
	rm -f *.o *.mod *~ forcing_generator

test :
	@echo FCOMP = $(FCOMP)
	@echo FFLAGS_OPT = $(FFLAGS_OPT)
	@echo LFLAGS_OPT = $(LFLAGS_OPT)
	@echo FLASHMAKEFILE_MISSING = $(FLASHMAKEFILE_MISSING)
