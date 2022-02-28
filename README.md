# Turbulence Driving Generator #

This tool generates a time sequence of random Fourier modes, used to drive turbulence in hydrodynamical simulation codes.

Turbulence driving based on this method is currently supported by implementations in AREPO, FLASH, GADGET, PHANTOM, PLUTO, QUOKKA.

Main files:
 - 'turbulence_generator.h' contains the main library functions and data for using the generator.
 - 'turbulence_generator.c' contains a basic example for how to include and use the generator.
 - 'turbulence_generator.inp' is the parameter file that controls the turbulence generator.

Directory 'plugins' contains some examples of how the turbulence_generator.h library is used in hydrodynamical codes.

### What are the first steps when working with this repository? ###

* Create a fork of this repo to your own bitbucket account.
* Clone the fork from your bitbucket account to your local computer.
* Make your own modifications, commit, push, etc.
* If you want (some of) your own changes to go into this main repo, please create a pull request.

### Copyright notes ###

If you use this code, please refer to and cite Federrath et al. (2010, A&A 512, A81):
https://ui.adsabs.harvard.edu/abs/2010A%26A...512A..81F/abstract

### Main contact ###

* Christoph Federrath (christoph.federrath@anu.edu.au)
