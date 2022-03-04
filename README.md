# Turbulence Generator #

This tool generates a time sequence of random Fourier modes, used to drive turbulence in hydrodynamical simulation codes. It can also generate single turbulent realisations.

Turbulence driving based on this method is currently supported by implementations in AREPO, FLASH, GADGET, PHANTOM, PLUTO, QUOKKA.

Main files:

* 'TurbGen.h' contains the main C++ class with functions and data structures used by the generator.
* 'TurbGenDemo.C' contains 3 basic examples for how to include and use the generator.
* 'turbulence_generator.inp' is the parameter file that controls the turbulence driving.

Directory 'plugins' contains some examples of how TurbGen.h is used in hydrodynamical codes.

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
