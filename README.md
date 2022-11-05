# TurbGen - Turbulence Generator #

This tool generates a time sequence of random Fourier modes via an Ornstein-Uhlenbeck (OU) process, used to drive turbulence in hydrodynamical simulation codes. It can also generate single turbulent realisations.

Turbulence driving based on this method is currently supported by implementations in AREPO, FLASH, GADGET, PHANTOM, PLUTO.

Main files:

* 'TurbGen.h' contains the main C++ class with functions and data structures used by the generator.
* 'TurbGen.cpp' is an MPI-parallelised program that computes turbulent field(s) with specified parameters and writes the field(s) to an HDF5 file.
* 'TurbGenDemo.cpp' contains 3 basic examples for how to include and use the generator, including the generation of driving via an OU process and the generation of single turbulent fields.
* 'turbulence_generator.inp' is the parameter file that controls the turbulence driving.

Directory 'plugins' contains examples of how TurbGen.h is used in hydrodynamical codes.

### What are the first steps when working with this repository? ###

* Create a fork of this repo to your own bitbucket account.
* Clone the fork from your bitbucket account to your local computer.
* Make your own modifications, commit, push, etc.
* If you want (some of) your own changes to go into this main repo, please create a pull request.

### Copyright notes ###

If you use this code, please refer to and cite Federrath et al. (2010, A&A 512, A81):
https://ui.adsabs.harvard.edu/abs/2010A%26A...512A..81F/abstract

@ARTICLE{2010A&A...512A..81F,<br/>
author = {{Federrath}, C. and {Roman-Duval}, J. and {Klessen}, R.~S. and {Schmidt}, W. and {Mac Low}, M. -M.},<br/>
title = "{Comparing the statistics of interstellar turbulence in simulations and observations. Solenoidal versus compressive turbulence forcing}",<br/>
journal = {\aap},<br/>
keywords = {hydrodynamics, ISM: clouds, ISM: kinematics and dynamics, methods: numerical, methods: statistical, turbulence, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},<br/>
year = 2010,<br/>
month = mar,<br/>
volume = {512},<br/>
eid = {A81},<br/>
pages = {A81},<br/>
doi = {10.1051/0004-6361/200912437},<br/>
archivePrefix = {arXiv},<br/>
eprint = {0905.1060},<br/>
primaryClass = {astro-ph.SR},<br/>
adsurl = {https://ui.adsabs.harvard.edu/abs/2010A&A...512A..81F },<br/>
adsnote = {Provided by the SAO/NASA Astrophysics Data System}}

### Main contact ###

* Christoph Federrath (christoph.federrath@anu.edu.au)
