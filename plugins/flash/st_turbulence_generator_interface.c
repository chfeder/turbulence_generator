#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mangle_names.h"
#include "st_turbulence_generator.h"

// Initialise the turbulence generator with input parameters from 'parameter_file'.
// If MPI is present, we pass the MPI rank to TurbGen_init_turbulence_generator;
// note that the MPI rank is simply to control printf to stdout, i.e., only rank 0 will print.
// This function should only be called once, i.e., to initialise the turbulence generator.
// Here we also return the delta time between updates of driving patterns, which can be
// used to constrain the simulation timestep (just to make sure the code goes through all
// driving patterns and doesn't skip any; hardly ever happens though, for a real application).
void FTOC(st_init_turbulence_generator_c)(char * parameter_file, double * dt_driv) {
  tgd.verbose = false; // for additional output during operations of the turbulence generator
  int st_turb_gen_MyPE = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &st_turb_gen_MyPE);
  if (TurbGen_init_turbulence_generator(parameter_file, st_turb_gen_MyPE) != 0) exit(-1);
  *dt_driv = tgd.decay / tgd.nsteps_per_turnover_time; // return time between pattern updates
}

// Function to update the turbulence driving mode coefficients.
// Based on input 'time', it checks if the pattern needs to be updated.
// If it was updated, return 1 in 'have_updated_pattern', else return 0.
// Note that this call to TurbGen_update(time) is very cheap and returns right away if current input
// time means that no update of the turbulence driving pattern is necessary; therefore, we can safely
// call this in every hydro timestep without wasting compute time.
void FTOC(st_check_for_update_of_turbulence_pattern_c)(double * time, int * have_updated_pattern) {
  if (TurbGen_check_for_update(*time)) {
    *have_updated_pattern = 1;
  } else {
    *have_updated_pattern = 0;
  }
}

// Function returns the turbulent vector field (vx,vy,vz) given start and end coordinates
// pos_beg and pos_end, on a uniform grid of size n[0]*n[1]*n[2]. Returns single precision (float).
void FTOC(st_get_turb_vector_unigrid_c)(const double pos_beg[3], const double pos_end[3],
                                        const int n[3], float * vx, float * vy, float * vz) {
  float * grid_out[3] = {vx, vy, vz};
  TurbGen_get_turb_vector_unigrid(pos_beg, pos_end, n, grid_out);
}
