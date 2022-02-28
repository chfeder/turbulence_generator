/**********************************************************
 *** Basic example code for using Turbulence Generator. ***
 *** Written by Christoph Federrath, 2022.              ***
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "turbulence_generator.h"

int main(void)
{
  // turbulence_generator.h defines global data structure named 'tgd' of type 'TurbGenData'
  tgd.verbose = false; // in case we want additional output during operations of the turbulence generator

  // First, initialise the turbulence generator with input parameters from 'parameter_file'.
  // If MPI is present, we pass the MPI rank to TurbGen_init_turbulence_generator;
  // note that the MPI rank is simply to control printf to stdout, i.e., only rank 0 will print.
  // This function should only be called once, i.e., to initialise the turbulence generator.
  char * parameter_file = "turbulence_generator.inp";
  int MPI_rank = 0;
  TurbGen_init_turbulence_generator(parameter_file, MPI_rank);

  // as an example, let's loop through a few time steps (e.g., of a hydro code)
  // to demonstrate how the turbulence generator is updated and how the physical turbulent
  // acceleration field is returned by "turbulence_generator.h

  // get the turbulent turnover time
  double t_turb = tgd.decay;
  // let's go through 5 turnover times based on the inputs in 'parameter_file' above
  double t_end = 5.0 * t_turb;

  // set some number of time steps that the hydro code would go through
  int nsteps = 123;
  // time step (this will normally be set by the hydro code;
  // ultimately, one simply feeds in the 'time' at which one wants the respective
  // physical turbulent acceleration field -- see below in call to TurbGen_check_for_update
  double dt = t_end / nsteps;

  // loop through time steps
  for (int step = 0; step <= nsteps; step++) {

    double time = step * dt;
    printf("turbulence_generator: time = %f\n", time);

    // Function to update the turbulence driving mode coefficients.
    // Based on input 'time', it checks if the pattern needs to be updated.
    // If it was updated, return true in 'have_updated_pattern', else return false.
    // Note that this call to TurbGen_update(time) is very cheap and returns right away if current input
    // time means that no update of the turbulence driving pattern is necessary; therefore, we can safely
    // call this in every hydro timestep without wasting compute time.
    bool have_updated_pattern = TurbGen_check_for_update(time);

    if (have_updated_pattern) {
      // say we want the turbulent vector field value at some physical position pos[3]
      double pos[3] = {0.1, -0.5, 0.3};
      double v[3]; // return vector
      TurbGen_get_turb_vector(pos, v); // get the vector
      printf("turbulence_generator: turbulent vector at position (%f %f %f) is vx, vy, vz = %f %f %f\n",
              pos[0], pos[1], pos[2], v[0], v[1], v[2]);
    }

  }

  printf("turbulence_generator: ==========================================================\n");

  // The above example returns the turbulent vector field at only a single position.
  // Here we want the turbulent vector field at all positions on a 3D uniform grid of size nx, nz, nz.
  // In grid codes, this would be function to call, because it is faster than calling TurbGen_get_turb_vector
  // at every requested position; here, it is done for an entire block (uniform grid) of data.

  // random example of physical coordinates of a target output grid (pos_beg must be < pos_end):
  double pos_beg[3] = {0.1, -0.5, 0.3}; // first cell coordinate in uniform grid (x,y,z)
  double pos_end[3] = {0.4, -0.3, 0.5}; // last  cell coordinate in uniform grid (x,y,z)
  int n[3] = {10, 20, 13}; // some random number of cells that define the output grid size
  float * grid_out[3]; // output grid (three cartesian components: vx, vy, vz), which receives the turbulent acceleration field
  // allocate
  grid_out[0] = malloc(n[0] * sizeof(float));
  grid_out[1] = malloc(n[1] * sizeof(float));
  grid_out[2] = malloc(n[2] * sizeof(float));
  // call to return 3 uniform grids with the 3 components of turbulent acceleration field at requested positions
  TurbGen_get_turb_vector_unigrid(pos_beg, pos_end, n, grid_out);
  // indices in x,y,z, i.e., i,j,k, and 1D index
  int i, j, k; long index;
  i = 3; j = 0; k = 10;
  index = k*n[1]*n[0] + j*n[0] + i;
  printf("turbulence_generator: e.g., x-component of turbulent vector field at index i,j,k = (%i %i %i) is %f\n", i,j,k, grid_out[0][index]);
  i = 9; j = 1; k = 4;
  index = k*n[1]*n[0] + j*n[0] + i;
  printf("turbulence_generator: e.g., z-component of turbulent vector field at index i,j,k = (%i %i %i) is %f\n", i,j,k, grid_out[2][index]);
  // clean up
  free(grid_out[0]); free(grid_out[1]); free(grid_out[2]);

  return 0;
}


