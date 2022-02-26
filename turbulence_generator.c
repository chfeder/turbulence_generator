
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "turbulence_generator.h"

int main(void)
{
  //TurbGen_generate_driving_file("forcing_generator.inp");
  //double time_in_file;
  //TurbGen_read_driving_file("turb_pl+2.50_v1.00e+00_zeta0.5_seed140281.dat", 2.5, &time_in_file);
  //printf("time_in_file = %f\n", time_in_file);

  if (TurbGen_init_turbulence_generator("forcing_generator.inp") != 0) return -1;

  TurbGen_print_info();

  TurbGen_update(0.1);

  double vx, vy, vz;
  TurbGen_get_turb_vector(0.0,   0.1, 0.1, -0.2, &vx, &vy, &vz);
  printf("vx, vy, vz = %f %f %f\n", vx, vy, vz);

  return 0;
}


