!!****if* source/physics/sourceTerms/Stir/StirFromFileMain/Stir_data
!!
!! NAME
!!  Stir_data
!!
!! SYNOPSIS
!!  Stir_data()
!!
!! DESCRIPTION
!!  Stores the local data for Source Term: StirFromFile
!!
!! AUTHOR
!!  Christoph Federrath, 2008-2022
!!
!!***

Module Stir_data

#include "Flash.h"

  character (len=80), save :: st_infilename
  logical, save  :: st_useStir, st_computeDt
  real, save :: st_stop_driving_time
  real, save :: st_x_amplitude_factor, st_y_amplitude_factor, st_z_amplitude_factor
  real(kind=8), save :: dt_update_accel

  ! local container of real kind=4 to receive the turbulent acceleration field
  real(kind=4), save, dimension(NXB, NYB, NZB) :: accx, accy, accz

end Module Stir_data
