!!****if* source/physics/sourceTerms/Stir/StirFromFileMain/Stir
!!
!! NAME
!!  Stir
!!
!! SYNOPSIS
!!  Stir(integer(in) :: blockCount,
!!       integer(in) :: blockList(blockCount),
!!       real(in)    :: dt)
!!
!! DESCRIPTION
!!   Apply the turbulence driving operator on the list of blocks provided as input
!!
!! ARGUMENTS
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to apply turbulence driving
!!   dt           : the current timestep
!!   pass         : optional: dt-pass, in case of split solver
!!
!! AUTHOR
!!   Christoph Federrath, 2008-2022
!!
!!    (Aug 2013: added a write-out of the correction for the center-of-mass motion)
!!    (2012: added a pre-proc statement to treat the acceleration field as a force for testing; not the default)
!!    (2011/2012: added a write-out of the energy injection rate)
!!    (2017: added dynamic allocation of st_nmodes arrays)
!!    (2022: use of turbulence_generator c library/header)
!!
!!***

subroutine Stir(blockCount, blockList, dt, pass)

  use Stir_data
  use Driver_interface, ONLY : Driver_getSimTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface,   ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr, &
                               Grid_getDeltas, Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords
  use Driver_data, ONLY : dr_restart, dr_nstep, dr_globalMe
  use IO_data, ONLY : io_integralFreq

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in)                        :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList
  real, intent(in)                           :: dt
  integer, intent(in), optional              :: pass

  logical, save              :: firstCall=.true., driving_stopped=.false.
  real, save                 :: last_time = HUGE(real(1.0))
  integer                    :: blockID, i, j, k, ii, jj, kk
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical                    :: update_accel = .true.
  integer                    :: update_accel_int, error
  real, dimension(MDIM)      :: del, blockSize, blockCenter ! MDIM is always 3
  real(kind=8)               :: pos_beg(MDIM), pos_end(MDIM)
  integer                    :: ncells(MDIM)
  real                       :: time, ekin_old, ekin_new, d_ekin
  real                       :: mass, xMomentum, yMomentum, zMomentum, xForce, yForce, zForce

  integer, parameter :: funit = 22
  character(len=80)  :: outfile = "stir.dat"
  logical            :: check_for_io

  integer, parameter :: nGlobalSum = 7                  ! Number of globally-summed quantities
  real(kind=8)       :: globalSumQuantities(nGlobalSum) ! Global summed quantities
  real(kind=8)       :: localSumQuantities(nGlobalSum)  ! Global summed quantities
  real(kind=8)       :: ekin_added, ekin_added_red, dvol, dmass, accel

  real, DIMENSION(:,:,:,:), POINTER :: solnData

! this is to determine whether we remove the average force and average momentum in every timestep
#define CORRECT_BULK_MOTION
#ifdef CORRECT_BULK_MOTION
  real(kind=8) :: correction
#else
  real(kind=8), parameter :: correction = 0.0
#endif

  if (firstCall) then
    if (dr_globalMe == MASTER_PE) then
      open(funit, file=trim(outfile), position='APPEND') ! write header
      write(funit,'(10(1X,A16))') '[00]time', '[01]dt', '[02]d(Ekin)', '[03]d(Ekin)/dt', &
                                  '[04]xForce', '[05]yForce', '[06]zForce', &
                                  '[07]xMomentum', '[08]yMomentum', '[09]zMomentum'
      close(funit)
    endif
    firstCall = .false.
  endif

  ! =====================================================================

  ! if not using stirring or stirring is turned off after some time (st_stop_driving_time), then return
  if ((.not. st_useStir) .or. driving_stopped) return

  call Driver_getSimTime(time)

  if (time .ge. st_stop_driving_time) then
    driving_stopped = .true.
    ! clear the acceleration field
    accx(:,:,:) = 0.0
    accy(:,:,:) = 0.0
    accz(:,:,:) = 0.0
#ifdef ACCX_VAR
    do blockID = 1, blockCount
      call Grid_getBlkPtr(blockList(blockID),solnData)
      solnData(ACCX_VAR,:,:,:) = 0.0
      solnData(ACCY_VAR,:,:,:) = 0.0
      solnData(ACCZ_VAR,:,:,:) = 0.0
      call Grid_releaseBlkPtr(blockList(blockID),solnData)
    enddo
#endif
    if (dr_globalMe == MASTER_PE) print *, 'Stir: TURBULENCE DRIVING STOPPED!'
    return
  endif

  call Timers_start("Stir")

  ! check if we need to update the turbulent acceleration field
  update_accel = .false.
  call st_check_for_update_of_turbulence_pattern_c(real(time,kind=8), update_accel_int)
  if (update_accel_int .ne. 0) update_accel = .true.

  globalSumQuantities(:) = 0.0
  localSumQuantities (:) = 0.0

#ifdef CORRECT_BULK_MOTION

  ! sum quantities over list of blocks (to determine global mean force and momentum)
  do blockID = 1, blockCount

    ! get the index limits of the block
    call Grid_getBlkIndexLimits(blockList(blockID), blkLimits, blkLimitsGC)

    ! get a pointer to the current block of data
    call Grid_getBlkPtr(blockList(blockID), solnData)

    ! getting the dx's
    call Grid_getDeltas(blocklist(blockID), del)

#if NDIM == 1
    dvol = del(IAXIS)
#endif
#if NDIM == 2
    dvol = del(IAXIS) * del(JAXIS)
#endif
#if NDIM == 3
    dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
#endif

    ! update turbulent acceleration field, otherwise use previous acceleration field
    if (update_accel) then
      call Grid_getBlkPhysicalSize(blockList(blockID), blockSize)
      call Grid_getBlkCenterCoords(blockList(blockID), blockCenter)
      pos_beg = blockCenter - 0.5*blockSize + del/2.0 ! first active cell coordinate in block (x,y,z)
      pos_end = blockCenter + 0.5*blockSize - del/2.0 ! last  active cell coordinate in block (x,y,z)
      ncells = blkLimits(HIGH,:)-blkLimits(LOW,:)+1 ! number of active cells in (x,y,z)
      call st_get_turb_vector_unigrid_c(pos_beg, pos_end, ncells, accx, accy, accz)
    endif

    ! loop over all grid cells and sum local contributions to global mean force and momentum
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

          ii = i-blkLimits(LOW,IAXIS)+1 ! x index of accx, accy, accz starts at 1 and goes to NXB
          jj = j-blkLimits(LOW,JAXIS)+1 ! y index of accx, accy, accz starts at 1 and goes to NYB
          kk = k-blkLimits(LOW,KAXIS)+1 ! z index of accx, accy, accz starts at 1 and goes to NZB

#ifdef ACCX_VAR
          ! if we use ACCX_VAR, ..Y, ..Z (usually when using AMR because of re-gridding), copy from accx, ..y, ..z
          if (update_accel) then
            solnData(ACCX_VAR,i,j,k) = accx(ii,jj,kk)
            solnData(ACCY_VAR,i,j,k) = accy(ii,jj,kk)
            solnData(ACCZ_VAR,i,j,k) = accz(ii,jj,kk)
          endif
#endif

          ! cell mass
          dmass = solnData(DENS_VAR,i,j,k)*dvol

          ! mass
#ifdef DENS_VAR
          localSumQuantities(1) = localSumQuantities(1) + dmass

          ! momentum
#ifdef VELX_VAR
          localSumQuantities(2) = localSumQuantities(2) + solnData(VELX_VAR,i,j,k)*dmass
#endif
#ifdef VELY_VAR
          localSumQuantities(3) = localSumQuantities(3) + solnData(VELY_VAR,i,j,k)*dmass
#endif
#ifdef VELZ_VAR
          localSumQuantities(4) = localSumQuantities(4) + solnData(VELZ_VAR,i,j,k)*dmass
#endif
          ! driving force
#ifdef ACCX_VAR
          localSumQuantities(5) = localSumQuantities(5) + solnData(ACCX_VAR,i,j,k)*dmass
#else
          localSumQuantities(5) = localSumQuantities(5) + accx(ii,jj,kk)*dmass
#endif
#ifdef ACCY_VAR
          localSumQuantities(6) = localSumQuantities(6) + solnData(ACCY_VAR,i,j,k)*dmass
#else
          localSumQuantities(6) = localSumQuantities(6) + accy(ii,jj,kk)*dmass
#endif
#ifdef ACCZ_VAR
          localSumQuantities(7) = localSumQuantities(7) + solnData(ACCZ_VAR,i,j,k)*dmass
#else
          localSumQuantities(7) = localSumQuantities(7) + accz(ii,jj,kk)*dmass
#endif
#endif
! ifdef DENS_VAR
        enddo ! i
      enddo ! j
    enddo ! k

    call Grid_releaseBlkPtr(blockList(blockID),solnData)

  enddo ! blocks

  ! now communicate all global summed quantities to all processors
  call MPI_AllReduce(localSumQuantities, globalSumQuantities, nGlobalSum, &
                      MPI_DOUBLE_PRECISION, MPI_Sum, MPI_Comm_World, error)

  mass      = globalSumQuantities(1)
  xMomentum = globalSumQuantities(2)
  yMomentum = globalSumQuantities(3)
  zMomentum = globalSumQuantities(4)
  xForce    = globalSumQuantities(5)
  yForce    = globalSumQuantities(6)
  zForce    = globalSumQuantities(7)

#endif
! ifdef CORRECT_BULK_MOTION

  ! set to zero for adding block and cell contributions below
  ekin_added = 0.0

  ! Loop over local blocks again to actually apply the turbulent acceleration
  do blockID = 1, blockCount

    ! get the index limits of the block
    call Grid_getBlkIndexLimits(blockList(blockID), blkLimits, blkLimitsGC)

    ! get a pointer to the current block of data
    call Grid_getBlkPtr(blockList(blockID), solnData)

    ! getting the dx's
    call Grid_getDeltas(blocklist(blockID), del)

#if NDIM == 1
    dvol = del(IAXIS)
#endif
#if NDIM == 2
    dvol = del(IAXIS) * del(JAXIS)
#endif
#if NDIM == 3
    dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
#endif

! update acceleration if CORRECT_BULK_MOTION is not used (otherwise, st_get_turb_vector_unigrid_c was called above)
#ifndef CORRECT_BULK_MOTION
    ! update turbulent acceleration field, otherwise use previous acceleration field
    if (update_accel) then
      call Grid_getBlkPhysicalSize(blockList(blockID), blockSize)
      call Grid_getBlkCenterCoords(blockList(blockID), blockCenter)
      pos_beg = blockCenter - 0.5*blockSize + del/2.0 ! first active cell coordinate in block (x,y,z)
      pos_end = blockCenter + 0.5*blockSize - del/2.0 ! last  active cell coordinate in block (x,y,z)
      ncells = blkLimits(HIGH,:)-blkLimits(LOW,:)+1 ! number of active cells in (x,y,z)
      call st_get_turb_vector_unigrid_c(pos_beg, pos_end, ncells, accx, accy, accz)
    endif
#endif
! ifndef CORRECT_BULK_MOTION

    ! loop over all grid cells and apply turbulent acceleration
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

          ii = i-blkLimits(LOW,IAXIS)+1 ! x index of accx, accy, accz starts at 1 and goes to NXB
          jj = j-blkLimits(LOW,JAXIS)+1 ! y index of accx, accy, accz starts at 1 and goes to NYB
          kk = k-blkLimits(LOW,KAXIS)+1 ! z index of accx, accy, accz starts at 1 and goes to NZB

#ifndef CORRECT_BULK_MOTION
#ifdef ACCX_VAR
          ! if we use ACCX_VAR, ..Y, ..Z (usually when using AMR because of re-gridding), copy from accx, ..y, ..z
          if (update_accel) then
            solnData(ACCX_VAR,i,j,k) = accx(ii,jj,kk)
            solnData(ACCY_VAR,i,j,k) = accy(ii,jj,kk)
            solnData(ACCZ_VAR,i,j,k) = accz(ii,jj,kk)
          endif
#endif
#endif

#ifdef VELX_VAR
          ekin_old = 0.5*(solnData(VELX_VAR,i,j,k)**2+solnData(VELY_VAR,i,j,k)**2+solnData(VELZ_VAR,i,j,k)**2)
#endif

#ifdef VELX_VAR
#ifdef ACCX_VAR
          accel = solnData(ACCX_VAR,i,j,k)
#else
          accel = accx(ii,jj,kk)
#endif
#ifdef CORRECT_BULK_MOTION
          correction = xForce/mass*dt + xMomentum/mass
#endif
          solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + accel*dt - correction
#endif

#if NDIM > 1
#ifdef VELY_VAR
#ifdef ACCY_VAR
          accel = solnData(ACCY_VAR,i,j,k)
#else
          accel = accy(ii,jj,kk)
#endif
#ifdef CORRECT_BULK_MOTION
          correction = yForce/mass*dt + yMomentum/mass
#endif
          solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + accel*dt - correction
#endif
#endif

#if NDIM > 2
#ifdef VELZ_VAR
#ifdef ACCZ_VAR
          accel = solnData(ACCZ_VAR,i,j,k)
#else
          accel = accz(ii,jj,kk)
#endif
#ifdef CORRECT_BULK_MOTION
          correction = zForce/mass*dt + zMomentum/mass
#endif
          solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) + accel*dt - correction
#endif
#endif

#ifdef VELX_VAR
          ekin_new = 0.5*(solnData(VELX_VAR,i,j,k)**2+solnData(VELY_VAR,i,j,k)**2+solnData(VELZ_VAR,i,j,k)**2)
#endif
          ! compute energy injection
          d_ekin = ekin_new - ekin_old
#ifdef ENER_VAR
          ! update the total energy
          solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + d_ekin
#endif
#ifdef DENS_VAR
          ! add up total injected kinetic energy
          ekin_added = ekin_added + d_ekin * solnData(DENS_VAR,i,j,k) * dvol
#endif
#ifdef INJR_VAR
          ! fill kinetic energy injection rate d(0.5 rho v^2) / dt; note that rho=const here
          solnData(INJR_VAR,i,j,k) = d_ekin * solnData(DENS_VAR,i,j,k) / dt
#endif
        enddo ! i
      enddo ! j
    enddo ! k

    call Grid_releaseBlkPtr(blockList(blockID),solnData)

  enddo ! loop over blocks

  ! write time evolution of ekin_added to file
  if (io_integralFreq .gt. 0) then
    ! synchronise output frequency with IO_output's use of io_integralFreq
    check_for_io = .true.
    if (PRESENT(pass)) then
      if (pass .eq. 1) check_for_io = .false.
    endif
    if ((check_for_io) .and. (mod(dr_nstep+1, io_integralFreq) .eq. 0)) then
      if (abs(time - last_time) .gt. TINY(time)) then
        last_time = time
        ! sum up injected kinetic energy contributions from all blocks and processors
        ekin_added_red = 0.0
        call MPI_Reduce(ekin_added, ekin_added_red, 1, MPI_DOUBLE_PRECISION, MPI_Sum, MASTER_PE, MPI_Comm_World, error)
        ekin_added = ekin_added_red
        ! only MASTER_PE writes
        if (dr_globalMe == MASTER_PE) then
          open(funit, file=trim(outfile), position='APPEND')
          write(funit,'(10(1X,ES16.9))') time, dt, ekin_added, ekin_added/dt, &
                                        xForce, yForce, zForce, xMomentum, yMomentum, zMomentum
          close(funit)
        endif
      endif
    endif
  endif

  call Timers_stop ("Stir")

  return

end subroutine Stir
