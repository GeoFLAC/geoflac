!-*- F90 -*-

! --------- Flac -------------------------

subroutine flac

use arrays
use params
#ifdef USE_NVTX
use nvtx_mod
#endif
implicit none
! Update Thermal State
! Skip the therm calculations if itherm = 3
#ifdef USE_NVTX
if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('fl_therm')
#endif
call fl_therm
#ifdef USE_NVTX
if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif

if (itherm /= 2) then
    ! Changing marker phases
    ! XXX: change_phase is slow, don't call it every loop
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('change_phase')
#endif
    if( mod(nloop, ifreq_rmasses).eq.0 ) call change_phase
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif

    ! Update stresses by constitutive law (and mix isotropic stresses)
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('fl_rheol')
#endif
    call fl_rheol
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif

    ! update stress boundary conditions
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('bc_update')
#endif
    if (nystressbc.eq.1) call bc_update
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif

    ! Calculations in a node: forces, balance, velocities, new coordinates
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('fl_node')
#endif
    call fl_node
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif

    ! New coordinates
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('fl_move')
#endif
    call fl_move
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif

    ! Adjust real masses due to temperature
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('rmasses')
#endif
    if( mod(nloop,ifreq_rmasses).eq.0 ) call rmasses
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif

    ! Adjust inertial masses or time step due to deformations
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxStartRange('dt_mass')
#endif
    if( mod(nloop,ifreq_imasses) .eq. 0 ) call dt_mass
#ifdef USE_NVTX
    if (mod(nloop, istart_profile) .lt. 10) call nvtxEndRange()
#endif
endif


return
end subroutine flac
