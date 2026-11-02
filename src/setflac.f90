
! Setup some parameters (rmass,amass,initial stress,vel,viscosity)

subroutine setflac
use arrays
use params
implicit none
nloop = 0
time = 0.d0
istart_profile = 50
!$ACC update device(nloop,time)

! Mesh generator
call init_cord

! Initial accumulated plastic strain
aps = 0

! Initial velocity
vel = 0

dvol = 0
strain = 0

! Phases in the mesh
call init_phase

! Setup markers
call init_marker

! Check if viscous rheology present
call check_visc_rheol

! Inverse Areas of triangles
call init_areas

! Minimum element area of the initial (t=0) grid; used later as the
! per-timestep cap on dike-intrusion eigenstrain release. Computed once
! here (never recomputed after remeshing) and persisted across restarts
! (see saveflac.f90/rsflac.f90), since a restart's current grid is no
! longer the initial one.
Amin = minval(0.5d0/area(:,:,1) + 0.5d0/area(:,:,2))

! Initiate temperature field
call init_temp

! Calculation of the initial STRESSES (as hydrostatic)
call init_stress

! Setup boundary conditions
call init_bc

temp0 = temp
xoriginal = cord(:,:,1)
zoriginal = cord(:,:,2)
shrheat = 0
sshrheat = 0
dtopo = 0
extrusion = 0
fmelt = 0
fmelt2 = 0
fmagma = 0
fmagma2 = 0
e2sr = 1d-16
se2sr = 1d-16
q_init = 0.d0

call update_acc

! Distribution of REAL masses to nodes
call rmasses

! Initialization of viscosity
if( ivis_present.eq.1 ) call init_visc

! Inertial masses and time steps (elastic and maxwell)
call dt_mass

! Initiate parameters for stress averaging
dtavg=0
nsrate=-1
! dtavg/nsrate are reset after update_acc's bulk push above (which must
! run before rmasses/init_visc/dt_mass, since those need nx, nz, and the
! other params scalars/small arrays already on device); push these two
! explicitly so the device copy doesn't retain a stale pre-reset value.
!$ACC update device(dtavg, nsrate) async(1)

return
end


subroutine update_acc
use arrays
use params
implicit none
!$ACC update device(nx,nz,nzonx,nzony,nelz_x(maxzone),nelz_y(maxzone), &
!$ACC     ny_rem,mode_rem,ntest_rem,ivis_shape, &
!$ACC     itype_melting,nelem_serp,nelem_dike,nmass_update,nopbmax,nydrsides,nystressbc, &
!$ACC     nofbc,nofside(maxbc),nbc1(maxbc),nbc2(maxbc),nbc(maxbc), &
!$ACC     mix_strain,mix_stress,lastsave,lastout, &
!$ACC     io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz, &
!$ACC     io_sxz,io_pres,io_temp,io_phase,io_visc,io_unused,io_density, &
!$ACC     io_src,io_diss,io_forc,io_hfl,io_topo, &
!$ACC     irphase,irtemp,ircoord, &
!$ACC     nphase,irheol(maxph), &
!$ACC     nmarkers, &
!$ACC     ix1(maxinh),ix2(maxinh),iy1(maxinh),iy2(maxinh),inphase(maxinh), &
!$ACC     igeom(maxinh),inhom,i_prestress, &
!$ACC     itherm,istress_therm,itemp_bc,ishearh, &
!$ACC     nzone_age,ixtb1(maxzone_age),ixtb2(maxzone_age), &
!$ACC     ictherm(maxzone_age), &
!$ACC     nph_layer(maxzone_age),iph_col(maxzone_age, maxzone_layer), &
!$ACC     iph_col_trans(maxzone_age), &
!$ACC     if_hydro,nyhydro,iphsub, &
!$ACC     movegrid,ndim,ifreq_visc,i_rey, &
!$ACC     incoming_left,incoming_right, &
!$ACC     ivis_present,idt_scale,ifreq_imasses,ifreq_rmasses, &
!$ACC     nloop,ifreq_avgsr,nsrate) async(1)

!$ACC update device(x0,z0,rxbo,rzbo,sizez_x(maxzone),sizez_y(maxzone), &
!$ACC     dx_rem,angle_rem,topo_kappa,fac_kappa, &
!$ACC     v_min,v_max,efoldc, &
!$ACC     dxmin,dzmin, &
!$ACC     prod_magma,rho_magma,prod_magma2, &
!$ACC     bca(maxbc),bcb(maxbc),bcc(maxbc),xReyn, &
!$ACC     bcd(maxbc),bce(maxbc),bcf(maxbc),bcg(maxbc),bch(maxbc),bci(maxbc), &
!$ACC     dt_scale,strain_inert,vbc,frac, &
!$ACC     dt_maxwell,fracm, &
!$ACC     dt_elastic,demf, &
!$ACC     dtout_screen,dtout_file,dtsave_file, &
!$ACC     visc(maxph),den(maxph),alfa(maxph),beta(maxph),pln(maxph), &
!$ACC     acoef(maxph),eactiv(maxph),vactiv(maxph),rl(maxph),rm(maxph), &
!$ACC     plstrain1(maxph),plstrain2(maxph),fric1(maxph),fric2(maxph), &
!$ACC     cohesion1(maxph),cohesion2(maxph), &
!$ACC     dilat1(maxph),dilat2(maxph), &
!$ACC     conduct(maxph),cp(maxph), &
!$ACC     ts(maxph),tl(maxph),tk(maxph),erodibility(maxph), &
!$ACC     ten_off,tau_heal,xinitaps(maxinh), &
!$ACC     extra_pres,t_top,t_bot,hs,hr,bot_bc, &
!$ACC     hc(maxzone_age,maxzone_layer), &
!$ACC     age_1(maxzone_age),tp1(maxzone_age),tp2(maxzone_age), &
!$ACC     g,pisos,drosub,damp_vis, &
!$ACC     angle_mzone,fmagma_max,ratio_mantle_mzone,ratio_crust_mzone,Amin, &
!$ACC     latent_heat_magma,lambda_freeze,lambda_freeze_tdep, &
!$ACC     weaken_ratio_plastic,weaken_ratio_viscous, &
!$ACC     dtavg, &
!$ACC     time,dt,time_max) async(1)
end
