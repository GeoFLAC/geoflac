! -*- F90 -*-
module params
implicit none

integer, parameter :: maxzone = 10   ! max nzonx and nzony
integer, parameter :: maxbc = 20   ! max # of bcs
integer, parameter :: maxzone_age = 32   ! max # of nzone_age
integer, parameter :: maxzone_layer = 10   ! max # of layers in a nzone_age column
integer, parameter :: maxph = 20   ! max # of phases
integer, parameter :: maxinh = 50   ! max # of inhomogeneities

real*8, parameter :: sec_year = 3.1558d+7  ! seconds in a year


integer :: nx,nz,nzonx,nzony,nelz_x(maxzone),nelz_y(maxzone), &
     ny_rem,mode_rem,ntest_rem,ivis_shape, &
     itype_melting,nelem_serp,nmass_update,nopbmax,nydrsides,nystressbc, &
     nofbc,nofside(maxbc),nbc1(maxbc),nbc2(maxbc),nbc(maxbc), &
     mix_strain,mix_stress,lastsave,lastout, &
     io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz, &
     io_sxz,io_pres,io_temp,io_phase,io_visc,io_unused,io_density, &
     io_src,io_diss,io_forc,io_hfl,io_topo, &
     irphase,irtemp,ircoord, &
     nphase,irheol(maxph), &
     nmarkers, &
     ix1(maxinh),ix2(maxinh),iy1(maxinh),iy2(maxinh),inphase(maxinh), &
     igeom(maxinh),inhom,i_prestress, &
     itherm,istress_therm,itemp_bc,ishearh, &
     nzone_age,ixtb1(maxzone_age),ixtb2(maxzone_age), &
     ictherm(maxzone_age), &
     nph_layer(maxzone_age),iph_col(maxzone_age, maxzone_layer), &
     iph_col_trans(maxzone_age), &
     if_hydro,nyhydro,iphsub, &
     movegrid,ndim,ifreq_visc,i_rey, &
     incoming_left,incoming_right, &
     ivis_present,idt_scale,ifreq_imasses,ifreq_rmasses, &
     nloop,ifreq_avgsr,nsrate, &
     ! magma migration related
     isolidus,igeotherm,icrust_melt,istart_profile,imagtop, &
     mzone_style

!$ACC declare create(nx,nz,nzonx,nzony,nelz_x(maxzone),nelz_y(maxzone), &
!$ACC     ny_rem,mode_rem,ntest_rem,ivis_shape, &
!$ACC     itype_melting,nelem_serp,nmass_update,nopbmax,nydrsides,nystressbc, &
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
!$ACC     nloop,ifreq_avgsr,nsrate, &
!$ACC     isolidus,igeotherm,icrust_melt,imagtop, &
!$ACC     mzone_style)

real*8 :: x0,z0,rxbo,rzbo,sizez_x(maxzone),sizez_y(maxzone), &
     dx_rem,angle_rem,topo_kappa,fac_kappa, &
     v_min,v_max,efoldc, &
     dxmin,dzmin, &
     prod_magma,rho_magma, &
     bca(maxbc),bcb(maxbc),bcc(maxbc),xReyn, &
     bcd(maxbc),bce(maxbc),bcf(maxbc),bcg(maxbc),bch(maxbc),bci(maxbc), &
     dt_scale,strain_inert,vbc,frac, &
     dt_maxwell,fracm, &
     dt_elastic,demf, &
     dtout_screen,dtout_file,dtsave_file, &
     visc(maxph),den(maxph),alfa(maxph),beta(maxph),pln(maxph), &
     acoef(maxph),eactiv(maxph),rl(maxph),rm(maxph), &
     plstrain1(maxph),plstrain2(maxph),fric1(maxph),fric2(maxph), &
     cohesion1(maxph),cohesion2(maxph), &
     dilat1(maxph),dilat2(maxph), &
     conduct(maxph),cp(maxph), &
     ts(maxph),tl(maxph),tk(maxph),fk(maxph), &
     ten_off,tau_heal,xinitaps(maxinh), &
     t_top,t_bot,hs,hr,bot_bc, &
     hc(maxzone_age,maxzone_layer), &
     age_1(maxzone_age),tp1(maxzone_age),tp2(maxzone_age), &
     g,pisos,drosub,damp_vis, &
     angle_mzone,fmagma_max,ratio_mantle_mzone, &
     latent_heat_magma,lambda_freeze,lambda_freeze_tdep, &
     weaken_ratio_plastic,weaken_ratio_viscous, &
     dtavg, &
     time,dt,time_max, &
     ! magma migration related
     temp_per,g_x0,g_amplitude,g_width,width_mzone,ratio_crust_mzone, &
     watercont,xmodalcpx,dumC(4),xDF(4),xDT(4)

!$ACC declare create(x0,z0,rxbo,rzbo,sizez_x(maxzone),sizez_y(maxzone), &
!$ACC     dx_rem,angle_rem,topo_kappa,fac_kappa, &
!$ACC     v_min,v_max,efoldc, &
!$ACC     dxmin,dzmin, &
!$ACC     prod_magma,rho_magma, &
!$ACC     bca(maxbc),bcb(maxbc),bcc(maxbc),xReyn, &
!$ACC     bcd(maxbc),bce(maxbc),bcf(maxbc),bcg(maxbc),bch(maxbc),bci(maxbc), &
!$ACC     dt_scale,strain_inert,vbc,frac, &
!$ACC     dt_maxwell,fracm, &
!$ACC     dt_elastic,demf, &
!$ACC     dtout_screen,dtout_file,dtsave_file, &
!$ACC     visc(maxph),den(maxph),alfa(maxph),beta(maxph),pln(maxph), &
!$ACC     acoef(maxph),eactiv(maxph),rl(maxph),rm(maxph), &
!$ACC     plstrain1(maxph),plstrain2(maxph),fric1(maxph),fric2(maxph), &
!$ACC     cohesion1(maxph),cohesion2(maxph), &
!$ACC     dilat1(maxph),dilat2(maxph), &
!$ACC     conduct(maxph),cp(maxph), &
!$ACC     ts(maxph),tl(maxph),tk(maxph),fk(maxph), &
!$ACC     ten_off,tau_heal,xinitaps(maxinh), &
!$ACC     t_top,t_bot,hs,hr,bot_bc, &
!$ACC     hc(maxzone_age,maxzone_layer), &
!$ACC     age_1(maxzone_age),tp1(maxzone_age),tp2(maxzone_age), &
!$ACC     g,pisos,drosub,damp_vis, &
!$ACC     angle_mzone,fmagma_max,ratio_mantle_mzone, &
!$ACC     latent_heat_magma,lambda_freeze,lambda_freeze_tdep, &
!$ACC     weaken_ratio_plastic,weaken_ratio_viscous, &
!$ACC     dtavg, &
!$ACC     time,dt,time_max, &
!$ACC     temp_per,g_x0,g_amplitude,g_width,width_mzone,ratio_crust_mzone, &
!$ACC     watercont,xmodalcpx,dumC(4),xDF(4),xDT(4))

character phasefile*20,tempfile*20,coordfile*20

end module params
