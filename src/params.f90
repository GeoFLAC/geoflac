! -*- F90 -*-
module params
implicit none

integer, parameter :: maxzone = 10   ! max nzonx and nzony
integer, parameter :: maxbc = 20   ! max # of bcs
integer, parameter :: maxph = 20   ! max # of phases
integer, parameter :: maxinh = 50   ! max # of inhomogeneities

integer :: nx,nz,nq,nzonx,nzony,nelz_x(maxzone),nelz_y(maxzone), &
     ny_rem,mode_rem,ntest_rem,ivis_shape,igeotherm, &
     ny_inject,nelem_inject,nmass_update,nopbmax,iinj, &
     nofbc,nofside(maxbc),nbc1(maxbc),nbc2(maxbc),nbc(maxbc), &
     mix_strain,mix_stress,lastsave,lastout, &
     io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz, &
     io_sxz,io_pres,io_temp,io_melt,io_visc,io_phas,io_mark, &
     io_src,io_diss,io_forc,io_hfl,io_topo, &
     irh_mark,irphase,irtemp,ircoord, &
     nphase,mphase,irheol(maxph), &
     ltop(20),lbottom(20),lphase(20), &
     imx1(20),imx2(20),imy1(20),imy2(20), &
     itx1(20),itx2(20),ity1(20),ity2(20), &
     nphasl,nzone_marker,nmarkers, iint_marker,iint_tracer,nzone_tracer, &
     ix1(maxinh),ix2(maxinh),iy1(maxinh),iy2(maxinh),inphase(maxinh), &
     igeom(maxinh),inhom, &
     itherm,istress_therm,initial_geoth,itemp_bc,ix1t,ix2t,iy1t,iy2t,ishearh, &
     ntherm,ixtb1(maxph),ixtb2(maxph),nzone_age,i_prestress, &
     iph_col1(maxph),iph_col2(maxph),iph_col3(maxph),iph_col4(maxph), &
     iph_col5(maxph),iph_col_trans(maxph), &
     if_intrus,if_hydro, &
     nyhydro,iphsub, &
     nloop_restarted,nloop,irestart,ifreq_avgsr,ireset, &
     devmax, dvmax

real*8 :: x0,z0,rxbo,rzbo,sizez_x(maxzone),sizez_y(maxzone), &
     dx_rem,angle_rem,anglemin1,anglemint,topo_kappa,fac_kappa, &
     velbc_l,velbc_r,v_min,v_max,efoldc, &
     g_x0,g_y0c,g_amplitude,g_width,ydrsides,ynstressbc, &
     rate_inject, &
     bca(maxbc),bcb(maxbc),bcc(maxbc),dx_init,xReyn, &
     bcd(maxbc),bce(maxbc),bcf(maxbc),bcg(maxbc),bch(maxbc),bci(maxbc), &
     dt_scale,strain_inert,vbc,amul,ratl,ratu,frac, &
     dtmax_therm,dt_maxwell,fracm,srate0, &
     dt_elastic,dt_elastic0,demf,boff, &
     ivis_present,n_boff_cutoff,idt_scale,ifreq_imasses,ifreq_rmasses, &
     movegrid,ndim,ifreq_visc,nmtracers,i_rey,iac_rem, &
     iph_int,incoming_left,incoming_right, &
     iynts,tbos, iax1,iay1,ibx1,iby1,icx1,icy1,idx1,idy1, &
     dtout_screen,dtout_file,dtsave_file, &
     visc(maxph),den(maxph),alfa(maxph),beta(maxph),pln(maxph), &
     acoef(maxph),eactiv(maxph),rl(maxph),rm(maxph), &
     plstrain1(maxph),plstrain2(maxph),fric1(maxph),fric2(maxph), &
     cohesion1(maxph),cohesion2(maxph), &
     dilat1(maxph),dilat2(maxph), &
     conduct(maxph),cp(maxph), &
     ts(maxph),tl(maxph),tk(maxph),fk(maxph), &
     ten_off,tau_heal,dt_outtracer,xinitaps(maxinh), &
     t_top,t_bot,hs,hr,temp_per,bot_bc,time_t,dt_therm, &
     hc1(maxph),hc2(maxph),hc3(maxph),hc4(maxph), &
     age_1(maxph),g,pisos,drosub,damp_vis, &
     sec_year,time,dt,time_max

character phasefile*20,tempfile*20,coordfile*20

end module params