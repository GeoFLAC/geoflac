!===============================================
!   Lithostatic pressure
function xLith_dp( j, i)
use arrays
use params
use phases
include 'precision.inc'

iph = iphase(j,i)
densT = Eff_dens(j,i)
dh1 = cord (j,i,2) - cord (j+1,i,2)
dh2 = cord (j,i+1,2) - cord (j+1,i+1,2)
dh  = 0.5d0 * (dh1+dh2)
dPT = densT * g * dh
xLith_dp = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )

return
end function xLith_dp

!==============================================
! Density
function Eff_dens( j, i)
!$ACC routine seq
use arrays
use params
use phases
include 'precision.inc'

zcord = 0.25d0*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
tmpr = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

! 660-km phase change with clapeyron slope 600C/50km
if (zcord < -660d3 + (tmpr - t_bot) * 500.d0/6) then
    Eff_dens = 3450 * (1 - 2d-5 * tmpr)  ! ~6% density jump
    return
endif

press = 0.

if (itype_melting.ne.2) then
    do ii = 1, 4
        press = press - (stress0(j,i,1,ii)+stress0(j,i,2,ii)+stress0(j,i,4,ii))
    enddo
    press = press / 12
else
    ! Lithostatic pressure
    rogh  = 0.
    tmpr1= 0.
    do jj = 1,j
        iph = iphase(jj,i)
        tmpr1 = 0.25*(temp(jj,i)+temp(jj+1,i)+temp(jj,i+1)+temp(jj+1,i+1))
        densT= den(iph) * (1.-alfa(iph)*tmpr1)
        dh1 = cord(jj,i,2)-cord(jj+1,i,2) 
        dh2 = cord(jj,i+1,2)-cord(jj+1,i+1,2) 
        dh  = 0.5*(dh1+dh2)
        dPT = densT*g*dh
        dP = dPT*(1.-beta(iph)*rogh)/(1.+beta(iph)/2.*dPT)
        press = rogh + 0.5*dP
        rogh = rogh + dP
    enddo
endif

Eff_dens = 0.d0
do k = 1, nphase
    ratio = phase_ratio(k,j,i)
    ! when ratio is small, it won't affect the density
    if(ratio .lt. 0.01d0) cycle

    dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

    if (itype_melting.ne.2) then
        press = dens*g*zcord
        dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )
    else
        if (k == ksed1 .or. k == ksed2) then
            if(k == ksed1) then
                sed_min_density = 2200.
                zefold = 6000.
            else if(k == ksed2) then
                sed_min_density = 2600.
                zefold = 1500.
            endif

            press = stressI(j,i)

            delta_den = den(k) - sed_min_density

            if (j==1 .and. cord(j,i,2)>0.) then
                ! sediment above ground is already compacted
                dens = den(k) * (1.-alfa(k)*tmpr)
            else
                ! sediment below sea level is uncompacted and less dense
                dens = (den(k) - delta_den*exp((zcord-0.5*(cord(1,i,2)+cord(1,i+1,2)))/zefold)) &
                        * ( 1 - alfa(k)*tmpr + beta(k)*press )
            endif
            if (dens < sed_min_density) dens = sed_min_density
        endif

        ! if(k == ksed1) then
        !     sed_min_density = 2200.
        !     delta_den = den(k) - sed_min_density
        !     zefold = 6000.
        !     if (j==1 .and. cord(j,i,2)>0.) then
        !         ! sediment above ground is already compacted
        !         dens = den(k) * (1.-alfa(k)*tmpr)
        !     else
        !         ! sediment below sea level is uncompacted and less dense
        !         dens = (den(k) - delta_den*exp((zcord-0.5*(cord(1,i,2)+cord(1,i+1,2)))/zefold)) &
        !                 * ( 1 - alfa(k)*tmpr + beta(k)*press )
        !     endif
        !     if (dens < sed_min_density) dens = sed_min_density
        ! endif
        
        ! if(k == ksed2) then
        !     sed_min_density = 2600.
        !     delta_den = den(k) - sed_min_density
        !     zefold = 1500.
        !     if (j==1 .and. cord(j,i,2)>0.) then
        !         ! sediment above ground is already compacted
        !         dens = den(k) * (1.-alfa(k)*tmpr)
        !     else
        !         ! sediment below sea level is uncompacted and less dense
        !         dens = (den(k) - delta_den*exp((zcord-0.5*(cord(1,i,2)+cord(1,i+1,2)))/zefold)) &
        !                 * ( 1 - alfa(k)*tmpr + beta(k)*press )
        !     endif
        !     if (dens < sed_min_density) dens = sed_min_density
        ! endif

        if ((xfmelt(j,i)  + Eff_melt(j,i)).gt.0.95) xfmelt(j,i) = 0.

        dens = dens * ( 1.- (xfmelt(j,i)  + Eff_melt(j,i))) + rho_magma*(xfmelt(j,i)  + Eff_melt(j,i))
    endif

    Eff_dens = Eff_dens + ratio*dens

enddo

return
end function Eff_dens



!=================================================
! Effective Heat Capacity incorporating latent heat
!=================================================
function Eff_cp( j, i )
!$ACC routine seq
use arrays
use params
implicit none

integer :: iph, j, i
double precision :: Eff_cp, HeatLatent, tmpr

HeatLatent = 420000.

iph = iphase(j,i)

if (itype_melting .ne. 2) then
    Eff_cp = cp(iph)
else
    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    if( tmpr .lt. ts(iph) ) then
        Eff_cp = cp(iph)
    elseif( tmpr .lt. tk(iph) ) then
        Eff_cp = cp(iph) + HeatLatent * fk(iph)/(tk(iph)-ts(iph))
    elseif( tmpr .lt. tl(iph) ) then
        Eff_cp = cp(iph) + HeatLatent * (1.-fk(iph))/(tl(iph)-tk(iph))
    else
        Eff_cp = cp(iph)
    endif
endif

! HOOK
! Intrusions - melting effect - see user_ab.f90
! if( if_intrus .eq. 1 ) then
!     HeatLatent = 420000.

!     tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
!     if( tmpr .lt. ts(iph) ) then
!         Eff_cp = cp(iph)
!     elseif( tmpr .lt. tl(iph)+1 ) then
!         Eff_cp = cp(iph) + HeatLatent/(tl(iph)-ts(iph))
!     else
!         Eff_cp = cp(iph)
!     endif
! endif

return
end function Eff_cp


!=================================================
! Effective Thermal Conductivity
!=================================================
function Eff_conduct( j, i )
!$ACC routine seq
use arrays
use params
implicit none

double precision, external :: HydroDiff
double precision :: Eff_conduct, cond
integer :: j, i, k, iph

if (itype_melting .ne. 2) then
    iph = iphase(j,i)
    Eff_conduct = conduct(iph)
else
    Eff_conduct = 0.
    do k = 1 , nphase

        ! when ratio is small, it won't affect the density
        if(phase_ratio(k,j,i) .lt. 0.01) cycle

        cond = conduct(k)

        ! HOOK
        ! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
        if( if_hydro .eq. 1 ) then
            cond = HydroDiff(j,i)*den(k)*cp(k)
        endif
        Eff_conduct = Eff_conduct + phase_ratio(k,j,i)*cond
    enddo
endif

return
end function Eff_conduct



!=================================================
! Non-Newtonian viscosity
!=================================================

! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function Eff_visc( j, i )
!$ACC routine seq
use arrays
use params
use phases
use marker_data
include 'precision.inc'

if (itype_melting .ne. 2) then
    do k = 1, inhom
        !Fixed weak zone
        if (igeom(k).eq.100) then
            if (i >= ix1(k) .and. i <= ix2(k) .and. &
                j >= iy1(k) .and. j <= iy2(k)) then
                Eff_visc = v_min
                return
            endif
        endif
    enddo
endif

Eff_visc = 0.d0
r = 8.31448d0
tmpr = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

if (itype_melting .ne. 2) then
    zcord = 0.25d0*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
    if (zcord < -660d3) then
        ! T-dep. viscosity
        ! eta=1e22 when T=1330 Celsius, inc. 10x when T is 200 deg lower
        Eff_visc = 1d22 * exp( -0.0115d0 * (tmpr - 1330))
        Eff_visc = min(v_max, max(v_min, Eff_visc))
        return
    elseif (zcord < -410d3) then
        Eff_visc = 1d21 * exp( -0.0115d0 * (tmpr - 1330))
        Eff_visc = min(v_max, max(v_min, Eff_visc))
        return
    endif
else
    ! Lithostatic pressure
    preslith = 0.
    rogh  = 0.
    do jj = 1,j
        iph = iphase(jj,i)
        tmpr1 = 0.25*(temp(jj,i)+temp(jj+1,i)+temp(jj,i+1)+temp(jj+1,i+1))
        densT= den(iph) * (1.-alfa(iph)*tmpr1)
        dh1 = cord(jj,i,2)-cord(jj+1,i,2) 
        dh2 = cord(jj,i+1,2)-cord(jj+1,i+1,2) 
        dh  = 0.5*(dh1+dh2)
        dPT = densT*g*dh
        dP = dPT*(1.-beta(iph)*rogh)/(1.+beta(iph)/2.*dPT)
        preslith = rogh + 0.5*dP
        rogh = rogh + dP
    enddo
endif

do k = 1, nphase
    if(phase_ratio(k,j,i) .lt. 0.01) cycle

    pow  =  1./pln(k) - 1.
    pow1 = -1./pln(k)

    if (itype_melting .ne. 2) then
        vis = 0.25d0 * srat**pow*(0.75d0*acoef(k))**pow1* &
            exp(eactiv(k)/(pln(k)*r*(tmpr+273.d0)))*1.d+6
    else
        vis = 0.25*srat**pow*(0.75*acoef(k))**pow1* &
            exp((eactiv(k)+27.e-6*preslith)/(pln(k)*r*(tmpr+273.)))*1.e+6
        if (k.eq.khtmsz) then
            xxgz = 15.*0.25*(stressII(j,i)/1.e6)**(-4./3.)
            vis=(0.75*7.7e-2)**(-1.)*(xxgz**(2.))*exp(290.e3/r/(tmpr+273.))*1.e6
        elseif (k.eq.kocean2 .or. k.eq.kocean0) then
            xlava_age = 0.
            ! Lava flow at the surface must harden over the heat diffusion time scale L=1000 m, diff = 1e-6 m^2s-1
            ! Td cooling lava over 1 km (element size) is 1e6/1e-6 = 1e12 s.
            Tcool = 1.e12
            ! we assume the lava flow arrive at a temperature of ~1000 C (it is more like 1250 C in lava tubes)
            ! The initial viscosity is vis_1000 which decays exponentially since its formation xlava_age:
            kinc = nmark_elem(1,i) ! number of marker in element 
            nl = 0  ! counts the lava markers
            do n = 1,kinc
                kid = mark_id_elem(n,j,i)
                if (mark_phase(kid).eq.7) then
                nl = nl + 1
                xlava_age = xlava_age + mark_age(kid)
                endif
            enddo
            ! age since extrusion
            xlava_age = xlava_age/nl
            !  if (tmpr.gt.1000.) tmpr = 1000. ! Make sure we are cooling off
            T_lava = max(1200.*exp(-xlava_age/Tcool),tmpr) 
            vis = 0.25*srat**pow*(0.75*acoef(k))**pow1* &
            exp((eactiv(k)+27.e-6*preslith)/(pln(k)*r*(T_lava+273.)))*1.e+6
        endif
            

        if (xfmelt(j,i)+Eff_melt(j,i).gt.0.5) xfmelt(j,i) = 0.

        vis = vis*dexp(-30.*(xfmelt(j,i)+Eff_melt(j,i))) 
    endif

    ! Final cut-off
    if (vis .lt. v_min) vis = v_min
    if (vis .gt. v_max) vis = v_max
    ! harmonic mean
    Eff_visc = Eff_visc + phase_ratio(k,j,i) / vis
enddo

Eff_visc = 1 / Eff_visc
if (itype_melting == 1) Eff_visc = Eff_visc * exp(weaken_ratio_viscous * fmagma(j,i) / fmagma_max)

! Final cut-off
Eff_visc = min(v_max, max(v_min, Eff_visc))

return
end function Eff_visc
