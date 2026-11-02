!==============================================
! Density
function Eff_dens( j, i)
  !$ACC routine seq
  use arrays
  use params
  use phases
  implicit none
  integer, intent(in) :: j, i
  double precision :: Eff_dens
  integer :: iph, ii, k
  double precision :: adi_grad, den_ramp_dist, amp410, amp660, zcord, tmpr0, zsurf, tmpr, deptmpr410, deptmpr420, den_amp, deptmpr660, deptmpr670, press, ratio, dens, fma
  ! Adiabatic gradient
  adi_grad = 3.d-1
  ! Density transition thickness
  den_ramp_dist = 30.d3
  ! Density enhancement (in percentage) from the upper mantle to the transition zone
  ! Based on PREM model
  amp410 = 0.05d0
  ! Density enhancement (in percentage) from the transition zone to the lower mantle
  ! Based on PREM model
  amp660 = 0.085d0

  zcord = 0.25d0*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
  tmpr0 = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
  zsurf = 0.50d0*(cord(1,i,2)+cord(1,i+1,2))
  tmpr = (zsurf - zcord)*1.d-3*adi_grad + tmpr0
  iph = iphase(j,i)

  ! Phase diagram taken from The Solid Earth by C.M.R.Fowler, 2nd edition
  ! Fixed points (996.031431 K, 318.012214 km) (2338.655701 K, 454.304204 km)
  deptmpr410 = -101.58d0*(tmpr+273.d0)-216828.d0
  deptmpr420 = deptmpr410-den_ramp_dist
  
  den_amp = 1.0d0
  if (deptmpr410 > zcord) then
      frac = (deptmpr410 - zcord) / (deptmpr410 - deptmpr420)
      if (frac < 0.d0) frac=0.d0
      if (frac > 1.d0) frac=1.d0
      den_amp=1.0d0 + amp410 * frac
  endif
  
  ! Phase diagram taken from The Solid Earth by C.M.R.Fowler, 2nd edition
  ! Fixed points (1021.586179 K, 720.499495 km) (2333.396576 K, 633.187440 km)
  deptmpr660 = 66.55d0*(tmpr+273.d0)-788494.d0
  deptmpr670 = deptmpr660-den_ramp_dist
  if (deptmpr660 > zcord) then
      frac = (deptmpr660 - zcord) / (deptmpr660 - deptmpr670)
      if (frac < 0.d0) frac=0.d0
      if (frac > 1.d0) frac=1.d0
      den_amp=1.05d0 + amp660 * frac
  endif

  press = 0
  do ii = 1, 4
      press = press - (stress0(j,i,1,ii)+stress0(j,i,2,ii)+stress0(j,i,4,ii))
  enddo
  press = press / 12

  Eff_dens = 0.d0
  do k = 1, nphase
    ratio = phase_ratio(k,j,i)
    ! when ratio is small, it won't affect the density
    if(ratio .lt. 0.01d0) cycle

    dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

    !press = 0
    !press = stressI(j,i)
    press = dens*g*zcord

    dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

    Eff_dens = Eff_dens + ratio*dens
  enddo

  if (any(iph == mantle_phases)) then
      Eff_dens = Eff_dens * den_amp
  endif
  fma = min(fmagma(j,i)+fmagma2(j,i),fmagma_max)
  Eff_dens = Eff_dens - fma * (Eff_dens - rho_magma)
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

integer, intent(in) :: j, i
integer :: iph
double precision :: Eff_cp

iph = iphase(j,i)
Eff_cp = cp(iph)

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

integer, intent(in) :: j, i
integer :: iph
double precision :: Eff_conduct

iph = iphase(j,i)
Eff_conduct = conduct(iph)

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
implicit none
integer, intent(in) :: j, i
double precision :: Eff_visc
integer :: k
double precision :: adi_grad, r, zcord, tmpr0, zsurf, tmpr, s11, s22, s33, pres, deptmpr410, deptmpr660, srat, pow, pow1, vis, fma
! Adiabatic gradient
adi_grad = 3.d-1

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

Eff_visc = 0.d0
r=8.31448d0

zcord = 0.25d0*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
tmpr0 = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
zsurf = 0.50d0*(cord(1,i,2)+cord(1,i+1,2))
tmpr = (zsurf - zcord)*1.d-3*adi_grad + tmpr0

s11 = 0.25d0 * (stress0(j,i,1,1)+stress0(j,i,1,2)+stress0(j,i,1,3)+stress0(j,i,1,4))
s22 = 0.25d0 * (stress0(j,i,2,1)+stress0(j,i,2,2)+stress0(j,i,2,3)+stress0(j,i,2,4))
s33 = 0.25d0 * (stress0(j,i,4,1)+stress0(j,i,4,2)+stress0(j,i,4,3)+stress0(j,i,4,4))
pres = -1*(s11+s22+s33)/3d0

! Phase diagram taken from The Solid Earth by C.M.R.Fowler, 2nd edition
! Fixed points (996.031431 K, 318.012214 km) (2338.655701 K, 454.304204 km)
deptmpr410 = -101.58d0*(tmpr+273.d0)-216828.d0
! Fixed points (1021.586179 K, 720.499495 km) (2333.396576 K, 633.187440 km)
deptmpr660 = 66.55d0*(tmpr+273.d0)-788494.d0

if (zcord < deptmpr660) then
    Eff_visc = 1d22 * exp( -0.0115d0 * (tmpr - t_bot))
    Eff_visc = min(v_max, max(v_min, Eff_visc))
    return
elseif (zcord < deptmpr410) then
    Eff_visc = 1d21 * exp( -0.0115d0 * (tmpr - t_bot))
    Eff_visc = min(v_max, max(v_min, Eff_visc))
    return
endif

srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

do k = 1, nphase
    if(phase_ratio(k,j,i) .lt. 0.01d0) cycle
    pow  =  1.d0/pln(k) - 1.d0
    pow1 = -1.d0/pln(k)

    vis = 0.25d0 * srat**pow*(0.75d0*acoef(k))**pow1* &
          exp((eactiv(k)+vactiv(k)*pres)/(pln(k)*r*(tmpr+273.d0)))*1.d+6

    if (vis .lt. v_min) vis = v_min
    if (vis .gt. v_max) vis = v_max

    ! harmonic mean
    Eff_visc = Eff_visc + phase_ratio(k,j,i) / vis
    !write(*,*) i,j, Eff_visc, vis, tmpr,phase_ratio(k,j,i)
enddo

Eff_visc = 1 / Eff_visc
fma = min(fmagma(j,i)+fmagma2(j,i),fmagma_max)
if (itype_melting == 1) Eff_visc = Eff_visc * exp(weaken_ratio_viscous * fma / fmagma_max)

! Final cut-off
Eff_visc = min(v_max, max(v_min, Eff_visc))

return
end function Eff_visc
