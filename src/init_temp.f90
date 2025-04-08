! Initiate temperature profile

subroutine init_temp
use arrays
use params
use phases
use marker_data
implicit none

integer :: i, j, n, i1, i2, k, kk, ixc, iwidth
double precision :: age_1n, tp1n, tp2n, ratio, amp, y, pert, pert2

double precision :: cond_c, cond_m, dens_c, dens_m, pi, diffusivity
double precision :: F, am, alc, xc, yc, xsfh, auc, xdz, zPotT
integer :: kl, dum
double precision :: shf, huc, hlc
double precision :: Q(nz,nx)

!  Read distribution of temperatures from the dat file
if (irtemp .gt. 0) then
    open( 1, file=tempfile, status='old', err=101 )
    do i = 1, nx
    do j = 1, nz
        read( 1, * ,err=102 ) temp(j,i)
    enddo
    enddo
    close(1)

    goto 10
    101 call SysMsg('INIT_TEMP: Cannot open file with temperatures initial distrib!')
    stop 21
    102 call SysMsg('INIT_TEMP: Error reading file with temperatures initial distrib!')
    stop 21
endif

!!  geotherm of a given age accross the box with variable age
do n = 1, nzone_age
    if (n /= 1) then
        if (iph_col_trans(n-1) == 1) cycle
    endif

    do i = ixtb1(n), ixtb2(n)
        if (iph_col_trans(n) == 1) then
            i1 = ixtb1(n)
            i2 = ixtb2(n)
            ratio = (cord(1,i,1) - cord(1,i1,1)) / (cord(1,i2,1) - cord(1,i1,1))
            age_1n = age_1(n) + (age_1(n+1) - age_1(n)) * ratio
            tp1n = tp1(n) + (tp1(n+1) - tp1(n)) * ratio
            tp2n = tp2(n) + (tp2(n+1) - tp2(n)) * ratio
        else
            age_1n = age_1(n)
            tp1n = tp1(n)
            tp2n = tp2(n)
        endif
        call init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
    enddo
enddo

10 continue

! DISTRIBUTE SOURCES in elements
do j = 1,nz-1
    y = -( cord(j+1,1,2)+cord(j,1,2) )/2 / 1000
    source(j,1:nx-1) = hs*exp(-y/hr)
end do

do i = 1, inhom
    ! Initial gaussian temperature perturbation

    ! vertical gaussian
    ! x between (ix1, ix2), y between (iy1, iy2)
    ! gaussian dist. in x direction
    ! linear dist in y direction
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.11) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do j = ix1(i),ix2(i)
            pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
            do k = iy1(i),iy2(i)
                pert2 = 1.0d0*(k-iy1(i)) / (iy2(i) - iy1(i))
                temp(k,j) = min(t_bot, temp(k,j)+pert*pert2)
            enddo
        enddo
    endif

    ! slant gaussian
    ! x between (ix1, ix2) at top, shift 1-grid to right for every depth grid
    ! z between (iy1, iy2)
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.13) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do k = iy1(i),iy2(i)
            kk = k - iy1(i)
            do j = ix1(i),ix2(i)
                pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
                temp(k,j+kk) = max(t_top, min(t_bot, temp(k,j+kk)+pert))
                !print *, k, j, pert
            enddo
        enddo
    endif

    ! slant gaussian
    ! x between (ix1, ix2) at top, shift 1-grid to left for every depth grid
    ! z between (iy1, iy2)
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.14) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do k = iy1(i),iy2(i)
            kk = k - iy1(i)
            do j = ix1(i),ix2(i)
                pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
                temp(k,j-kk) = max(t_top, min(t_bot, temp(k,j-kk)+pert))
                !print *, k, j, pert
            enddo
        enddo
    endif
enddo

if (igeotherm .eq. 1) then
    !!  geotherm of a given age accross the box with variable age
    cond_c = 2.3
    cond_m = 3.3
    dens_c = 2800.
    dens_m = 3300.
    pi = 3.14159
    diffusivity = 1.e-6

    shf = 70.
    huc = 6.
    hlc = 2.
    F = 0.74 ! Partition coefficient
    alc = 0.4e-6 ! LowerCrust HG
    am = 0.02e-6 ! Mantle HG

    !! Continental geotherm
    ! After HasterokChapman 2011
    do i = 1, nx
        kl = 0
        do j = 1,nz-1
            if (i.lt.nx) then
                xc = 0.25*(cord (j,i  ,1) + cord(j+1,i  ,1) + &
                        cord (j,i+1,1) + cord(j+1,i+1,1))
                yc = 0.25*(cord (j,i  ,2) + cord(j+1,i  ,2) + &
                        cord (j,i+1,2) + cord(j+1,i+1,2))
            endif
            xsfh  = shf + g_amplitude/1e3*exp(-((xc-g_x0)/g_width)**2.)
            auc = (1. - F)*xsfh*1.e-3/(huc*1000.) ! UpperCrust heat generation
            ! Bootstrap through temps w depth
            temp(1,i) = t_top 
            Q(1,i) = xsfh*1.e-3 
            
            ! Depth in km 
            y = (cord(1,i,2)-cord(j,i,2))*1.e-3
            ! Distance between elements
            xdz = abs(cord(j+1,i,2)-cord(j,i,2))
            
            if (itemp_bc.eq.3) then
                zPotT = bot_bc*exp(9.81*y*0.00004)
            else 
                zPotT = 1300.*exp(9.81*y*0.00004)
            end if

            ! Upper crust
            if (y.lt.huc) then
                temp(j+1,i) = temp(j,i) + (Q(j,i)/cond_c)*(xdz) - (auc)/(2.*cond_c)*(xdz**2)
                Q(j+1,i) = Q(j,i) - auc*(xdz)
                if (i .lt. nx) source(j,i) = auc*1.e-3
            endif

            ! Lower crust
            if (y.ge.huc.and.y.lt.(hlc+huc)) then
                temp(j+1,i) = temp(j,i) + (Q(j,i)/cond_c)*(xdz) - (alc)/(2.*cond_c)*(xdz**2)
                Q(j+1,i) = Q(j,i) - alc*(xdz)
                if (i .lt. nx) source(j,i) = alc*1.e-3
            endif
                
                ! Mantle Lithosphere
            if (y.ge.(hlc+huc).and.kl.eq.0) then
                temp(j+1,i) = temp(j,i) + (Q(j,i)/cond_m)*(xdz) - (am)/(2.*cond_m)*(xdz**2)
                Q(j+1,i) = Q(j,i) - am*(xdz) 
                if (i .lt. nx) source(j,i) = am*1.e-3
                if(temp(j+1,i).ge.1300.) then
                    if (i .lt. nx) then
                        iphase(j,i) = kmant2
                        call newphase2marker(j,j,i,i,kmant2)
                    endif
                endif
                if(temp(j+1,i).ge.zPotT) then
                        kl =1
                        dum = 0
                endif
            endif

                ! Asthenosphere
            if (kl.eq.1) then
                y = (cord(1,i,2)-cord(j+1,i,2))*1.e-3
                
                if (itemp_bc.eq.3) then
                        zPotT=bot_bc*exp(9.81*y*0.00004)
                else
                        zPotT=1300.*exp(9.81*y*0.00004)
                endif

                temp(j+1,i) = zPotT
                if (i .lt. nx) then
                    iphase(j,i) = kmant2
                    call newphase2marker(j,j,i,i,kmant2)
                endif
                if (i .lt. nx) source(j,i) = 0.

                if (dum.eq.0) then
                        xlab(i) = y
                        dum = 1
                end if
            endif          
        enddo
    enddo
endif

!call RedefineTemp

return
end subroutine init_temp


subroutine init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
!$ACC routine vector
use arrays
use params
include 'precision.inc'

double precision, parameter :: pi = 3.14159265358979323846d0


if (ictherm(n)==1) then
    !! Oceanic geotherm (half space cooling model, T&S 3rd ed. Eq(4.113))
    diffusivity = 1.d-6
    age = age_1n
    !$ACC loop
    do j = 1,nz
        f = (cord(1,i,2)-cord(j,i,2)) / sqrt(4 * diffusivity * age * 1.d6 * sec_year)
        temp(j,i) = t_top + (t_bot - t_top) * erf(f)
        !print *, j, age, -cord(j,i,2), temp(j,i)
    enddo
elseif (ictherm(n)==2) then
    !! Oceanic geotherm (plate cooling cooling model, T&S 3rd ed. Eq(4.130))
    diffusivity = 1.d-6
    age = age_1n
    yL0 = tp1n    ! plate thickness in km
    age_init = age*sec_year*1d6
    tau_d = yL0*yL0*1d6 / (pi*pi*diffusivity)
    !$ACC loop
    do j = 1,nz
        ! depth in km
        y = (cord(1,i,2) - cord(j,i,2)) * 1d-3
        tss = t_top + (t_bot - t_top) * y / yL0
        tt = 0.d0
        do k = 1,100
            tt = tt + 1.d0/k * exp(-k*k*age_init/tau_d) * sin(pi*k*y/yL0)
        enddo
        temp(j,i) = tss + 2.d0/pi*(t_bot-t_top)*tt
        if(temp(j,i)>t_bot .or. y>yL0) temp(j,i) = t_bot
    enddo
elseif (ictherm(n)==12) then
    !! Continental geotherm (plate cooling model with radiogenic heating)
    !
    ! Starting from the steady state (ss) solution as in T&S 3rd ed. Eq(4.30)
    ! Let the ss moho temperature be tm and ss heatflux be qm.
    ! qm = cond * (t_bot - tm) / (yL0 - ymoho)
    ! Substitute qm to 4.30 to solve for tm
    age = age_1n
    yL0 = tp1n    ! plate thickness in km
    ymoho = tp2n  ! crust thickness in km
    cond = 3.3d0
    dens_c = 2700.d0
    diffusivity = 1.d-6
    !   write(*,*) rzbo, hs, hr
    age_init = age*sec_year*1d6
    tau_d = yL0*yL0*1d6 / (pi*pi*diffusivity)
    rr = ymoho / (yL0 - ymoho)
    tm = (t_top + rr*t_bot + dens_c*hs*hr*hr*1d6/cond*(1d0-exp(-ymoho/hr))) / (1 + rr)
    qm = cond * (t_bot - tm) / (yL0 - ymoho)
    !$ACC loop
    do j = 1,nz
        ! depth in km
        y = (cord(1,i,2) - cord(j,i,2)) * 1d-3

        ! ss part with radiogenic heat
        if (y <= ymoho) then
            tss = t_top + qm/cond*y + (dens_c*hs*hr*hr*1.d+6/cond)*(1d0-exp(-exp(-y/hr)))
        elseif (y <= yL0) then ! below moho, inside lithosphere
            tss = tm + qm/cond*(y-ymoho)
        else
            tss = t_bot
        endif
        ! time-dependent part
        ! see T&S 3rd ed. Eq(4.130)
        tt = 0.d0
        do k = 1,100
            tt = tt + 1.d0/k * exp(-k*k*age_init/tau_d) * sin(pi*k*y/yL0)
        enddo
        temp(j,i) = tss + 2.d0/pi*(t_bot-t_top)*tt
        if(temp(j,i)>t_bot .or. y>yL0) temp(j,i) = t_bot
        !       write(*,*) j,y,tss,yL0,tt
    enddo
elseif (ictherm(n)==21) then
    !! Constant geotherm gradient at top layer, then T=t_bot all the way to the bottom
    bot_dep = age_1n   ! bottom of the top layer
    do j = 1,nz
        y = (cord(1,i,2)-cord(j,i,2))*1.d-3
        temp(j,i) = y * (t_bot-t_top) / bot_dep
        if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
    enddo
elseif (ictherm(n)==22) then
    !! Constant geotherm gradient at top two layers, then T=t_bot all the way to the bottom
    temp1 = age_1n
    ylayer1 = tp1n  ! depth in km
    ylayer2 = tp2n  ! depth in km
    !$ACC loop
    do j = 1,nz
        y = (cord(1,i,2)-cord(j,i,2))*1.d-3
        if (y <= ylayer1) then
            temp(j,i) = t_top + (temp1 - t_top) * y / ylayer1
        else
            temp(j,i) = temp1 + (t_bot - temp1) * (y - ylayer1) / (ylayer2 - ylayer1)
        endif
        if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
    enddo
else
    !call sysmsg('init_temp: ictherm not supported!')
    stop 1
endif

return
end subroutine init_geotherm_profile


subroutine sidewalltemp(i1, i2)
!$ACC routine(init_geotherm_profile) vector

use arrays, only : temp, cord, source
use params
implicit none
! This subroutine is intended for remeshing.

integer, intent(in) :: i1, i2
double precision :: age_1n, tp1n, tp2n
integer :: n, i, j, k
double precision :: cond_c, cond_m, dens_c, dens_m, pi, diffusivity
double precision :: shf, huc, hlc
double precision :: tr, q_m, tm, age_init, diff_m, tau_d, y, tss, tt, pp, an
double precision :: xhc, age_1_new

cond_c = 2.2
cond_m = 3.3
dens_c = 2700.
dens_m = 3300.
pi = 3.14159
diffusivity = 1.e-6
shf = 70.
huc = 6.
hlc = 2.
age_1_new = 0.

xhc = 0.

if(i1 == 1) then
    ! left sidewall
    n = 1
else
    ! right sidewall
    n = nzone_age
endif
age_1n = age_1(n)
tp1n = tp1(n)
tp2n = tp2(n)

if (igeotherm .eq. 1) then
    tr= dens_c*hs*hr*hr*1.e+6/cond_c*exp(1.-exp(-(huc+hlc)/hr))
    q_m = (t_bot-t_top-tr)/(((huc+hlc)*1000.)/cond_c+((200.e3-(xhc)*1000.))/cond_m)
    tm  = t_top + (q_m/cond_c)*(huc+hlc)*1000. + tr
    age_init = age_1_new*3.14*1.e+7*1.e+6 + time
    diff_m = cond_m/1000./dens_m
    tau_d = 200.e3*200.e3/(pi*pi*diff_m)

    do i = i1, i2
        do j = 1,nz
            ! depth in km
            y = (cord(1,i,2)-cord(j,i,2))*1.e-3
            !  steady state part
            if (y.le.(huc+hlc)) then
                tss = t_top+(q_m/cond_c)*y*1000.+(dens_c*hs*hr*hr*1.e+6/cond_c)*exp(1.-exp(-y/hr))
            else
                tss = tm + (q_m/cond_m)*1000.*(y-xhc)
            endif

            ! time-dependent part
            tt = 0.
            pp =-1.
            do k = 1,100
                an = 1.*k
                pp = -pp
                tt = tt +pp/(an)*exp(-an*an*age_init/tau_d)*dsin(pi*k*(200.e3-y*1000.)/(200.e3))
            enddo
            temp(j,i) = tss +2./pi*(t_bot-t_top)*tt
            if(temp(j,i).gt.1330.or.y.gt.200.) temp(j,i)= 1330.
            if (j.eq.1) temp(j,i) = t_top
        enddo
    enddo
else
!$ACC parallel loop async(1)
    do i = i1, i2
        call init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
    enddo
endif

if(i1 == 1) then
    !$ACC parallel loop async(1)
    do i = i1, i2-1
        source(1:nz-1,i) = source(1:nz-1,i2)
    enddo
else
    !$ACC parallel loop async(1)
    do i = i1, i2-1
        source(1:nz-1,i) = source(1:nz-1,i1-1)
    enddo
endif

return
end subroutine sidewalltemp
