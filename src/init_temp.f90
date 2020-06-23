! Initiate temperature profile

subroutine init_temp
use arrays
use params
use phases
implicit none

integer :: i, j, i1, i2, iwidth, ixc, k, kk, n
double precision :: xc, yc, geoth, age, amp, &
                    cond_c, cond_m, dens_c, dens_m, diffusivity, &
                    pert, pert2, pi, y

!  Read distribution of temperatures from the dat file
if (irtemp .gt. 0) then
    open( 1, file=tempfile, status='old', err=101 )
    do i = 1, nx
    do j = 1, nz
        read( 1, * ,err=102 ) temp(j,i)
!     if(temp(j,i).ge.1000.d0) temp(j,i) = 1000.d0
    enddo
    enddo
    close(1)

   goto 10
    101 call SysMsg('INIT_TEMP: Cannot open file with temperatures initial distrib!')
    stop 21
    102 call SysMsg('INIT_TEMP: Error reading file with temperatures initial distrib!')
    stop 21
endif

select case(iynts)
case (1)
    call sysmsg('init_temp: iynts=1 is not supported anymore!')
    stop 1
case (2)
    !!  geotherm of a given age accross the box with variable age
    cond_c = 2.2d0
    cond_m = 3.3d0
    dens_c = 2700.d0
    dens_m = 3300.d0
    pi = 3.14159d0
    diffusivity = 1.d-6
    do n = 1, nzone_age
        if (n /= 1) then
            if (iph_col_trans(n-1) == 1) cycle
        endif

!!$        if(iph_col1(n)==kocean1 .or. iph_col1(n)==kocean2   &
!!$            .or. iph_col2(n)==kocean1 .or. iph_col2(n)==kocean2  &
!!$            .or. iph_col3(n)==kocean1 .or. iph_col3(n)==kocean2) then
            !! Oceanic geotherm (half space cooling model)
            do i = ixtb1(n), ixtb2(n)
                age = age_1(n)
                if (iph_col_trans(n) == 1) then
                    i1 = ixtb1(n)
                    i2 = ixtb2(n)
                    age = age_1(n) + (age_1(n+1) - age_1(n)) * (cord(1,i,1) - cord(1,i1,1)) / (cord(1,i2,1) - cord(1,i1,1))
                endif
                do j = 1,nz
                    ! depth in km
                    y = (cord(1,i,2)-cord(j,i,2)) / sqrt(4 * diffusivity * age * 1.d6 * sec_year)
                    temp(j,i) = t_top + (t_bot - t_top) * erf(y)
                    !print *, j, age, -cord(j,i,2), temp(j,i)
                enddo
            enddo
!!$        else
!!$            !! Continental geotherm
!!$            tr= dens_c*hs*hr*hr*1.d+6/cond_c*exp(1.d0-exp(-hc3(n)/hr))
!!$            q_m = (t_bot-t_top-tr)/((hc3(n)*1000.d0)/cond_c+((200.d3-(hc3(n))*1000.d0))/cond_m)
!!$            tm  = t_top + (q_m/cond_c)*hc3(n)*1000.d0 + tr
!!$            !   write(*,*) rzbo, tr, hs, hr, hc3(n), q_m, tm
!!$            diff_m = cond_m/1000.d0/dens_m
!!$            tau_d = 200.d3*200.d3/(pi*pi*diff_m)
!!$            do i = ixtb1(n), ixtb2(n)
!!$                age = age_1(n)
!!$                if (iph_col_trans(n) == 1) then
!!$                    i1 = ixtb1(n)
!!$                    i2 = ixtb2(n)
!!$                    age = age_1(n) + (age_1(n+1) - age_1(n)) * (cord(1,i,1) - cord(1,i1,1)) / (cord(1,i2,1) - cord(1,i1,1))
!!$                endif
!!$                age_init = age*3.14d0*1.d+7*1.d+6
!!$                do j = 1,nz
!!$                    ! depth in km
!!$                    y = (cord(1,i,2)-cord(j,i,2))*1.d-3
!!$                    !  steady state part
!!$                    if (y.le.hc3(n)) tss = t_top+(q_m/cond_c)*y*1000.d0+(dens_c*hs*hr*hr*1.d+6/cond_c)*exp(1.d0-exp(-y/hr))
!!$                    if (y.gt.hc3(n)) tss = tm + (q_m/cond_m)*1000.d0*(y-hc3(n))
!!$
!!$                    ! time-dependent part
!!$                    tt = 0.d0
!!$                    pp =-1.d0
!!$                    do k = 1,100
!!$                        an = 1.d0*k
!!$                        pp = -pp
!!$                        tt = tt +pp/(an)*exp(-an*an*age_init/tau_d)*dsin(pi*k*(200.d3-y*1000.d0)/(200.d3))
!!$                    enddo
!!$                    temp(j,i) = tss +2.d0/pi*(t_bot-t_top)*tt
!!$                    if(temp(j,i).gt.1330.or.y.gt.200.d0) temp(j,i)= 1330.d0
!!$                    if (j.eq.1) temp(j,i) = t_top
!!$                    !       write(*,*) tss,tm,q_m,cond_m,hc3(n),y,tt
!!$                enddo
!!$            enddo
!!$        endif
    enddo

case default
    ! estimate initial temperature as linear (for first approx. of conductivities)
    do j = 1,nz
        temp(j,1:nx) = (t_bot-t_top)/abs(rzbo)*abs(cord(j,1,2)-z0) + t_top
    end do
end select

10 continue

! DISTRIBUTE SOURCES in elements
do j = 1,nz-1
    y = -( cord(j+1,1,2)+cord(j,1,2) )/2 / 1000
    source(j,1:nx-1) = hs*exp(-y/hr)
end do

! Initial rectangular temperature perturbation
if( temp_per.ne.0.d0 ) then
    temp(iy1t:iy2t,ix1t:ix2t) = temp(iy1t:iy2t,ix1t:ix2t) + temp_per
endif              

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

return
end subroutine init_temp


subroutine sidewalltemp(i1, i2)
  !$ACC routine seq
  use arrays, only : temp, cord, source
  use params
  use phases
  implicit none

  integer :: i1, i2, n, i, j
  double precision :: cond_c, cond_m, dens_c, dens_m, pi, diffusivity, y

  ! This subroutine is intended for remeshing.
  cond_c = 2.2d0
  cond_m = 3.3d0
  dens_c = 2700.d0
  dens_m = 3300.d0
  pi = 3.14159d0
  diffusivity = 1.d-6

  if(nzone_age < 1) then
      stop 'nzone_age < 1, cannot determine temperature of incoming material'
  endif

  if(i1 == 1) then
      ! left sidewall
      n = 1
  else
      ! right sidewall
      n = nzone_age
  endif

!!$  if(iph_col3(n)==kocean1 .or. iph_col3(n)==kocean2) then
      !! Oceanic geotherm (half space cooling model)
      do i = i1, i2
          do j = 1,nz
              ! depth in km
              y = (cord(1,i,2)-cord(j,i,2)) / sqrt(4 * diffusivity * age_1(n) * 1.d6 * sec_year)
              temp(j,i) = t_top + (t_bot - t_top) * erf(y)
              !print *, j, age_1(n), -cord(j,i,2), temp(j,i)
          enddo
      enddo
!!$  else
!!$      !! Continental geotherm
!!$      tr= dens_c*hs*hr*hr*1.d+6/cond_c*exp(1.-exp(-hc3(n)/hr))
!!$      q_m = (t_bot-t_top-tr)/((hc3(n)*1000.d0)/cond_c+((200.d3-(hc3(n))*1000.d0))/cond_m)
!!$      tm  = t_top + (q_m/cond_c)*hc3(n)*1000.d0 + tr
!!$      !   write(*,*) rzbo, tr, hs, hr, hc3(n), q_m, tm
!!$      age_init = age_1(n)*3.14d0*1.d+7*1.d+6 + time
!!$      diff_m = cond_m/1000.d0/dens_m
!!$      tau_d = 200.d3*200.d3/(pi*pi*diff_m)
!!$
!!$      do i = i1, i2
!!$          do j = 1,nz
!!$              ! depth in km
!!$              y = (cord(1,i,2)-cord(j,i,2))*1.d-3
!!$              !  steady state part
!!$              if (y.le.hc3(n)) tss = t_top+(q_m/cond_c)*y*1000.d0+(dens_c*hs*hr*hr*1.d+6/cond_c)*exp(1.d0-exp(-y/hr))
!!$              if (y.gt.hc3(n)) tss = tm + (q_m/cond_m)*1000.d0*(y-hc3(n))
!!$
!!$              ! time-dependent part
!!$              tt = 0.d0
!!$              pp =-1.d0
!!$              do k = 1,100
!!$                  an = 1.d0*k
!!$                  pp = -pp
!!$                  tt = tt +pp/(an)*exp(-an*an*age_init/tau_d)*dsin(pi*k*(200.d3-y*1000.d0)/(200.d3))
!!$              enddo
!!$              temp(j,i) = tss +2.d0/pi*(t_bot-t_top)*tt
!!$              if(temp(j,i).gt.1330.or.y.gt.200.d0) temp(j,i)= 1330.d0
!!$              if (j.eq.1) temp(j,i) = t_top
!!$              !       write(*,*) tss,tm,q_m,cond_m,hc3(n),y,tt
!!$          enddo
!!$      enddo
!!$  endif

  if(i1 == 1) then
      do i = i1, i2
          source(1:nz-1,i) = source(1:nz-1,i2+1)
      enddo
  else
      do i = i1, i2
          source(1:nz-1,i) = source(1:nz-1,i1-1)
      enddo
  endif
  return
end subroutine sidewalltemp
