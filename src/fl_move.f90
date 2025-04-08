
! Move grid and adjust stresses due to rotation

subroutine fl_move
use arrays
use params
include 'precision.inc'

double precision :: dh(nx)

! Move Grid
if (movegrid .eq. 0) return

! UPDATING COORDINATES

!$OMP parallel private(i)
!$OMP do
!$ACC parallel loop async(1)
do i = 1,nx
!    write(*,*) cord(j,i,1),cord(j,i,2),vel(j,i,1),vel(j,i,2),dt
    cord(:,i,1) = cord(:,i,1) + vel(:,i,1)*dt
    cord(:,i,2) = cord(:,i,2) + vel(:,i,2)*dt
!    write(*,*) cord(j,i,1),cord(j,i,2)
enddo
!$OMP end do
!$OMP end parallel

dh = 0.
! Diffuse topography
if( topo_kappa.gt.0.d0) call diff_topo(dh)

if (itype_melting .eq. 1) then
    call arc_extru(dh)
else if (itype_melting .eq. 2) then
    call mor_melting(dh)
end if

! adjust markers
if (topo_kappa > 0 .or. itype_melting .ge. 1) then
    call correct_surface_marker(dh)

    if(mod(nloop, ifreq_avgsr) .eq. 0) then
        call resurface
    end if
endif


!$OMP parallel private(i,j,x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  det,dw12,s11,s22,s12)
!$OMP do
!--- Adjusting Stresses And Updating Areas Of Elements
!$ACC parallel loop collapse(2) async(1)
do  i = 1,nx-1
    do  j = 1,nz-1

        ! Coordinates
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)

        ! Velocities
        vx1 = vel (j  ,i  ,1)
        vy1 = vel (j  ,i  ,2)
        vx2 = vel (j+1,i  ,1)
        vy2 = vel (j+1,i  ,2)
        vx3 = vel (j  ,i+1,1)
        vy3 = vel (j  ,i+1,2)
        vx4 = vel (j+1,i+1,1)
        vy4 = vel (j+1,i+1,2)

        ! (1) Element A:
        det=((x2*y3-y2*x3)-(x1*y3-y1*x3)+(x1*y2-y1*x2))
        dvol(j,i,1) = det*area(j,i,1) - 1
        area(j,i,1) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx1*(x3-x2)+vx2*(x1-x3)+vx3*(x2-x1) - &
            vy1*(y2-y3)-vy2*(y3-y1)-vy3*(y1-y2))/det*dt
        s11 = stress0(j,i,1,1)
        s22 = stress0(j,i,2,1)
        s12 = stress0(j,i,3,1)
        stress0(j,i,1,1) = s11 + s12*2*dw12
        stress0(j,i,2,1) = s22 - s12*2*dw12
        stress0(j,i,3,1) = s12 + dw12*(s22-s11)

        ! rotate strains 
        s11 = strain(j,i,1)
        s22 = strain(j,i,2)
        s12 = strain(j,i,3)
        strain(j,i,1) = s11 + s12*2*dw12
        strain(j,i,2) = s22 - s12*2*dw12
        strain(j,i,3) = s12 + dw12*(s22-s11)

        ! (2) Element B:
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        dvol(j,i,2) = det*area(j,i,2) - 1
        area(j,i,2) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx3*(x4-x2)+vx2*(x3-x4)+vx4*(x2-x3) - &
           vy3*(y2-y4)-vy2*(y4-y3)-vy4*(y3-y2))/det*dt
        s11 = stress0(j,i,1,2)
        s22 = stress0(j,i,2,2)
        s12 = stress0(j,i,3,2)
        stress0(j,i,1,2) = s11 + s12*2*dw12
        stress0(j,i,2,2) = s22 - s12*2*dw12
        stress0(j,i,3,2) = s12 + dw12*(s22-s11)

        ! (3) Element C:
        det=((x2*y4-y2*x4)-(x1*y4-y1*x4)+(x1*y2-y1*x2))
        dvol(j,i,3) = det*area(j,i,3) - 1
        area(j,i,3) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx1*(x4-x2)+vx2*(x1-x4)+vx4*(x2-x1) - &
           vy1*(y2-y4)-vy2*(y4-y1)-vy4*(y1-y2))/det*dt
        s11 = stress0(j,i,1,3)
        s22 = stress0(j,i,2,3)
        s12 = stress0(j,i,3,3)
        stress0(j,i,1,3) = s11 + s12*2*dw12
        stress0(j,i,2,3) = s22 - s12*2*dw12
        stress0(j,i,3,3) = s12 + dw12*(s22-s11)

        ! (4) Element D:
        det=((x4*y3-y4*x3)-(x1*y3-y1*x3)+(x1*y4-y1*x4))
        dvol(j,i,4) = det*area(j,i,4) - 1
        area(j,i,4) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx1*(x3-x4)+vx4*(x1-x3)+vx3*(x4-x1) - &
            vy1*(y4-y3)-vy4*(y3-y1)-vy3*(y1-y4))/det*dt
        s11 = stress0(j,i,1,4)
        s22 = stress0(j,i,2,4)
        s12 = stress0(j,i,3,4)
        stress0(j,i,1,4) = s11 + s12*2*dw12
        stress0(j,i,2,4) = s22 - s12*2*dw12
        stress0(j,i,3,4) = s12 + dw12*(s22-s11)

        if (any(area(j,i,:) <= 0)) then
            ! write(333, *) 'area', j, i, nloop
            ! write(333, *) area(j,i,:)
            ! write(333, *) 'cord:'
            ! write(333, *) cord(j  ,i  ,:)
            ! write(333, *) cord(j  ,i+1,:)
            ! write(333, *) cord(j+1,i  ,:)
            ! write(333, *) cord(j+1,i+1,:)
            ! write(333, *) 'vel:'
            ! write(333, *) vel(j  ,i  ,:)
            ! write(333, *) vel(j  ,i+1,:)
            ! write(333, *) vel(j+1,i  ,:)
            ! write(333, *) vel(j+1,i+1,:)
            ! flush(333)
            ! call SysMsg('Negative area!')
            stop 40
        endif
    enddo
enddo
!$OMP end do
!$OMP end parallel

return
end subroutine fl_move

subroutine check_chamber
    use arrays
    use params
    use phases
    implicit none
    double precision :: critmeltfrac(nz-1,nx-1), tmpr, critmelt
    integer :: i,j,k,ibug,imagcolumn(nx-1),jmagcolumn(nx-1),jmagtop,magtopcount,imag_min,imag_max


    imagtop = int((nx-1)/2)
    imagcolumn(:) = 0
    jmagcolumn(:) = nz-1
    ibug = 0
    
    critmeltfrac = 0.
    
    ! find top magma chambers
    do i = 1, nx-1
        k = jmoho(i)    
        do j = k , nz-1
            tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
            if (tmpr.ge.1100.) then
                critmeltfrac(j,i) = Eff_melt(j,i)
            endif
        enddo
    enddo
    
    critmelt = maxval(critmeltfrac)
    do i = 1, nx-1
        k = jmoho(i)
        do j = k , nz-1
            if (critmeltfrac(j,i) .ge. critmelt) then
                jmagcolumn(i) = j
                imagcolumn(i) = i
                ibug = ibug + 1
                cycle
            endif
        enddo
    enddo
    
    jmagtop = minval(jmagcolumn)
    magtopcount = 0
    if (ibug .gt. 0 ) then
        do i = 1, nx-1
            if (jmagcolumn(i) .eq. jmagtop.and.magtopcount.eq.0) then
                magtopcount = magtopcount + 1
                imag_min = imagcolumn(i)
            endif
            if (jmagcolumn(i) .eq. jmagtop.and. magtopcount .gt. 0) then
                magtopcount = magtopcount + 1
                imag_max = imagcolumn(i)
            endif
        enddo
        if (imag_max.eq.0) then
            imagtop = imag_min
        else
            imagtop = int(0.5*(imag_max+imag_min))
        endif
    endif
    
return
end subroutine check_chamber


subroutine mor_melting(dh)
use arrays
use params
use phases
use marker_data
implicit none

double precision, external :: dlmin_prop
double precision :: totalmelt, totalarea, dh(nx)
double precision :: mor_extrusion_rate, mor_dike_rate, pi, dlmin
double precision :: control_vol_ch, vc_rate, xintr, xmu, xsigma, vol_ch, str_ch, vol_tot
double precision :: unit_xmelt_migrated, xl_vol, quad_area, xdike_migrated
double precision :: extru_limit, intru_unit_limit, total_extru_strain, tmp, extru_melt
integer :: i, j, ii, jj, kinc, n_to_add, kk, ihalfwidth_mzone
integer :: iextru_start, iextru_end, nintru, iintru_bot, ibasement
    
dlmin = dlmin_prop()

av_intrusion = 0.
new_intrusion = 0.
ii = imagtop
totalmelt = 0.
totalarea = 0.
ibasement = 2

mor_extrusion_rate = 1.d0 - ratio_crust_mzone - ratio_mantle_mzone
mor_dike_rate = ratio_crust_mzone

iintru_bot = jmoho(ii) - 1
nintru = iintru_bot - ibasement + 1

ihalfwidth_mzone = int(width_mzone / 2 / dxmin)
iextru_start = max(1,ii-2*ihalfwidth_mzone)
iextru_end = min(nx-1,ii+2*ihalfwidth_mzone)

!$OMP parallel do private(i,j) collapse(2)
do i = 1, nx-1
    do j = 1, nz-1
        dummye(j,i) = 0.5d0/area(j,i,1) + 0.5d0/area(j,i,2)
    enddo
enddo

! volume of the melt in this column
!$OMP parallel do private(i,j,quad_area) reduction(+:totalmelt,totalarea) collapse(2)
do i = iextru_start, iextru_end
    do j = 1, nz-1
        if (Eff_melt(j,i).gt.0.0) then ! only the melt accumulated
            quad_area = dummye(j,i)
            totalarea = totalarea + quad_area
            if (fmelt(j,i).gt.0. .and. Eff_melt(j,i).gt.0.03) then  ! only the melt produced that can move
                totalmelt = totalmelt + quad_area * fmelt(j,i)
            endif
        endif
    enddo
enddo

pi = sqrt(2.*3.14159265358979323846)

! calculate the volume change and magma for extrusives
! We need to let the melt spread across the basin (we have to be symmetric)
! we impose a max of 2 particles per element for the volume change
! We use the same width has where the melt is collected

intru_unit_limit = 1.*2.*vbc*dlmin ! per element
! 1/4 of the total volume of intrusion
extru_limit = nintru * intru_unit_limit * mor_extrusion_rate / mor_dike_rate

!Diking 
! calculate the volume change for each dike elements
! The dike does not get in the extrusives

xintr = mor_dike_rate * totalmelt / nintru
control_vol_ch = intru_unit_limit * dt  ! in m^2

!$OMP parallel do private(jj,vol_tot,vc_rate,vol_ch,str_ch)
do jj = ibasement, iintru_bot
    av_intrusion(jj,ii) = xintr ! in m^2
    vol_tot = av_intrusion(jj,ii) + stored_intrusion(jj,ii) ! in m^2
    vc_rate = vol_tot / dt ! in m^2/s

    if (vc_rate.gt.intru_unit_limit) then  ! in m^2/s
        vol_ch = control_vol_ch
        stored_intrusion(jj,ii) = vol_tot - vol_ch  ! in m^2
    else
        vol_ch = vol_tot
        stored_intrusion(jj,ii) = 0.  ! in m^2
    endif

    dv_intr(jj,ii) = dv_intr(jj,ii) + vol_ch

    str_ch = vol_ch / dummye(jj,ii) ! in volumic strain
    new_intrusion(jj,ii) = str_ch
    fmagma(jj,ii) = fmagma(jj,ii) + str_ch
enddo

xdike_migrated = 0.
! Calculate the amount of melt intruded
!$OMP parallel do private(jj) reduction(+:xdike_migrated)
do jj = ibasement, iintru_bot
    xdike_migrated = xdike_migrated + new_intrusion(jj,ii)
enddo

xdike_migrated = xdike_migrated / totalarea
!$OMP parallel do private(i,j) collapse(2)
do i = 1,nx-1
do j = 1,nz-1
    if (Eff_melt(j,i).ge.0.03) then
        Eff_melt(j,i) = Eff_melt(j,i) - xdike_migrated * dummye(j,i)
    endif
enddo
enddo

! add basalt in intruded element
!$OMP parallel do private(i,jj,kinc,xl_vol,quad_area,n_to_add,kk) collapse(2)
do i = ii-3, ii+3
    do jj = ibasement, iintru_bot
        kinc = nmark_elem(jj,i)
        xl_vol = dv_intr(jj,i)
        quad_area = dummye(jj,i)
        if (xl_vol * kinc >= quad_area .and. kinc .ne. max_markers_per_elem) then
            ! intrusion, add a mafic marker
            n_to_add = min(ceiling((xl_vol / quad_area) * kinc), max_markers_per_elem - kinc)
            do kk = 1, n_to_add
                call add_marker_dike(jj,i, 1d0, time, nloop+i+kk, kocean0)
            enddo
            dv_intr(jj,i) = 0.

            ! recalculate phase ratio
            call count_phase_ratio(jj,i)
        endif
    enddo
enddo



! REVISE  Need to be able to inject in more elements laterally
! Volcanic flow rate ? From sesismic 1e-6 m^2/s  500 km^2/4 myr
extru_melt = mor_extrusion_rate * totalmelt  ! in m^2
vc_rate = (extru_melt + stored_intrusion(1,ii)) / dt ! in m^2/s the full rate of extension rate
control_vol_ch = extru_limit * dt  ! in m^2
quad_area = dummye(1,ii)

if (vc_rate.gt.extru_limit) then
    stored_intrusion(1,ii) = stored_intrusion(1,ii) + extru_melt - control_vol_ch ! in m^2
    total_extru_strain = control_vol_ch / quad_area ! in volumic strain
else
    total_extru_strain = (extru_melt + stored_intrusion(1,ii)) / quad_area
    stored_intrusion(1,ii) = 0.
endif

xintr = total_extru_strain
xsigma = dfloat(ihalfwidth_mzone/2)
! Distribute average intrusion on a normal distribution
!$OMP parallel do private(i,xmu,tmp)
do i = iextru_start, iextru_end
    xmu = dfloat(ii - i)
    tmp = (xintr/pi/xsigma) * dexp(-(xmu)**2./2./xsigma**2.)
    av_intrusion(1,i) = tmp
    fmagma(1,i) = fmagma(1,i) + tmp

    ! height of extrusion in this column
    extrusion(i) = tmp * quad_area / (cord(1,i+1,1) - cord(1,i,1)) / 2.
    extr_acc(i) = extr_acc(i) + extrusion(i)
    !$ACC atomic update
    !$OMP atomic update
    cord(1,i,2) = cord(1,i,2) + extrusion(i)
    !$ACC atomic update
    !$OMP atomic update
    cord(1,i+1,2) = cord(1,i+1,2) + extrusion(i)

    !$ACC atomic update
    !$OMP atomic update
    dh(i) = dh(i) + extrusion(i)
    !$ACC atomic update
    !$OMP atomic update
    dh(i+1) = dh(i+1) + extrusion(i)
enddo

! Ratio of amount of melt removed the mantle per asthenosphere element
unit_xmelt_migrated = total_extru_strain / totalarea
! Remove this melt on average from the asthenosphere
!$OMP parallel do private(i,j) collapse(2)
do j = 1,nz-1
    do i = 1,nx-1
        if (Eff_melt(j,i).ge.0.03) then
            Eff_melt(j,i) = Eff_melt(j,i) - unit_xmelt_migrated * dummye(j,i)
        endif
    enddo
enddo
! add basalt in extrusion elements  

return
end subroutine mor_melting


subroutine arc_extru(dh)
use arrays
use params
include 'precision.inc'

double precision :: dh(nx)

! magma extrusion
arc_extrusion_rate = 1.d0 - ratio_mantle_mzone
if (arc_extrusion_rate > 0) then
    !$ACC parallel loop async(1)
    do i = 2, nx-2  ! avoid edge elements, which should not contain arc magma
        totalmelt = 0
        !$ACC loop reduction(+:totalmelt)
        do j = 1, nz-1
            quad_area = 0.5d0/area(j,i,1) + 0.5d0/area(j,i,2)
            ! volume of the melt in this column
            totalmelt = totalmelt + quad_area * fmelt(j,i)
        enddo
        ! height of extrusion in this column
        extrusion(i) = arc_extrusion_rate * dt * totalmelt * prod_magma &
            / (cord(1,i+1,1) - cord(1,i,1) + 0.5d0 * (cord(1,i,1) - cord(1,i-1,1) + cord(1,i+2,1) - cord(1,i+1,1)))
        !print *, i, extrusion(i), totalmelt
        extr_acc(i) = extr_acc(i) + extrusion(i)
        !$ACC atomic update
        cord(1,i,2) = cord(1,i,2) + extrusion(i)
        !$ACC atomic update
        cord(1,i+1,2) = cord(1,i+1,2) + extrusion(i)
        !$ACC atomic update
        dh(i) = dh(i) + extrusion(i)
        !$ACC atomic update
        dh(i+1) = dh(i+1) + extrusion(i)
    enddo
endif

return
end subroutine arc_extru


!============================================================
! Diffuse topography
!============================================================
subroutine diff_topo(dh)
use arrays
use params
use phases
include 'precision.inc'

double precision :: dh(nx)

!EROSION PROCESSES
if( topo_kappa .gt. 0.d0 ) then

    !$ACC parallel loop async(2)
    do i = 1, nx
        stmpn(i) = topo_kappa ! elevation-dep. topo diffusivity
    enddo

    topomean = 0
    !$ACC parallel loop reduction(+:topomean) async(1)
    do i = 1, nx
        topomean = topomean + cord(1,i,2) / nx
    enddo

    !$ACC wait(2)

    ! !$ACC parallel loop async(1)
    ! do i = 1, nx-1
    !     ! higher erosion for sediments above mean topo
    !     if (iphase(1,i) == ksed2 .and. cord(1,i,2) > topomean) stmpn(i) = stmpn(i) * 10
    ! enddo


    !$ACC parallel loop async(1)
    do i = 2, nx-1
        snder = ( stmpn(i+1)*(cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
            stmpn(i-1)*(cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
        dtopo(i) = dt * snder
    end do

    !$ACC serial async(1)
    dtopo(1) = dtopo(2)
    dtopo(nx) = dtopo(nx-1)
    !$ACC end serial

    !$ACC parallel loop async(1)
    do i = 1, nx
        ! erosion cannot erode over 0.5x element height
        dtopo(i) = min(dtopo(i), 0.5*(cord(1,i,2)-cord(2,i,2)))
    enddo

    ! accumulated topo change since last resurface
    !$ACC wait(1)
    !$ACC parallel loop async(2)
    do i = 1, nx-1
        dhacc(i) = dhacc(i) + 0.5d0 * (dtopo(i) + dtopo(i + 1))
    enddo

    !$ACC parallel loop async(1)
    do i = 1, nx
        cord(1,i,2) = cord(1,i,2) + dtopo(i)
        dh(i) = dh(i) + dtopo(i)
    enddo
endif

return
end subroutine diff_topo


subroutine resurface
!$ACC routine(bar2xy) seq
!$ACC routine(shape_functions) seq
!$ACC routine(add_marker_at_top) seq
use marker_data
use arrays
use params
use phases
include 'precision.inc'

!$ACC serial private(dz_ratio) async(1)
do i = 1, nx-1
    ! averge thickness of this element
    elz = 0.5d0 * (cord(1,i,2) - cord(2,i,2) + cord(1,i+1,2) - cord(2,i+1,2))
    ! change in topo
    chgtopo = dhacc(i)
    ! # of markers in this element
    kinc = nmark_elem(1,i)

    ichanged = 0


    if (chgtopo * kinc >= elz .and. kinc .ne. max_markers_per_elem) then
        ! sedimentation, add a sediment marker
        !write(6,*) 'sediment', i, chgtopo, elz
        n_to_add = min(ceiling(chgtopo / elz * kinc), max_markers_per_elem - kinc)
        dz_ratio = min(chgtopo / elz, 1.0d0)
        do ii = 1, n_to_add
            call add_marker_at_top(i, dz_ratio, time, nloop+i+ii, ksed2)
        enddo

        dhacc(i) = 0
        ichanged = 1
    endif

    ! change in topo due to volcanism
    chgtopo2 = extr_acc(i)
    if (chgtopo2 * kinc >= elz .and. kinc .ne. max_markers_per_elem) then
        ! extrusion, add an arc marker
        n_to_add = min(ceiling(chgtopo2 / elz * kinc), max_markers_per_elem - kinc)
        dz_ratio = min(chgtopo2 / elz, 1.0d0)
        !write(6,*) 'arc', i, chgtopo2, elz, n_to_add, dz_ratio
        if (itype_melting .eq. 1) then
            kind = karc1
        else if (itype_melting .eq. 2) then
            kind = kocean2
        endif

        do ii = 1, n_to_add
            call add_marker_at_top(i, dz_ratio, time, nloop+i+ii, kind)
        enddo

        extr_acc(i) = 0
        ichanged = 1
    endif

    if (ichanged == 1) then
        ! recalculate phase ratio
        call count_phase_ratio(1,i)
    endif
end do
!$ACC end serial
!$ACC update self(nmarkers) async(1)

return
end subroutine resurface


subroutine correct_surface_marker(dh)
!$ACC routine seq
!$ACC routine(euler2bar) seq
use marker_data
use params
use arrays
include 'precision.inc'

double precision :: dh(nx), ichanged(nx-1)
double precision :: bar(3), xx(3), yy(3)

do i = 1, nx-1
    ichanged(i) = 0
    dh1 = dh(i)
    dh2 = dh(i+1)
    k = 1
    do while (k .le. nmark_elem(1, i))
        n = mark_id_elem(k,1,i)

        if ( n .eq. 0 .or. k .gt. nmark_elem(1, i)) then
            print*, 'something wrong in correct top marker'
            exit
        end if
        ntr = mark_ntriag(n)
        bar(1) = mark_a1(n)
        bar(2) = mark_a2(n)
        bar(3) = 1.d0 - bar(1) -bar(2)

        kk = MOD(ntr - 1, 2) + 1
        jj = MOD((ntr - kk) / 2, nz-1) + 1
        ii = (ntr - kk) / 2 / (nz - 1) + 1

        if (kk .eq. 1) then
            xx(1) = cord(1 ,i  ,1)
            xx(2) = cord(2 ,i  ,1)
            xx(3) = cord(1 ,i+1,1)
            yy(1) = cord(1 ,i  ,2) -dh1
            yy(2) = cord(2 ,i  ,2)
            yy(3) = cord(1 ,i+1,2) -dh2
        else
            xx(1) = cord(1 ,i+1,1)
            xx(2) = cord(2 ,i  ,1)
            xx(3) = cord(2 ,i+1,1)
            yy(1) = cord(1 ,i+1,2) -dh2
            yy(2) = cord(2 ,i  ,2)
            yy(3) = cord(2 ,i+1,2)
        endif

        xxx=sum(bar*xx)
        yyy=sum(bar*yy)
        jtop = 1
        itop = i
        call euler2bar(xxx,yyy,bar(1),bar(2),ntr,itop,jtop,inc)

        ! marker out of the mesh -> remover the marker
        if (inc .eq. 0) then
            ! replace marker k with the last marker
            mark_id_elem(k, 1, i) = mark_id_elem(nmark_elem(1, i), 1, i)
            mark_id_elem(nmark_elem(1, i), 1, i) = 0
            ! delete marker
            mark_dead(n) = 0
            nmark_elem(1, i) = nmark_elem(1, i) - 1
            ichanged(i) = 1
        else
            mark_a1(n) = bar(1)
            mark_a2(n) = bar(2)
            mark_ntriag(n) = ntr

            if (itop .eq. i .and. jtop .eq. 1) then
                k = k + 1
            else
                kinc_next = nmark_elem(1, itop)
                if (kinc_next < max_markers_per_elem) then
                    nmark_elem(1, itop) = nmark_elem(1, itop) + 1
                    mark_id_elem(kinc_next+1, 1, itop) = n
                else
                    mark_dead(n) = 0
                end if

                mark_id_elem(k, 1, i) = mark_id_elem(nmark_elem(1, i), 1, i)
                mark_id_elem(nmark_elem(1, i), 1, i) = 0
                nmark_elem(1, i) = nmark_elem(1, i) - 1
                ichanged(i) = 1
            end if
        end if
    end do
end do

! recalculate the phase ratio
!$OMP parallel do private(i)
do i = 1, nx-1
    if (ichanged(i) .eq. 1) call count_phase_ratio(1,i)
end do

return
end subroutine correct_surface_marker
    

subroutine add_marker_at_top(i, dz_ratio, time, loop, kph)
!$ACC routine seq
!$ACC routine(add_marker) seq
use myrandom_mod
use marker_data
use arrays
include 'precision.inc'

iseed = loop + i
icount = 0
do while(.true.)
    call myrandom(iseed, r1)
    call myrandom(iseed, r2)
    j = 1

    ! (x1, y1) and (x2, y2)
    x1 = cord(j  ,i,1)*(1-r1) + cord(j  ,i+1,1)*r1
    y1 = cord(j  ,i,2)*(1-r1) + cord(j  ,i+1,2)*r1
    x2 = cord(j+1,i,1)*(1-r1) + cord(j+1,i+1,1)*r1
    y2 = cord(j+1,i,2)*(1-r1) + cord(j+1,i+1,2)*r1

    ! connect the above two points
    ! (this point is not uniformly distributed within the element area
    ! and is biased against the thicker side of the element, but this
    ! point is almost gauranteed to be inside the element)
    r2 = r2 * dz_ratio
    xx = x1*(1-r2) + x2*r2
    yy = y1*(1-r2) + y2*r2

    call add_marker(xx, yy, kph, zpressm(j,i), Eff_melt(j,i), time, 1, i, inc)
    if(inc==1 .or. inc==-1) exit
    icount = icount + 1
    if(icount > 100) stop 134
    !write(333,*) 'add_marker_at_top failed: ', xx, yy, rx, elz, kph
    !write(333,*) '  ', cord(1,i,:)
    !write(333,*) '  ', cord(1,i+1,:)
    !write(333,*) '  ', cord(2,i,:)
    !write(333,*) '  ', cord(2,i+1,:)
    !call SysMsg('Cannot add marker.')
end do

return
end subroutine add_marker_at_top


subroutine add_marker_dike(j,i, vol_ratio, time, loop, kph)
!$ACC routine seq
!$ACC routine(add_marker) seq
use myrandom_mod
use marker_data
use arrays
include 'precision.inc'    


iseed = loop + i
icount = 0
do while(.true.)
    call myrandom(iseed, r1)
    call myrandom(iseed, r2)


    ! (x1, y1) and (x2, y2)
    x1 = cord(j  ,i,1)*(1-r1) + cord(j  ,i+1,1)*r1
    y1 = cord(j  ,i,2)*(1-r1) + cord(j  ,i+1,2)*r1
    x2 = cord(j+1,i,1)*(1-r1) + cord(j+1,i+1,1)*r1
    y2 = cord(j+1,i,2)*(1-r1) + cord(j+1,i+1,2)*r1

    ! connect the above two points
    ! (this point is not uniformly distributed within the element area
    ! and is biased against the thicker side of the element, but this
    ! point is almost gauranteed to be inside the element)
    r2 = r2 * vol_ratio
    xx = x1*(1-r2) + x2*r2
    yy = y1*(1-r2) + y2*r2

    call add_marker(xx, yy, kph, zpressm(j,i), Eff_melt(j,i), time, j, i, inc)
    if(inc==1 .or. inc==-1) exit
    icount = icount + 1
    write(*,*) inc,icount
    if(icount > 100) stop 134
end do

return
end subroutine add_marker_dike
