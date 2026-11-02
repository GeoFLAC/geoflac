
! Calculate thermal field in explict form

subroutine fl_therm
!$ACC routine(Eff_cp) seq
!$ACC routine(Eff_conduct) seq
!$ACC routine(add_marker) seq
!$ACC routine(count_phase_ratio) seq
use arrays
use params
use phases
use marker_data
implicit none
double precision :: D(3,3)  ! diffusion operator
double precision :: diff_elem(nz-1, nx-1)
double precision :: inv_area
double precision :: x1, y1, x2, y2, x3, y3, x4, y4
double precision :: shpdx_1, shpdx_2, shpdx_3, shpdz_1, shpdz_2, shpdz_3
integer :: i, j, iph, jm, ii, ihalf, ii_start, ii_end, ncol, jj, ihalfwidth_mzone
double precision :: tan_mzone, cp_eff, cond_eff, dissip, quad_area, fr_lambda, delta_fmagma, deltaT, real_area13, area_n, rhs, area_ratio, z_moho, z_melt, x_melt, h, x, z, tmpr, Eff_cp, Eff_conduct
double precision :: delta_fmagma2, deltaT1, deltaT2
double precision :: x_sum, z_sum, w_sum, w_barrier, weight, x_center, z_center
double precision :: x_wide_l, x_wide_r, w_l, w_r, h2, z_surf, tan_mzone2_l, tan_mzone2_r
integer :: i_center, j_center, ihalfwidth_mzone2_l, ihalfwidth_mzone2_r
logical :: found, mor_zone_valid
! dike-intrusion eigenstrain/marker bookkeeping (arc and MOR pathways)
integer :: jtop, kinc, n_to_add, kk, inc
double precision :: elem_area, pending, released, xc, yc, marker_vol_unit

tan_mzone = tan(0.5d0 * angle_mzone * 3.14159265358979323846d0 / 180.d0)
! max. width of the magma zone @ moho (as if melting occurs at 200 km)
ihalfwidth_mzone = ceiling(tan_mzone * 200e3 / dxmin)
w_barrier = 1.0d-4

! real_area = 0.5* (1./area(n,t))
! Calculate Fluxes in every triangle
!      flux (num_triangle, direction(x,y), j, i)
!
!  1 - 3
!  |   |
!  2 - 4
!
!  diagonal / :
!
!   A:        B:
!
!  1---3         1
!  | /         / |
!  2         2---3

! saving old temperature
!$ACC kernels async(1)
temp0(:,:) = temp(:,:)
!$ACC end kernels

if (itype_melting == 1) then
    !$OMP Parallel private(i,j,cp_eff,tmpr,fr_lambda,delta_fmagma,delta_fmagma2,deltaT,deltaT1,deltaT2)
    ! M: fmegma, magma fraction in the element
    ! dM/dt = P - M * fr_lambda
    ! Rearrange after forward Euler for dM/dt
    ! M(t+dt) = M(t) + P * dt - M(t) * fr_lambda * dt
    !   where P is the production rate and the second term is the freezing rate
    !   fr_lambda is temperature dependent: lower T freezes faster
    ! fr_lambda = lambda_freeze * exp(-lambda_freeze_tdep * (T - Ttop))
    ! P = fmelt * prod_magma * R_zone *  (A_zone / A_melt)
    !   fmelt: melt fraction, calculated in chage_phase()
    !   prod_magma: magma production rate, how rapid a melt migrate away
    !               from the active melting region.
    !   A_melt: the area of the melting element
    !   R_zone, A_zone: the magma distribution ratio and the area of the two
    !                   zones (crust and mantle). R_crust + R_mantle == 1.0
    !   area_ratio = A_zone / A_melt
    !
    ! Magma freezing will release latent heat to the element
    ! rho * cp * deltaT = rho * latent_heat * delta_fmagma
    ! Rearrage to: deltaT = delta_fmagma * latent_heat / cp
    ! This heat is distributed to the 4 corners evenly
    !
    !$OMP do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        do j = 1,nz-1
            cp_eff = Eff_cp( j,i )
            tmpr = 0.25d0*(temp0(j,i)+temp0(j+1,i)+temp0(j,i+1)+temp0(j+1,i+1))
            fr_lambda = lambda_freeze * exp(-lambda_freeze_tdep * (tmpr-t_top))
            delta_fmagma = min(fmagma(j,i), fmagma(j,i) * dt * fr_lambda)
            delta_fmagma2 = min(fmagma2(j,i), fmagma2(j,i) * dt * fr_lambda)
            fmagma(j,i) = fmagma(j,i) - delta_fmagma
            fmagma2(j,i) = fmagma2(j,i) - delta_fmagma2

            ! latent heat released by freezing magma
            deltaT1 = delta_fmagma * latent_heat_magma / cp_eff / 4
            deltaT2 = delta_fmagma2 * latent_heat_magma / cp_eff / 4
            deltaT = deltaT1 + deltaT2
            !$OMP atomic update
            !$ACC atomic update
            temp(j  ,i  ) = temp(j  ,i  ) + deltaT
            !$OMP atomic update
            !$ACC atomic update
            temp(j  ,i+1) = temp(j  ,i+1) + deltaT
            !$OMP atomic update
            !$ACC atomic update
            temp(j+1,i  ) = temp(j+1,i  ) + deltaT
            !$OMP atomic update
            !$ACC atomic update
            temp(j+1,i+1) = temp(j+1,i+1) + deltaT
        enddo
    enddo
    !$OMP end do
         
    temp = min(temp, t_bot)
    !$OMP end parallel     

    x_sum = 0.d0
    z_sum = 0.d0
    w_sum = 0.d0
    x_wide_l = 0.d0
    x_wide_r = 0.d0
    w_l = 0.d0
    w_r = 0.d0
    ! Whether i_center, j_center, h2, tan_mzone2_l/r, ihalfwidth_mzone2_l/r
    ! were actually computed this call. Reset every call (these are plain
    ! locals, not SAVEd) so a stale value from a previous timestep can
    ! never leak into the MOR-basalt migration block below.
    mor_zone_valid = .false.

    !$OMP Parallel private(i,j,weight,x)
    !Calculate melting center and magma width
    !center = sum( (x,z)*fmelt) / sum( fmelt )
    !width = 2 * ( sum( melt*(x-x_center)**2 ) / sum( melt ) )**0.5
    !_________________________surface
    !          /|\
    !         / | \
    !        /  |  \
    !       /   |   \
    !      /    |    \       * melting center
    !     /     |h    \
    !    /      |      \
    !   /       |       \
    !  /        |        \
    ! /         |         \
    ! ----------*----------
    !   wide_l     wide_r

    !$OMP do reduction(+:x_sum,z_sum,w_sum)
    do i = 1,nx-1
        do j = 1,nz-1
            weight = fmelt2(j,i)
            if (weight > w_barrier) then
                x_sum = x_sum + 0.25d0 * (cord(j,i,1)+cord(j,i+1,1)+cord(j+1,i,1)+cord(j+1,i+1,1)) * weight 
                z_sum = z_sum + 0.25d0 * (cord(j,i,2)+cord(j,i+1,2)+cord(j+1,i,2)+cord(j+1,i+1,2)) * weight
                w_sum = w_sum + weight
            endif
        enddo
    enddo
    !$OMP end do
    !$OMP single
    if (w_sum > 0.0d0) then
        x_center = x_sum / w_sum
        z_center = z_sum / w_sum
        found = .false.
        do i = 1,nx-1
            do j = 1,nz-1
                if (cord(j,i,1) <= x_center .and. cord(j,i+1,1) >= x_center) then
                    if (cord(j,i,2) >= z_center .and. cord(j+1,i,2) <= z_center) then
                        i_center = i
                        j_center = j
                        found = .true.
                        exit
                    endif
                endif
            enddo
            if (found) exit
        enddo
        ! Only now are i_center/j_center guaranteed to be set. Guard the
        ! rest of the melting-zone calculation (and the migration block
        ! that uses it, below) on this rather than on w_sum > 0.0d0 alone:
        ! w_sum > 0 doesn't guarantee the point-in-cell search above
        ! actually found a containing element.
        mor_zone_valid = found
    endif
    !$OMP end single
    !$OMP barrier
    if (mor_zone_valid) then
        !$OMP do reduction(+:x_wide_l,x_wide_r,w_l,w_r)
        do i = 1,nx-1
            do j = 1,nz-1
                weight = fmelt2(j,i)
                x = 0.25d0 * (cord(j,i,1)+cord(j,i+1,1)+cord(j+1,i,1)+cord(j+1,i+1,1))  
                if (weight > w_barrier .and. x < x_center) then
                    x_wide_l = x_wide_l + (x_center - x)**2 * weight
                    w_l = w_l + weight
                elseif (weight > w_barrier .and. x > x_center) then
                    x_wide_r = x_wide_r + (x - x_center)**2 * weight
                    w_r = w_r + weight
                endif
            enddo
        enddo
        !$OMP end do
        !$OMP single
        x_wide_l = 1.8d0*sqrt(x_wide_l/w_l)
        x_wide_r = 1.8d0*sqrt(x_wide_r/w_r)
        if (x_wide_l < dxmin) x_wide_l = dxmin
        if (x_wide_r < dxmin) x_wide_r = dxmin
        if (x_wide_l /= x_wide_l) x_wide_l = dxmin
        if (x_wide_r /= x_wide_r) x_wide_r = dxmin
        !write(333,*)'kk',i_l,i_r,x_wide_l,x_wide_r,time/sec_year/1d6
        z_surf = 0.5d0 * (cord(1,i_center,2) + cord(1,i_center+1,2))
        h2 = z_surf - z_center
        tan_mzone2_l = x_wide_l/h2
        tan_mzone2_r = x_wide_r/h2
        ihalfwidth_mzone2_l = ceiling(tan_mzone2_l * 150e3 / dxmin)
        ihalfwidth_mzone2_r = ceiling(tan_mzone2_r * 150e3 / dxmin)
        !write(333,*) 'oo',i_center,j_center,x_center,z_center
        !$OMP end single
    endif
    !$OMP end parallel

    ! Zero the this-timestep dike-migration accumulator (arc pathway); it's
    ! consumed by the post-pass below into dike_backlog/dike_released and
    ! must not carry over any value from a previous call.
    !$ACC kernels async(1)
    dike_added = 0.d0
    !$ACC end kernels

    !$OMP Parallel private(i,j,jm,quad_area,area_ratio,ii,jj,jtop,z_moho,z_melt,x_melt,h,x,z, &
    !$OMP                  ihalf,ii_start,ii_end,ncol)
    !$OMP do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        ! XXX: Assume melting cannot happen above the moho. (j > jmoho(i)) is always true
        ! starting from j=1, not jmoho(i), so that the loop can be collapsed and faster in computation
        do j = 1,nz-1
            jm = jmoho(i)
            quad_area = 0.5d0/area(j,i,1) + 0.5d0/area(j,i,2) ! area of this element
            if (j>jm .and. fmelt(j,i) > 0) then
                ! This element is under melting. The melt will migrate to the
                ! mantle and crust above. Treat the migration as instantaneous.

                ! Within crust, melts migrate by diking, propagate upward vertically over nelem_dike width
                ! area_ratio: the area of the crust column / the area of the melting element
                !     ~ the thickness of the crust column / the thickness of the melting element
                area_ratio = (cord(1,i,2)+cord(1,i+1,2)-cord(jm,i,2)-cord(jm,i+1,2)) / &
                    (cord(j,i,2)+cord(j,i+1,2)-cord(j+1,i,2)-cord(j+1,i+1,2))

                ihalf = max(0, (nelem_dike - 1) / 2)
                ii_start = max(1, i - ihalf)
                ii_end = min(nx-1, ii_start + max(1, nelem_dike) - 1)
                ii_start = max(1, ii_end - max(1, nelem_dike) + 1)
                ncol = ii_end - ii_start + 1

                do ii = ii_start, ii_end
                    ! Dikes stop at the topmost non-sediment element: don't
                    ! intrude into loose/compacting sediment (ksed1/ksed2)
                    ! sitting at the very top of the column.
                    jtop = 1
                    do while (jtop < jmoho(ii) .and. &
                              (iphase(jtop,ii) == ksed1 .or. iphase(jtop,ii) == ksed2))
                        jtop = jtop + 1
                    enddo
                    do jj = jtop, jmoho(ii)
                        !$OMP atomic update
                        !$ACC atomic update
                        fmagma(jj,ii) = fmagma(jj,ii) + fmelt(j,i) * ratio_crust_mzone * area_ratio * prod_magma * dt / ncol
                        !$OMP atomic update
                        !$ACC atomic update
                        dike_added(jj,ii) = dike_added(jj,ii) + fmelt(j,i) * ratio_crust_mzone * area_ratio * prod_magma * dt / ncol
                    enddo
                enddo

                ! Within mantle, melts migrate by percolation, propagate upward slantly
                z_moho = 0.5d0 * (cord(jm,i,2) + cord(jm,i+1,2))
                z_melt = 0.5d0 * (cord(j+1,i,2) + cord(j+1,i+1,2))
                x_melt = 0.5d0 * (cord(j+1,i,1) + cord(j+1,i+1,1))
                h = z_moho - z_melt
                ! area_ratio: the area of the mantle triangle / the area of the melting element
                area_ratio = h * h * tan_mzone / quad_area
                ! ii: the potential region of magma distribution zone at moho
                !write(333,*) x_center, z_center, i_center, j_center
                do ii = max(1,i-ihalfwidth_mzone), min(nx-1,i+ihalfwidth_mzone)
                    do jj = jmoho(ii)+1, j
                        x = 0.5d0 * (cord(jj,ii,1) + cord(jj,ii+1,1))
                        z = 0.5d0 * (cord(jj,ii,2) + cord(jj,ii+1,2))
                        if (abs(x - x_melt) <= tan_mzone * (z - z_melt)) then
                            !$OMP atomic update
                            !$ACC atomic update
                            fmagma(jj,ii) = fmagma(jj,ii) + fmelt(j,i) * ratio_mantle_mzone * area_ratio * prod_magma * dt
                        endif
                    enddo
                enddo
            endif

            ! MOR basalt melting. Also requires mor_zone_valid: i_center,
            ! j_center, h2, tan_mzone2_l/r, ihalfwidth_mzone2_l/r are only
            ! meaningful when the melting-center search above actually
            ! located a containing element this timestep.
            if (fmelt2(j,i) > 0 .and. mor_zone_valid) then
                ! This element is under melting.
                ! Within mantle, melts migrate by percolation, propagate upward slantly
                area_ratio = 0.5d0 * (h2 * h2 * tan_mzone2_l / quad_area + h2 * h2 * tan_mzone2_r / quad_area)
                ! ii: the potential region of magma distribution zone at moho
                do ii = max(1,i_center-ihalfwidth_mzone2_l), min(nx-1,i_center+ihalfwidth_mzone2_r)
                    do jj = 1, j_center
                        !write(333,*)ii,jj
                        x = 0.5d0 * (cord(jj,ii,1) + cord(jj,ii+1,1))
                        z = 0.5d0 * (cord(jj,ii,2) + cord(jj,ii+1,2))
                        !write(333,*)i_center,ihalfwidth_mzone2,j_center,x,x_center,tan_mzone2,z_surf,z
                        if (x < x_center .and. abs(x_center - x) <= tan_mzone2_l * (z_surf - z)) then
                            !$OMP atomic update
                            !$ACC atomic update
                            fmagma2(jj,ii)=fmagma2(jj,ii) + fmelt2(j,i) * area_ratio * prod_magma2 * dt
                            !$OMP atomic update
                            !$ACC atomic update
                            mor_marker_vol(jj,ii) = mor_marker_vol(jj,ii) + fmelt2(j,i) * area_ratio * prod_magma2 * dt
                        endif
                        if (x >= x_center .and. abs(x - x_center) <= tan_mzone2_r * (z_surf - z)) then
                            !$OMP atomic update
                            !$ACC atomic update
                            fmagma2(jj,ii)=fmagma2(jj,ii) + fmelt2(j,i) * area_ratio * prod_magma2 * dt
                            !$OMP atomic update
                            !$ACC atomic update
                            mor_marker_vol(jj,ii) = mor_marker_vol(jj,ii) + fmelt2(j,i) * area_ratio * prod_magma2 * dt
                        endif
                    enddo
                enddo
            endif
            fmagma(j,i) = min(fmagma(j,i), fmagma_max)
            fmagma2(j,i) = min(fmagma2(j,i), fmagma_max)
        enddo
    enddo
    !$OMP end do
    !$OMP end parallel

    ! Rate-cap this timestep's dike-migration accumulation (arc pathway) at
    ! 0.01% of Amin (the initial grid's minimum element area) per timestep,
    ! banking whatever exceeds that in dike_backlog for release later. The
    ! released amount feeds fl_rheol.f90's eigenstrain term (dv1p) and also
    ! accumulates toward the next new marker. The MOR pathway (fmagma2) has
    ! no rate cap; its marker accumulator just tracks fmagma2's own
    ! additions directly.
    !$OMP parallel do private(ii,jj,elem_area,pending,released,kinc,marker_vol_unit,n_to_add,kk,xc,yc,inc) collapse(2)
    !$ACC parallel loop collapse(2) async(1)
    do ii = 1, nx-1
        do jj = 1, nz-1
            elem_area = 0.5d0/area(jj,ii,1) + 0.5d0/area(jj,ii,2)

            pending = dike_added(jj,ii)*elem_area + dike_backlog(jj,ii)
            released = min(pending, 1.0d-4*Amin)
            dike_backlog(jj,ii) = pending - released
            dike_released(jj,ii) = released / elem_area
            dike_marker_vol(jj,ii) = dike_marker_vol(jj,ii) + released

            ! New marker, arc pathway: once enough released dike volume has
            ! accumulated to represent at least one marker's average share
            ! of this element, add fresh karc1 markers at the centroid.
            kinc = nmark_elem(jj,ii)
            if (kinc > 0 .and. kinc < max_markers_per_elem) then
                marker_vol_unit = elem_area / kinc
                if (dike_marker_vol(jj,ii) >= marker_vol_unit) then
                    n_to_add = min(ceiling(dike_marker_vol(jj,ii)/marker_vol_unit), max_markers_per_elem - kinc)
                    xc = 0.25d0*(cord(jj,ii,1)+cord(jj+1,ii,1)+cord(jj,ii+1,1)+cord(jj+1,ii+1,1))
                    yc = 0.25d0*(cord(jj,ii,2)+cord(jj+1,ii,2)+cord(jj,ii+1,2)+cord(jj+1,ii+1,2))
                    do kk = 1, n_to_add
                        call add_marker(xc, yc, karc1, time, jj, ii, inc)
                    enddo
                    dike_marker_vol(jj,ii) = 0.d0
                    call count_phase_ratio(jj,ii)
                endif
            endif

            ! New marker, MOR pathway: same idea, using mor_marker_vol and
            ! kocean1, no rate cap (mirrors fl_move.f90's existing fmagma2
            ! extrusion-conversion mechanism, but for the migrated-into-crust
            ! portion rather than what erupts at the surface).
            kinc = nmark_elem(jj,ii)
            if (kinc > 0 .and. kinc < max_markers_per_elem) then
                marker_vol_unit = elem_area / kinc
                if (mor_marker_vol(jj,ii) >= marker_vol_unit) then
                    n_to_add = min(ceiling(mor_marker_vol(jj,ii)/marker_vol_unit), max_markers_per_elem - kinc)
                    xc = 0.25d0*(cord(jj,ii,1)+cord(jj+1,ii,1)+cord(jj,ii+1,1)+cord(jj+1,ii+1,1))
                    yc = 0.25d0*(cord(jj,ii,2)+cord(jj+1,ii,2)+cord(jj,ii+1,2)+cord(jj+1,ii+1,2))
                    do kk = 1, n_to_add
                        call add_marker(xc, yc, kocean1, time, jj, ii, inc)
                    enddo
                    mor_marker_vol(jj,ii) = 0.d0
                    call count_phase_ratio(jj,ii)
                endif
            endif
        enddo
    enddo
endif

!$ACC parallel loop collapse(2) async(1)
!$OMP Parallel private(i,j,iph,cp_eff,cond_eff,dissip,inv_area,real_area13,area_n,rhs, &
!$OMP                  x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  shpdx_1,shpdx_2,shpdx_3,shpdz_1,shpdz_2,shpdz_3)
!$OMP do
do i = 1,nx-1
    do j = 1,nz-1

        iph = iphase(j,i)

        ! Calculating effective material properties
        cp_eff = Eff_cp( j,i )
        cond_eff = Eff_conduct( j,i )

        ! if shearh-heating flag is true
        if( ishearh.eq.1 .and. itherm.ne.2 ) then
            dissip = shrheat(j,i)
        else
            dissip = 0
        endif

        ! diffusivity
        diff_elem(j,i) = cond_eff/den(iph)/cp_eff

        ! Additional sources - radiogenic and shear heating
        dummye(j,i) = ( source(j,i) + dissip/den(iph) ) / cp_eff

    end do
end do    
!$OMP end do

!$ACC parallel loop collapse(2) async(1)
!$OMP do
do i = 1,nx
    do j = 1,nz

        rhs = 0.d0
        area_n = 0.d0

        ! Element (j-1,i-1). Triangle B
        if ( j.ne.1 .and. i.ne.1 ) then
            inv_area = 1.d0/area(j-1,i-1,2)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            y2 = cord(j  ,i-1,2)
            y3 = cord(j-1,i  ,2)
            y4 = cord(j  ,i  ,2)
            x2 = cord(j  ,i-1,1)
            x3 = cord(j-1,i  ,1)
            x4 = cord(j  ,i  ,1)
            shpdx_1 = (y2 - y4) * area(j-1,i-1,2)
            shpdx_2 = (y4 - y3) * area(j-1,i-1,2)
            shpdx_3 = (y3 - y2) * area(j-1,i-1,2)
            shpdz_1 = (x4 - x2) * area(j-1,i-1,2)
            shpdz_2 = (x3 - x4) * area(j-1,i-1,2)
            shpdz_3 = (x2 - x3) * area(j-1,i-1,2)
            rhs = rhs + dummye(j-1,i-1) * real_area13 - 0.5d0 * diff_elem(j-1,i-1) * inv_area * ( &
                temp(j-1,i  ) * (shpdx_1 * shpdx_3 + shpdz_1 * shpdz_3) + &
                temp(j  ,i-1) * (shpdx_2 * shpdx_3 + shpdz_2 * shpdz_3) + &
                temp(j  ,i  ) * (shpdx_3 * shpdx_3 + shpdz_3 * shpdz_3) )
        endif

        ! Element (j-1,i). Triangles A,B
        if ( j.ne.1 .and. i.ne.nx ) then
            ! triangle A
            inv_area = 1.d0/area(j-1,i,1)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            y1 = cord(j-1,i  ,2)
            y2 = cord(j  ,i  ,2)
            y3 = cord(j-1,i+1,2)
            x1 = cord(j-1,i  ,1)
            x2 = cord(j  ,i  ,1)
            x3 = cord(j-1,i+1,1)
            shpdx_1 = (y2 - y3) * area(j-1,i,1)
            shpdx_2 = (y3 - y1) * area(j-1,i,1)
            shpdx_3 = (y1 - y2) * area(j-1,i,1)
            shpdz_1 = (x3 - x2) * area(j-1,i,1)
            shpdz_2 = (x1 - x3) * area(j-1,i,1)
            shpdz_3 = (x2 - x1) * area(j-1,i,1)
            rhs = rhs + dummye(j-1,i) * real_area13 - 0.5d0 * diff_elem(j-1,i) * inv_area * ( &
                temp(j-1,i  ) * (shpdx_1 * shpdx_2 + shpdz_1 * shpdz_2) + &
                temp(j  ,i  ) * (shpdx_2 * shpdx_2 + shpdz_2 * shpdz_2) + &
                temp(j-1,i+1) * (shpdx_3 * shpdx_2 + shpdz_3 * shpdz_2) )

            ! triangle B
            inv_area = 1.d0/area(j-1,i,2)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            y2 = cord(j  ,i  ,2)
            y3 = cord(j-1,i+1,2)
            y4 = cord(j  ,i+1,2)
            x2 = cord(j  ,i  ,1)
            x3 = cord(j-1,i+1,1)
            x4 = cord(j  ,i+1,1)
            shpdx_1 = (y2 - y4) * area(j-1,i,2)
            shpdx_2 = (y4 - y3) * area(j-1,i,2)
            shpdx_3 = (y3 - y2) * area(j-1,i,2)
            shpdz_1 = (x4 - x2) * area(j-1,i,2)
            shpdz_2 = (x3 - x4) * area(j-1,i,2)
            shpdz_3 = (x2 - x3) * area(j-1,i,2)
            rhs = rhs + dummye(j-1,i) * real_area13 - 0.5d0 * diff_elem(j-1,i) * inv_area * ( &
                temp(j-1,i+1) * (shpdx_1 * shpdx_2 + shpdz_1 * shpdz_2) + &
                temp(j  ,i  ) * (shpdx_2 * shpdx_2 + shpdz_2 * shpdz_2) + &
                temp(j  ,i+1) * (shpdx_3 * shpdx_2 + shpdz_3 * shpdz_2) )
        endif
        
        ! Element (j,i-1). Triangles A,B
        if ( j.ne.nz .and. i.ne.1 ) then
            ! triangle A
            inv_area = 1.d0/area(j,i-1,1)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            y1 = cord(j  ,i-1,2)
            y2 = cord(j+1,i-1,2)
            y3 = cord(j  ,i  ,2)
            x1 = cord(j  ,i-1,1)
            x2 = cord(j+1,i-1,1)
            x3 = cord(j  ,i  ,1)
            shpdx_1 = (y2 - y3) * area(j,i-1,1)
            shpdx_2 = (y3 - y1) * area(j,i-1,1)
            shpdx_3 = (y1 - y2) * area(j,i-1,1)
            shpdz_1 = (x3 - x2) * area(j,i-1,1)
            shpdz_2 = (x1 - x3) * area(j,i-1,1)
            shpdz_3 = (x2 - x1) * area(j,i-1,1)
            rhs = rhs + dummye(j,i-1) * real_area13 - 0.5d0 * diff_elem(j,i-1) * inv_area * ( &
                temp(j  ,i-1) * (shpdx_1 * shpdx_3 + shpdz_1 * shpdz_3) + &
                temp(j+1,i-1) * (shpdx_2 * shpdx_3 + shpdz_2 * shpdz_3) + &
                temp(j  ,i  ) * (shpdx_3 * shpdx_3 + shpdz_3 * shpdz_3) )

            ! triangle B
            inv_area = 1.d0/area(j,i-1,2)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            y2 = cord(j+1,i-1,2)
            y3 = cord(j  ,i  ,2)
            y4 = cord(j+1,i  ,2)
            x2 = cord(j+1,i-1,1)
            x3 = cord(j  ,i  ,1)
            x4 = cord(j+1,i  ,1)
            shpdx_1 = (y2 - y4) * area(j,i-1,2)
            shpdx_2 = (y4 - y3) * area(j,i-1,2)
            shpdx_3 = (y3 - y2) * area(j,i-1,2)
            shpdz_1 = (x4 - x2) * area(j,i-1,2)
            shpdz_2 = (x3 - x4) * area(j,i-1,2)
            shpdz_3 = (x2 - x3) * area(j,i-1,2)
            rhs = rhs + dummye(j,i-1) * real_area13 - 0.5d0 * diff_elem(j,i-1) * inv_area * ( &
                temp(j  ,i  ) * (shpdx_1 * shpdx_1 + shpdz_1 * shpdz_1) + &
                temp(j+1,i-1) * (shpdx_2 * shpdx_1 + shpdz_2 * shpdz_1) + &
                temp(j+1,i  ) * (shpdx_3 * shpdx_1 + shpdz_3 * shpdz_1) )
        endif

        ! Element (j,i). Triangle A
        if ( j.ne.nz .and. i.ne.nx ) then
            inv_area = 1.d0/area(j,i,1)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            y1 = cord(j  ,i  ,2)
            y2 = cord(j+1,i  ,2)
            y3 = cord(j  ,i+1,2)
            x1 = cord(j  ,i  ,1)
            x2 = cord(j+1,i  ,1)
            x3 = cord(j  ,i+1,1)
            shpdx_1 = (y2 - y3) * area(j,i,1)
            shpdx_2 = (y3 - y1) * area(j,i,1)
            shpdx_3 = (y1 - y2) * area(j,i,1)
            shpdz_1 = (x3 - x2) * area(j,i,1)
            shpdz_2 = (x1 - x3) * area(j,i,1)
            shpdz_3 = (x2 - x1) * area(j,i,1)
            rhs = rhs + dummye(j,i) * real_area13 - 0.5d0 * diff_elem(j,i) * inv_area * ( &
                temp(j  ,i  ) * (shpdx_1 * shpdx_1 + shpdz_1 * shpdz_1) + &
                temp(j+1,i  ) * (shpdx_2 * shpdx_1 + shpdz_2 * shpdz_1) + &
                temp(j  ,i+1) * (shpdx_3 * shpdx_1 + shpdz_3 * shpdz_1) )
        endif

        ! Update Temperature by Eulerian method 
        temp(j,i) = temp(j,i)+rhs*dt/area_n
    end do
end do
!$OMP end do

! Boundary conditions (top and bottom)
!$ACC parallel loop async(1)
!$OMP do
do i = 1,nx

    temp(1,i) = t_top

    if( itemp_bc.eq.1 ) then
        temp(nz,i) = bot_bc
    elseif( itemp_bc.eq.2 ) then
        cond_eff = Eff_conduct( nz-1, min(i,nx-1) )
        temp(nz,i) = temp(nz-1,i)  +  bot_bc * ( cord(nz-1,i,2)-cord(nz,i,2) ) / cond_eff
    endif

end do
!$OMP end do

! Boundary conditions: dt/dx =0 on left and right  
!$ACC parallel loop async(1)
!$OMP do
do j = 1,nz
    temp(j ,1)  = temp(j,2)
    temp(j, nx) = temp(j,nx-1)
end do
!$OMP end do
!$OMP end parallel

return
end 
