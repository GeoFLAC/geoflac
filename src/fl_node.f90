
!  Calculations of forces from stresses
subroutine fl_node
use arrays
use params
implicit none
integer, parameter :: kk(3, 4) = reshape((/ 2,3,4, 1,2,3, 1,2,4, 1,3,4 /), (/3, 4/))
integer, parameter :: mm(3, 4) = reshape((/ 3,3,2, 2,2,2, 3,1,3, 1,1,1 /), (/3, 4/))

integer :: i, j, k, je, ie, iunknown
double precision :: x1, y1, x2, y2, x3, y3, x4, y4, drat, fx, fy, factor, p_est, rosubg, press_norm_l, dlx_l, dly_l, press_norm_r, dlx_r, dly_r, rho_water_g, water_depth
double precision :: shpdx_loc(3, 4), shpdz_loc(3, 4)

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
!
!  diagonal \ :
!
!   C:        D:
!
!  1          1---3
!  | \         \  |
!  2---3          2
!
!    assemblage of forces is COUNTRE CLOCK-WISE !
!

drat = dt / dt_elastic
if (drat .lt. 1.d0) drat = 1.d0

!$OMP parallel private(i, j, k, fx, fy, factor, &
!$OMP                  p_est, rosubg, &
!$OMP                  press_norm_l, dlx_l, dly_l, &
!$OMP                  press_norm_r, dlx_r, dly_r, &
!$OMP                  iunknown, rho_water_g, water_depth, &
!$OMP                  je, ie, x1, y1, x2, y2, x3, y3, x4, y4, &
!$OMP                  shpdx_loc, shpdz_loc)
!
!$OMP do
!$ACC parallel loop collapse(2) async(1)
do i = 1,nx
    do j = 1,nz
        if(nystressbc.eq.0) then
           force(j,i,1) = 0
           force(j,i,2) = 0
        endif
        ! REGULAR PART - forces from stresses
        
        ! Element (j-1,i-1). Triangles B,C,D
        if ( j.ne.1 .and. i.ne.1 ) then
            je = j-1
            ie = i-1
            x1 = cord(je  ,ie  ,1)
            y1 = cord(je  ,ie  ,2)
            x2 = cord(je+1,ie  ,1)
            y2 = cord(je+1,ie  ,2)
            x3 = cord(je  ,ie+1,1)
            y3 = cord(je  ,ie+1,2)
            x4 = cord(je+1,ie+1,1)
            y4 = cord(je+1,ie+1,2)

            ! Triangle B (k=2)
            shpdx_loc(1, 2) = (y2 - y4) * area(je, ie, 2)
            shpdx_loc(2, 2) = (y4 - y3) * area(je, ie, 2)
            shpdx_loc(3, 2) = (y3 - y2) * area(je, ie, 2)
            shpdz_loc(1, 2) = (x4 - x2) * area(je, ie, 2)
            shpdz_loc(2, 2) = (x3 - x4) * area(je, ie, 2)
            shpdz_loc(3, 2) = (x2 - x3) * area(je, ie, 2)

            ! Triangle C (k=3)
            shpdx_loc(1, 3) = (y2 - y4) * area(je, ie, 3)
            shpdx_loc(2, 3) = (y4 - y1) * area(je, ie, 3)
            shpdx_loc(3, 3) = (y1 - y2) * area(je, ie, 3)
            shpdz_loc(1, 3) = (x4 - x2) * area(je, ie, 3)
            shpdz_loc(2, 3) = (x1 - x4) * area(je, ie, 3)
            shpdz_loc(3, 3) = (x2 - x1) * area(je, ie, 3)

            ! Triangle D (k=4)
            shpdx_loc(1, 4) = (y4 - y3) * area(je, ie, 4)
            shpdx_loc(2, 4) = (y3 - y1) * area(je, ie, 4)
            shpdx_loc(3, 4) = (y1 - y4) * area(je, ie, 4)
            shpdz_loc(1, 4) = (x3 - x4) * area(je, ie, 4)
            shpdz_loc(2, 4) = (x1 - x3) * area(je, ie, 4)
            shpdz_loc(3, 4) = (x4 - x1) * area(je, ie, 4)

            do k = 1, 3
                factor = 0.25d0 / area(je,ie,kk(k,1))
                fx = factor * (stress0(je,ie,1,kk(k,1)) * shpdx_loc(mm(k,1),kk(k,1)) + &
                               stress0(je,ie,3,kk(k,1)) * shpdz_loc(mm(k,1),kk(k,1)))
                fy = factor * (stress0(je,ie,3,kk(k,1)) * shpdx_loc(mm(k,1),kk(k,1)) + &
                               stress0(je,ie,2,kk(k,1)) * shpdz_loc(mm(k,1),kk(k,1)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! Element (j-1,i). Triangles A,B,C.
        if ( j.ne.1 .and. i.ne.nx ) then
            je = j-1
            ie = i
            x1 = cord(je  ,ie  ,1)
            y1 = cord(je  ,ie  ,2)
            x2 = cord(je+1,ie  ,1)
            y2 = cord(je+1,ie  ,2)
            x3 = cord(je  ,ie+1,1)
            y3 = cord(je  ,ie+1,2)
            x4 = cord(je+1,ie+1,1)
            y4 = cord(je+1,ie+1,2)

            ! Triangle A (k=1)
            shpdx_loc(1, 1) = (y2 - y3) * area(je, ie, 1)
            shpdx_loc(2, 1) = (y3 - y1) * area(je, ie, 1)
            shpdx_loc(3, 1) = (y1 - y2) * area(je, ie, 1)
            shpdz_loc(1, 1) = (x3 - x2) * area(je, ie, 1)
            shpdz_loc(2, 1) = (x1 - x3) * area(je, ie, 1)
            shpdz_loc(3, 1) = (x2 - x1) * area(je, ie, 1)

            ! Triangle B (k=2)
            shpdx_loc(1, 2) = (y2 - y4) * area(je, ie, 2)
            shpdx_loc(2, 2) = (y4 - y3) * area(je, ie, 2)
            shpdx_loc(3, 2) = (y3 - y2) * area(je, ie, 2)
            shpdz_loc(1, 2) = (x4 - x2) * area(je, ie, 2)
            shpdz_loc(2, 2) = (x3 - x4) * area(je, ie, 2)
            shpdz_loc(3, 2) = (x2 - x3) * area(je, ie, 2)

            ! Triangle C (k=3)
            shpdx_loc(1, 3) = (y2 - y4) * area(je, ie, 3)
            shpdx_loc(2, 3) = (y4 - y1) * area(je, ie, 3)
            shpdx_loc(3, 3) = (y1 - y2) * area(je, ie, 3)
            shpdz_loc(1, 3) = (x4 - x2) * area(je, ie, 3)
            shpdz_loc(2, 3) = (x1 - x4) * area(je, ie, 3)
            shpdz_loc(3, 3) = (x2 - x1) * area(je, ie, 3)

            do k = 1, 3
                factor = 0.25d0 / area(je,ie,kk(k,2))
                fx = factor * (stress0(je,ie,1,kk(k,2)) * shpdx_loc(mm(k,2),kk(k,2)) + &
                               stress0(je,ie,3,kk(k,2)) * shpdz_loc(mm(k,2),kk(k,2)))
                fy = factor * (stress0(je,ie,3,kk(k,2)) * shpdx_loc(mm(k,2),kk(k,2)) + &
                               stress0(je,ie,2,kk(k,2)) * shpdz_loc(mm(k,2),kk(k,2)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! Element (j,i-1). Triangles A,B,D
        if ( j.ne.nz .and. i.ne.1 ) then
            je = j
            ie = i-1
            x1 = cord(je  ,ie  ,1)
            y1 = cord(je  ,ie  ,2)
            x2 = cord(je+1,ie  ,1)
            y2 = cord(je+1,ie  ,2)
            x3 = cord(je  ,ie+1,1)
            y3 = cord(je  ,ie+1,2)
            x4 = cord(je+1,ie+1,1)
            y4 = cord(je+1,ie+1,2)

            ! Triangle A (k=1)
            shpdx_loc(1, 1) = (y2 - y3) * area(je, ie, 1)
            shpdx_loc(2, 1) = (y3 - y1) * area(je, ie, 1)
            shpdx_loc(3, 1) = (y1 - y2) * area(je, ie, 1)
            shpdz_loc(1, 1) = (x3 - x2) * area(je, ie, 1)
            shpdz_loc(2, 1) = (x1 - x3) * area(je, ie, 1)
            shpdz_loc(3, 1) = (x2 - x1) * area(je, ie, 1)

            ! Triangle B (k=2)
            shpdx_loc(1, 2) = (y2 - y4) * area(je, ie, 2)
            shpdx_loc(2, 2) = (y4 - y3) * area(je, ie, 2)
            shpdx_loc(3, 2) = (y3 - y2) * area(je, ie, 2)
            shpdz_loc(1, 2) = (x4 - x2) * area(je, ie, 2)
            shpdz_loc(2, 2) = (x3 - x4) * area(je, ie, 2)
            shpdz_loc(3, 2) = (x2 - x3) * area(je, ie, 2)

            ! Triangle D (k=4)
            shpdx_loc(1, 4) = (y4 - y3) * area(je, ie, 4)
            shpdx_loc(2, 4) = (y3 - y1) * area(je, ie, 4)
            shpdx_loc(3, 4) = (y1 - y4) * area(je, ie, 4)
            shpdz_loc(1, 4) = (x3 - x4) * area(je, ie, 4)
            shpdz_loc(2, 4) = (x1 - x3) * area(je, ie, 4)
            shpdz_loc(3, 4) = (x4 - x1) * area(je, ie, 4)

            do k = 1, 3
                factor = 0.25d0 / area(je,ie,kk(k,3))
                fx = factor * (stress0(je,ie,1,kk(k,3)) * shpdx_loc(mm(k,3),kk(k,3)) + &
                               stress0(je,ie,3,kk(k,3)) * shpdz_loc(mm(k,3),kk(k,3)))
                fy = factor * (stress0(je,ie,3,kk(k,3)) * shpdx_loc(mm(k,3),kk(k,3)) + &
                               stress0(je,ie,2,kk(k,3)) * shpdz_loc(mm(k,3),kk(k,3)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! Element (j,i). Triangles A,C,D
        if ( j.ne.nz .and. i.ne.nx ) then
            je = j
            ie = i
            x1 = cord(je  ,ie  ,1)
            y1 = cord(je  ,ie  ,2)
            x2 = cord(je+1,ie  ,1)
            y2 = cord(je+1,ie  ,2)
            x3 = cord(je  ,ie+1,1)
            y3 = cord(je  ,ie+1,2)
            x4 = cord(je+1,ie+1,1)
            y4 = cord(je+1,ie+1,2)

            ! Triangle A (k=1)
            shpdx_loc(1, 1) = (y2 - y3) * area(je, ie, 1)
            shpdx_loc(2, 1) = (y3 - y1) * area(je, ie, 1)
            shpdx_loc(3, 1) = (y1 - y2) * area(je, ie, 1)
            shpdz_loc(1, 1) = (x3 - x2) * area(je, ie, 1)
            shpdz_loc(2, 1) = (x1 - x3) * area(je, ie, 1)
            shpdz_loc(3, 1) = (x2 - x1) * area(je, ie, 1)

            ! Triangle C (k=3)
            shpdx_loc(1, 3) = (y2 - y4) * area(je, ie, 3)
            shpdx_loc(2, 3) = (y4 - y1) * area(je, ie, 3)
            shpdx_loc(3, 3) = (y1 - y2) * area(je, ie, 3)
            shpdz_loc(1, 3) = (x4 - x2) * area(je, ie, 3)
            shpdz_loc(2, 3) = (x1 - x4) * area(je, ie, 3)
            shpdz_loc(3, 3) = (x2 - x1) * area(je, ie, 3)

            ! Triangle D (k=4)
            shpdx_loc(1, 4) = (y4 - y3) * area(je, ie, 4)
            shpdx_loc(2, 4) = (y3 - y1) * area(je, ie, 4)
            shpdx_loc(3, 4) = (y1 - y4) * area(je, ie, 4)
            shpdz_loc(1, 4) = (x3 - x4) * area(je, ie, 4)
            shpdz_loc(2, 4) = (x1 - x3) * area(je, ie, 4)
            shpdz_loc(3, 4) = (x4 - x1) * area(je, ie, 4)

            do k = 1, 3
                factor = 0.25d0 / area(je,ie,kk(k,4))
                fx = factor * (stress0(je,ie,1,kk(k,4)) * shpdx_loc(mm(k,4),kk(k,4)) + &
                               stress0(je,ie,3,kk(k,4)) * shpdz_loc(mm(k,4),kk(k,4)))
                fy = factor * (stress0(je,ie,3,kk(k,4)) * shpdx_loc(mm(k,4),kk(k,4)) + &
                               stress0(je,ie,2,kk(k,4)) * shpdz_loc(mm(k,4),kk(k,4)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! GRAVITY FORCE
        force(j,i,2) = force(j,i,2) - rmass(j,i)*g

  enddo
enddo
!$OMP end do

! BOUNDARY CONDITIONS
if(nyhydro.gt.0) then
    !$OMP do
    !$ACC parallel loop async(1)
    do i=1,nx

        ! pressure from water sea on top
        rho_water_g = 1030.d0 * g
        if(i.lt.nx) then
            water_depth = 0.5d0*(cord(1,i+1,2)+cord(1,i,2))
        else
            water_depth = 0.5d0*(cord(1,i-1,2)+cord(1,i,2))
        endif

        if (water_depth.lt.0.d0) then ! No water (above sea level)
            if(i.eq.1) then
                press_norm_l = 0
                dlx_l = 0
                dly_l = 0
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.d0)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            elseif(i.eq.nx) then
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.d0)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = 0
                dlx_r = 0
                dly_r = 0
            else
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.d0)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.d0)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            endif
            force(1,i,1) = force(1,i,1)-0.5d0*press_norm_l*dly_l-0.5d0*press_norm_r*dly_r
            force(1,i,2) = force(1,i,2)+0.5d0*press_norm_l*dlx_l+0.5d0*press_norm_r*dlx_r
        endif
    enddo
    !$OMP end do

    !$OMP do
    !$ACC parallel loop async(1)
    do i=1,nx

        ! bottom support - Archimed force (normal to the surface, shear component = 0)
        p_est = pisos + 0.5d0*(den(iphsub)+drosub)*g*(cord(nz,i,2)-rzbo)
        rosubg = g * (den(iphsub)+drosub) * (1-alfa(iphsub)*temp(nz,i)+beta(iphsub)*p_est)

        if(i.eq.1) then
            press_norm_l = 0
            dlx_l = 0
            dly_l = 0

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        elseif(i.eq.nx) then
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = 0
            dlx_r = 0
            dly_r = 0
        else
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        endif
            
        force(nz,i,1) = force(nz,i,1)-0.5d0*press_norm_l*dly_l-0.5d0*press_norm_r*dly_r
        force(nz,i,2) = force(nz,i,2)+0.5d0*press_norm_l*dlx_l+0.5d0*press_norm_r*dlx_r

        !write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r

    enddo
    !$OMP end do
endif

!$OMP do
!$ACC parallel loop collapse(2) async(1)
do i=1,nx
    do j=1,nz

        ! DAMPING
        if( iand(ncod(j,i,1),1).ne.1 .and. abs(vel(j,i,1)).gt.1.0d-13 ) then
            force(j,i,1) = force(j,i,1) - demf*sign(force(j,i,1),vel(j,i,1))
        endif

        if( iand(ncod(j,i,2),2).ne.2 .and. abs(vel(j,i,2)).gt.1.0d-13 ) then
            force(j,i,2) = force(j,i,2) - demf*sign(force(j,i,2),vel(j,i,2))
        endif

        ! VELOCITIES FROM FORCES
        if( ncod(j,i,1) .eq. 1 ) then
            vel(j,i,1) = bc(j,i,1) 
        else
            vel(j,i,1) = vel(j,i,1) + dt*force(j,i,1)/(amass(j,i)*drat*drat)
        endif
        if( ncod(j,i,2) .eq. 1 ) then
            vel(j,i,2) = bc(j,i,2)
        else
            vel(j,i,2) = vel(j,i,2) + dt*force(j,i,2)/(amass(j,i)*drat*drat)
        endif
    end do
end do
!$OMP end do
!$OMP end parallel
if (i_prestress.eq.1.and.time.lt.600.d3*sec_year) then
     !$OMP parallel do private(i,k)
     !$ACC parallel loop collapse(2) async(1)
     do k = 1,2
        do i = 1, nx
            vel(nz,i,k) = 0.d0
        enddo
     enddo
     !$OMP end parallel do

     !$OMP parallel do private(j,k)
     !$ACC parallel loop collapse(2) async(1)
     do k = 1,2
        do j = 1, nz
            vel(j,1,k) = 0.d0
            vel(j,nx,k) = 0.d0
        enddo
     enddo
     !$OMP end parallel do
endif
return
end subroutine fl_node
