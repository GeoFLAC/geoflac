
! Calculate thermal field in explict form

subroutine fl_therm
use arrays
use params
use matprops
implicit none

double precision :: flux(nz-1,nx-1,2,2), add_source(nz-1,nx-1)
double precision, parameter :: heat_latent_magma = 4.2d5  ! J/kg, latent heat of freezing magma

integer :: i,j,iph
double precision :: cp_eff,cond_eff,dissip,diff,quad_area, &
                  x1,x2,x3,x4,y1,y2,y3,y4,t1,t2,t3,t4,tmpr, &
                  qs,real_area13,area_n,rhs

! real_area = 0.5* (1./area(n,t))
! Calculate Fluxes in every triangle
!      flux (j,i,num_triangle, direction(x,y)
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


!!ntherm = ntherm+1
!!$ACC update device(ntherm)

! saving old temperature
!$DIR PREFER_PARALLEL
!$ACC parallel 
if (istress_therm.gt.0) temp0(1:nz,1:nx) = temp(1:nz,1:nx)
!$ACC end parallel

!$OMP Parallel private(i,j,iph,cp_eff,cond_eff,dissip,diff,quad_area, &
!$OMP                  x1,x2,x3,x4,y1,y2,y3,y4,t1,t2,t3,t4,tmpr, &
!$OMP                  qs,real_area13,area_n,rhs)
!$OMP do
!$ACC parallel loop private(cp_eff, quad_area)
do i = 1,nx-1
    j = 1  ! top
    !iph = iphase(j,i)
    cp_eff = Eff_cp( j,i )

    ! area(j,i) is INVERSE of "real" DOUBLE area (=1./det)
    quad_area = 1.d0/(area(j,i,1)+area(j,i,2))

    temp(j,i  ) = temp(j,i  ) + andesitic_melt_vol(i) * heat_latent_magma / quad_area / cp_eff
    temp(j,i+1) = temp(j,i+1) + andesitic_melt_vol(i) * heat_latent_magma / quad_area / cp_eff
end do
!$OMP end do
!$ACC end parallel

!$ACC data create(add_source,flux)
!$ACC parallel loop collapse(2) private(iph, cp_eff, cond_eff, dissip, diff, &
!$ACC                                   x1, x2, x3, x4, y1, y2, y3, y4, tmpr) 
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
        diff = cond_eff/den(iph)/cp_eff

        ! Calculate fluxes in two triangles
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)
        t1 = temp (j   ,i  )
        t2 = temp (j+1 ,i  )
        t3 = temp (j   ,i+1)
        t4 = temp (j+1 ,i+1)

        ! Additional sources - radiogenic and shear heating
        tmpr = 0.25d0*(t1 + t2 + t3 + t4)
        !add_source(j,i) = ( source(j,i) + dissip/den(iph) - 600.d0*cp_eff*Eff_melt(iph,tmpr)) / cp_eff
        add_source(j,i) = ( source(j,i) + dissip/den(iph) ) / cp_eff

        ! (1) A element:
        flux(j,i,1,1) = -diff * ( t1*(y2-y3)+t2*(y3-y1)+t3*(y1-y2) ) * area(j,i,1)
        flux(j,i,1,2) = -diff * ( t1*(x3-x2)+t2*(x1-x3)+t3*(x2-x1) ) * area(j,i,1)
 
        ! (2) B element: Interchange of numeration: (1 -> 3,  3 -> 4)
        flux(j,i,2,1) = -diff * ( t3*(y2-y4)+t2*(y4-y3)+t4*(y3-y2) ) * area(j,i,2)
        flux(j,i,2,2) = -diff * ( t3*(x4-x2)+t2*(x3-x4)+t4*(x2-x3) ) * area(j,i,2)

    end do
end do    
!$OMP end do
!$ACC end parallel

!$ACC parallel loop collapse(2) private(rhs, area_n, rhs,real_area13)
!$OMP do
do i = 1,nx
    do j = 1,nz

        rhs = 0
        area_n = 0

        ! Element (j-1,i-1). Triangle B
        if ( j.ne.1 .and. i.ne.1 ) then

            ! side 2-3
            qs = flux(j-1,i-1,2,1) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 flux(j-1,i-1,2,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(j-1,i-1,2,1) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 flux(j-1,i-1,2,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j-1,i-1,2)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + add_source(j-1,i-1)*real_area13

        endif

        ! Element (j-1,i). Triangles A,B
        if ( j.ne.1 .and. i.ne.nx ) then

            ! triangle A
            ! side 1-2
            qs = flux(j-1,i  ,1,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 flux(j-1,i  ,1,2) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            rhs = rhs + 0.5d0*qs

            ! side 2-3
            qs = flux(j-1,i  ,1,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 flux(j-1,i  ,1,2) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j-1,i  ,1)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + add_source(j-1,i  )*real_area13

            ! triangle B
            ! side 1-2
            qs = flux(j-1,i  ,2,1) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 flux(j-1,i  ,2,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            rhs = rhs + 0.5d0*qs

            ! side 2-3
            qs = flux(j-1,i  ,2,1) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 flux(j-1,i  ,2,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j-1,i  ,2)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + add_source(j-1,i  )*real_area13

        endif
        
        ! Element (j,i-1). Triangles A,B
        if ( j.ne.nz .and. i.ne.1 ) then

            ! triangle A
            ! side 2-3
            qs = flux(j  ,i-1,1,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 flux(j  ,i-1,1,2) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(j  ,i-1,1,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 flux(j  ,i-1,1,2) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j  ,i-1,1)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + add_source(j  ,i-1)*real_area13

            ! triangle B
            ! side 1-2
            qs = flux(j  ,i-1,2,1) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 flux(j  ,i-1,2,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(j  ,i-1,2,1) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 flux(j  ,i-1,2,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j  ,i-1,2)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + add_source(j  ,i-1)*real_area13

        endif

        ! Element (j,i). Triangle A
        if ( j.ne.nz .and. i.ne.nx ) then

            ! side 1-2
            qs = flux(j  ,i  ,1,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 flux(j  ,i  ,1,2) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(j  ,i  ,1,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 flux(j  ,i  ,1,2) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j  ,i  ,1)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + add_source(j  ,i  )*real_area13

        endif

        ! Update Temperature by Eulerian method 
        temp(j,i) = temp(j,i)+rhs*dt/area_n
    end do
end do
!$OMP end do
!$ACC end parallel
!$ACC end data

! Boundary conditions (top and bottom)
!$ACC parallel loop
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
!$ACC end parallel

! Boundary conditions: dt/dx =0 on left and right  
!$ACC parallel loop
!$OMP do
do j = 1,nz
    temp(j ,1)  = temp(j,2)
    temp(j, nx) = temp(j,nx-1)
end do
!$OMP end do
!$OMP end parallel
!$ACC end parallel

! ! HOOK
! ! Intrusions - see user_ab.f90
! if( if_intrus.eq.1 ) call MakeIntrusions

return
end 
