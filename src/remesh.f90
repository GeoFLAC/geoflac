subroutine remesh
use arrays
use params
use marker_data
implicit none

integer :: nzt,nxt
integer :: i, idist, ii, j, k, l, m, n, iph, ju, jl
double precision :: xl, xr, densT, dh, dh1, dh2, dp, dpt, &
                    press, rogh, tmpr

! Save old mesh for interpolations
!$ACC kernels async(1)
cordo(:,:,:) = cord(:,:,:)
!$ACC end kernels

! Create The New grid (cord) using cordo(nz,i,2) for the bottom and cordo(1,i,2) for the surface
call rem_cord



!$ACC kernels async(1)
!===========================================================================!
! Length-weighted interpolation of accumulated topo change on the surface
dhnew(:) = 0
extnew(:) = 0

i = 1 ! index of new mesh
ii = 1 ! index of old mesh
xl = cord(1,i,1)
! loop over each old element
do while (max(i,ii) < nx)
    xr = min(cordo(1,ii+1,1), cord(1,i+1,1))
    if (xl < xr) then
        ! integration
        dhnew(i) = dhnew(i) + dhacc(ii) * (xr - xl)
        extnew(i) = extnew(i) + extr_acc(ii) * (xr - xl)
    endif

    if (xr == cordo(1,ii+1,1)) then
        ii = ii + 1
    else
        i = i + 1
    endif
    if (max(i,ii) < nx) then
        xl = xr
        xr = min(cordo(1,ii+1,1), cord(1,i+1,1))
    endif
enddo
! divided by segment length
dhacc(1:nx-1) = dhnew / (cord(1,2:nx,1) - cord(1,1:nx-1,1))
extr_acc(1:nx-1) = extnew / (cord(1,2:nx,1) - cord(1,1:nx-1,1))
!$ACC end kernels


!===========================================================================!
! REMESHING FOR ELEMENT-WISE PROPERTIES
! Linear interpolation in baricentric coordinates defined as centers of old mesh
nxt = nx-1
nzt = nz-1

! Old mesh - old-element centers
! New mesh - new-element centers
!$OMP parallel do
!$ACC parallel loop collapse(2) async(1)
do i = 1, nx-1
    do j = 1, nz-1
        cold(j,i,1) = 0.25d0*( cordo(j,i,1)+cordo(j+1,i,1)+cordo(j,i+1,1)+cordo(j+1,i+1,1) )
        cold(j,i,2) = 0.25d0*( cordo(j,i,2)+cordo(j+1,i,2)+cordo(j,i+1,2)+cordo(j+1,i+1,2) )
        cnew(j,i,1) = 0.25d0*( cord(j,i,1)+cord(j+1,i,1)+cord(j,i+1,1)+cord(j+1,i+1,1) )
        cnew(j,i,2) = 0.25d0*( cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2) )
    enddo
enddo
!$OMP end parallel do

! Calculate parameters of old-mesh triangles
call rem_trpars(nzt, nxt)

! Baricentric coordinates of new-elements centers
call rem_barcord(nzt, nxt)

! Do interpolations

! Interpolate Stress (in quadralaterals)
do k = 1,4
    do l = 1,4
        call rem_interpolate( nzt, nxt, dummye, stress0(:,:,k,l) )
    end do
end do

! HOOK
! Remeshing mode 3 - see user_luc.f90
!if( mode_rem .eq. 3 ) call rem_stress_alt()


! Interpolate strains
do k = 1, 3
    call rem_interpolate( nzt, nxt, dummye, strain(:,:,k) )
end do


! plastic strain
call rem_interpolate( nzt, nxt, dummye, aps )

! Magma fraction
if (itype_melting == 1) then
    call rem_interpolate( nzt, nxt, dummye, fmagma )
    call rem_interpolate( nzt, nxt, dummye, fmagma2 )
    ! Dike-intrusion eigenstrain/marker tracking (persistent parts only;
    ! dike_added/dike_released are same-timestep scratch, recomputed fresh
    ! every fl_therm.f90 call and never needed across a remesh event)
    call rem_interpolate( nzt, nxt, dummye, dike_backlog )
    call rem_interpolate( nzt, nxt, dummye, dike_marker_vol )
    call rem_interpolate( nzt, nxt, dummye, mor_marker_vol )
endif


!$OMP parallel do
!$ACC parallel loop collapse(2) async(1)
do i = 1, nxt
    do j = 1, nzt
        if( aps(j,i) .lt. 0.d0 ) then
            aps(j,i) = 0.d0
        endif
!       write(*,*) i,j,aps(j,i)
    end do
end do
!$OMP end parallel do

! phases
! XXX: Assuming material is coming from the left or right boundary
! and is within 3-element distance from the boundary and has the
! same layered structure as the 4th element away from the boundary

idist = 2
ju = 1
jl = 1
!$ACC data copyin(idist,ju,jl)
if(incoming_left==1) then
    !$ACC kernels async(1)
    aps(1:nz-1, 1:1+idist) = 0.0d0
    !$ACC end kernels

    n = 1
    i = 1
    do m = 1, nph_layer(n)-1
        !$ACC serial async(1)
        ju = jl
        do j = jl, nz-1
            if((cord(1,i,2) - 0.5d0*(cord(j,i,2)+cord(j+1,i,2))) > hc(n,m)*1d3) then
                jl = j
                exit
            endif
        enddo
        !$ACC end serial

        !$ACC parallel async(1)
        call newphase2marker(ju, jl-1, i, i+idist, iph_col(n,m))
        !$ACC end parallel
    enddo

    !$ACC parallel async(1)
    call newphase2marker(jl, nz-1, i, i+idist, iph_col(n,nph_layer(n)))
    !$ACC end parallel
endif

if(incoming_right==1) then
    !$ACC kernels async(1)
    aps(1:nz-1, nx-1-idist:nx-1) = 0.0d0
    !$ACC end kernels

    n = nzone_age
    i = nx - 1
    do m = 1, nph_layer(n)-1
        !$ACC serial async(1)
        ju = jl
        do j = jl, nz-1
            if((cord(1,i+1,2) - 0.5d0*(cord(j,i+1,2)+cord(j+1,i+1,2))) > hc(n,m)*1d3) then
                jl = j
                exit
            endif
        enddo
        !$ACC end serial

        !$ACC parallel async(1)
        call newphase2marker(ju, jl-1, i-idist, i, iph_col(n,m))
        !$ACC end parallel
    enddo

    !$ACC parallel async(1)
    call newphase2marker(jl, nz-1, i-idist, i, iph_col(n,nph_layer(n)))
    !$ACC end parallel
endif
!$ACC end data

!!$! XXX: the bottom elements must be mantle material, otherwise
!!$! too much deformation can occur(?)
!!$aps(nz-3:nz-1, 1:nx-1) = 0.0d0
!!$call newphase2marker(nz-3, nz-1, 1, nx-1, kmant1)


! sources
call rem_interpolate( nzt, nxt, dummye, source )


! REMESHING FOR NODE-WISE PROPERTIES
! Linear interpolation in baricentric coordinates of old mesh
nxt = nx
nzt = nz

!$OMP task
!$ACC kernels async(1)
! Old mesh - old coordinates points
cold(1:nz,1:nx,1:2) = cordo(1:nz,1:nx,1:2)
temp0(1:nz,1:nx) = temp(1:nz,1:nx)
! New mesh - new coordinates points
cnew(1:nz,1:nx,1:2) = cord(1:nz,1:nx,1:2)
!$ACC end kernels
!$OMP end task

! Calculate parameters of triangles of this mesh
call rem_trpars(nzt,nxt)

! Baricentric coordinates of new-elements centers
call rem_barcord(nzt,nxt)

! Do node-wise interpolations

! Velocities (in nodes)
do k = 1, 2
    call rem_interpolate( nzt, nxt, dummyn, vel(:,:,k) )
end do

! Temperatures (in nodes) 
call rem_interpolate( nzt, nxt, dummyn, temp )

! Changing the temperature of left-/right- most elements
! in accordance to initial temperature
if(incoming_left==1) call sidewalltemp(1,1+idist)
if(incoming_right==1) call sidewalltemp(nx-idist,nx)

call rem_interpolate( nzt, nxt, dummyn, xoriginal )
call rem_interpolate( nzt, nxt, dummyn, zoriginal )

! AFTER INTERPOLATIONS - RECALCULATE SOME DEPENDENT VARIABLES
if( ivis_present.eq.1 ) call init_visc

! Calculation of areas of triangle
call init_areas

! Distribution of masses in nodes
call rmasses

! 1) Determine Inertial Masses with a given dt_scale (idt_scale=1)
! 2) dt_mech with a given  Real Masses (idt_scale = 0)
call dt_mass

return
end




!===============================================
! parameters of triangles of a grid
!===============================================
subroutine rem_trpars(nzt,nxt)
use arrays
implicit none
integer, intent(in) :: nzt,nxt
integer :: i, j, k, n
double precision :: x1, x2, x3, y1, y2, y3, det

!$OMP parallel do private(x1,x2,x3,y1,y2,y3,det,n)
!$ACC parallel loop collapse(3) async(1)
do i = 1,nxt-1
    do j = 1,nzt-1
        do k = 1,2

            !  diagonal / :
            !   ii=1     ii=2 

            !  1---3         1
            !  | /         / |
            !  2         2---3

            if (k.eq.1) then
                x1 = cold(j  ,i  ,1)
                x2 = cold(j+1,i  ,1)
                x3 = cold(j  ,i+1,1)
                y1 = cold(j  ,i  ,2)
                y2 = cold(j+1,i  ,2)
                y3 = cold(j  ,i+1,2)
            else  !if (k.eq.2) then 
                x1 = cold(j  ,i+1,1) 
                x2 = cold(j+1,i  ,1) 
                x3 = cold(j+1,i+1,1)
                y1 = cold(j  ,i+1,2) 
                y2 = cold(j+1,i  ,2) 
                y3 = cold(j+1,i+1,2) 
            endif 
 
            det=( (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2) )

            n = 2*( (nzt-1)*(i-1)+j-1 ) + k

            !Find the parameters ONLY for 2 vertices
            pt(n,1,1)=(x2*y3-y2*x3)/det
            pt(n,1,2)=(y2-y3)/det
            pt(n,1,3)=(x3-x2)/det
            pt(n,2,1)=(x3*y1-y3*x1)/det
            pt(n,2,2)=(y3-y1)/det
            pt(n,2,3)=(x1-x3)/det
        end do
    end do
end do
!$OMP end parallel do

return
end



!===============================================
! baricentric coordinates of new mesh in old triangles
!===============================================
subroutine rem_barcord(nzt,nxt)
use arrays
implicit none
integer, intent(in) :: nzt,nxt
integer :: i, j, k, l, lt, io, jo, m, n, mm, nmin, numqu
double precision :: perr, xx, yy, a1, a2, a3, amod, amodmin, dist1, dist2, dist3


perr = 1.d-4
mm = 4  ! search range, might be extended if find grid / coarse grid is large

!$OMP parallel do private(xx,yy,l, lt,io,jo,k,n,a1,a2,a3,amod,amodmin,nmin, &
!$OMP                     numqu,dist1,dist2,dist3)
!$ACC parallel loop collapse(2) async(1)
do i = 1, nxt
    do j = 1, nzt
        xx = cnew(j,i,1)
        yy = cnew(j,i,2)

        amodmin = 1.d+10

        numtr(j,i) = 0
        search_loop: do io = max(1, i-mm), min(nxt-1, i+mm)
            do jo = max(1, j-mm), min(nzt-1, j+mm)
                do k = 1, 2
                    n = 2*( (nzt-1)*(io-1)+jo-1 ) + k
                    a1 = pt(n,1,1) + xx*pt(n,1,2) + yy*pt(n,1,3)
                    a2 = pt(n,2,1) + xx*pt(n,2,2) + yy*pt(n,2,3)
                    a3 = 1.d0 - a1 - a2
                    amod = a1*a1 + a2*a2 + a3*a3
                    if( amod .lt. amodmin ) then
                        amodmin = amod
                        nmin = n
                    endif
                    if ( a1.le.-perr .or. a1.ge.(1.0d0+perr) ) cycle
                    if ( a2.le.-perr .or. a2.ge.(1.0d0+perr) ) cycle
                    if ( a3.le.-perr .or. a3.ge.(1.0d0+perr) ) cycle
                    numtr(j,i) = n
                    barcord(j,i,1) = a1
                    barcord(j,i,2) = a2
                    barcord(j,i,3) = a3
                    exit search_loop
                end do

            end do
        end do search_loop
    
    ! if a point is outside of the mesh
    if( numtr(j,i) .eq. 0 ) then
        numtr(j,i) = nmin
        numqu = (numtr(j,i) + 1)/2          
        io = (numqu-1)/(nzt-1) + 1
        jo = numqu - (io-1)*(nzt-1)

        if( mod(numtr(j,i),2).eq.0 ) then
            dist1 = sqrt((xx-cold(jo,io+1,1))**2 + (yy-cold(jo,io+1,2))**2)
            dist2 = sqrt((xx-cold(jo+1,io,1))**2 + (yy-cold(jo+1,io,2))**2)
            dist3 = sqrt((xx-cold(jo+1,io+1,1))**2 + (yy-cold(jo+1,io+1,2))**2)
        else
            dist1 = sqrt((xx-cold(jo,io,1))**2 + (yy-cold(jo,io,2))**2)
            dist2 = sqrt((xx-cold(jo+1,io,1))**2 + (yy-cold(jo+1,io,2))**2)
            dist3 = sqrt((xx-cold(jo,io+1,1))**2 + (yy-cold(jo,io+1,2))**2)
        endif
        
        ! select two nearest points for interpolation
        if( dist1 .gt. dist2 ) then
            if( dist2 .gt. dist3 ) then
                ! 1,2,3
                barcord(j,i,1) = 0.d0
                barcord(j,i,2) = 1.d0/dist2/(1.d0/dist2+1.d0/dist3)
                barcord(j,i,3) = 1.d0/dist3/(1.d0/dist2+1.d0/dist3)
            else
                if( dist1 .gt. dist3 ) then
                    ! 1,3,2
                    barcord(j,i,1) = 0.d0
                    barcord(j,i,2) = 1.d0/dist2/(1.d0/dist2+1.d0/dist3)
                    barcord(j,i,3) = 1.d0/dist3/(1.d0/dist2+1.d0/dist3)
                else
                    ! 3,1,2
                    barcord(j,i,1) = 1.d0/dist1/(1.d0/dist1+1.d0/dist2)
                    barcord(j,i,2) = 1.d0/dist2/(1.d0/dist1+1.d0/dist2)
                    barcord(j,i,3) = 0.d0
                endif
            endif
        else
            if( dist3 .gt. dist1 ) then
                if( dist2 .gt. dist3 ) then
                    ! 2,3,1
                    barcord(j,i,1) = 1.d0/dist1/(1.d0/dist1+1.d0/dist3)
                    barcord(j,i,2) = 0.d0
                    barcord(j,i,3) = 1.d0/dist3/(1.d0/dist1+1.d0/dist3)
                else
                    ! 3,2,1
                    barcord(j,i,1) = 1.d0/dist1/(1.d0/dist1+1.d0/dist2)
                    barcord(j,i,2) = 1.d0/dist2/(1.d0/dist1+1.d0/dist2)
                    barcord(j,i,3) = 0.d0
                endif
            else
                ! 2,1,3
                barcord(j,i,1) = 1.d0/dist1/(1.d0/dist1+1.d0/dist3)
                barcord(j,i,2) = 0.d0
                barcord(j,i,3) = 1.d0/dist3/(1.d0/dist1+1.d0/dist3)
            endif
        endif

    endif

    end do
end do
!$OMP end parallel do

return
end


!===============================================
! interpolation
!===============================================
subroutine rem_interpolate( nzt, nxt, dummy, arr )
use arrays
implicit none
integer, intent(in) :: nzt,nxt
double precision, intent(out) :: dummy(nzt,nxt)
double precision, intent(inout) :: arr(nzt,nxt)
integer :: i, j, io, jo, numq
double precision :: f1, f2, f3

!$ACC kernels async(1)
dummy = arr
!$ACC end kernels

!$OMP parallel do private(numq,io,jo,f1,f2,f3)
!$ACC parallel loop collapse(2) async(1)
do i = 1, nxt
    do j = 1, nzt

        numq = (numtr(j,i)+1) / 2
        io = (numq-1)/(nzt-1) + 1
        jo = numq - (io-1)*(nzt-1)

        !  diagonal / :
        !   kk=1       kk=2
        !
        !  1---3         1
        !  | /         / |
        !  2         2---3

        if( mod(numtr(j,i),2).eq.0 ) then
            f1 = dummy(jo  ,io+1)
            f2 = dummy(jo+1,io  )
            f3 = dummy(jo+1,io+1)
        else
            f1 = dummy(jo  ,io  )
            f2 = dummy(jo+1,io  )
            f3 = dummy(jo  ,io+1)
        endif

        arr(j,i) = barcord(j,i,1)*f1 + barcord(j,i,2)*f2 + barcord(j,i,3)*f3
    end do
end do
!$OMP end parallel do

return

end
