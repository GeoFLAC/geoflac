! ---------------------------------------------------------------
! Calculation of coordinates at the start
! ---------------------------------------------------------------

subroutine init_cord
use arrays
use params
implicit none
integer :: i, j, u, ios
character(len=256) :: iomsg
double precision :: rmesh1(nx+nz), x, y, amp

if (ircoord .gt. 0) then
    ! read the coordinates from file
    open(newunit=u, file=coordfile, status='old', iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
        call SysMsg('INIT_CORD: Cannot open file with initial coordinates: ' // trim(iomsg))
        stop 21
    endif

    do i = 1,nx
        do j = 1,nz
            read(u, *, iostat=ios) x, y
            if (ios /= 0) then
                call SysMsg('INIT_CORD: Error reading file with initial coordinates!')
                stop 21
            endif
            cord(j,i,1) = x
            cord(j,i,2) = y
        end do
    end do
    close(u)
else
    ! Check dimensions for the mesh (number of elements and sizes)  
    call test_grid
    !  X - component 
    call mesh1 (x0,x0+rxbo,rmesh1,nx,nzonx,nelz_x,sizez_x) 

    do i = 1,nx
        cord (:,i,1) =  rmesh1(i)
    end do

    !  Z - component
    do j = 1,nz
        call mesh1 (z0,z0+rzbo,rmesh1,nz,nzony,nelz_y,sizez_y)
        cord (j,:,2) =  rmesh1(j)
    end do

    ! topo perturbation
    do i = 1, inhom
        if (igeom(i).eq.20) then
            ! table mountain
            amp = xinitaps(i)
            cord(1, ix1(i):ix2(i), 2) = cord(1, ix1(i):ix2(i), 2) + amp
        elseif (igeom(i).eq.21) then
            ! trapzoidal mountain
            amp = xinitaps(i)
            do j = ix1(i), ix2(i)-1
                cord(1,j,2) = cord(1,j,2) + amp*real(j-ix1(i))/(ix2(i)-ix1(i))
            enddo
            cord(1, ix2(i):iy1(i), 2) = cord(1, ix2(i):iy1(i), 2) + amp
            do j = iy1(i)+1, iy2(i)
                cord(1,j,2) = cord(1,j,2) + amp*real(iy2(i)-j)/(iy2(i)-iy1(i))
            enddo 
        endif
    enddo
    ! Re-distribute every column from (new top) to (fixed bottom) so the
    ! initial vertical mesh matches what rem_cord produces at runtime.
    ! Without this, the first rem_cord call after t=0 stretches the column
    ! and drags depth-based initial phase layers up by ~amp.
    do i = 1, nx
        call mesh1(cord(1,i,2), cord(nz,i,2), rmesh1, nz, nzony, nelz_y, sizez_y)
        do j = 1, nz
            cord(j,i,2) = rmesh1(j)
        end do
    end do
endif

! min element width and thickness
dxmin = minval(cord(1,2:nx,1) - cord(1,1:nx-1,1))
dzmin = minval(cord(1:nz-1,1,2) - cord(2:nz,1,2))

!$ACC kernels async(1)
dhacc(1:nx-1) = 0.d0
dhacc_correction(1:nx) = 0.d0
extr_acc(1:nx-1) = 0.d0
!$ACC end kernels

return
end



!====================================================================
!------Test for mesh dimensions (check size and num. of elem)

subroutine test_grid
use params
implicit none
integer :: i, nelsum
!- X direction
nelsum =0
do i = 1,nzonx
    nelsum =  nelsum  + nelz_x(i)
end do

if ( nzonx.ne.0 .and. mod(nzonx,2) .eq. 0 ) then
    call SysMsg('INIT_CORD: Number of zones is not odd.(X-direction)')
    stop 16
endif

if (nelsum .ne. (nx-1)) then
    call SysMsg('INIT_CORD: Nelem in zones is not correct.(X-direction)')
    stop 17
endif


!- Y direction
nelsum =0
do i = 1,nzony
    nelsum =  nelsum  + nelz_y(i)
end do

if ( nzony.ne.0 .and. mod(nzony,2) .eq. 0 ) then
    call SysMsg('INIT_CORD: Number of zones is not odd.(Y-direction)')
    stop 18
endif

if (nelsum .ne. (nz-1)) then
    call SysMsg('INIT_CORD: Nelem in zones is not correct.(Y-direction)')
    stop 19
endif

return
end 


!==========================================================================
! ONE-D mesh
!==========================================================================
subroutine mesh1 (x1,x2,xmesh,n,nzon,nelz,sizez)
!$ACC routine seq
implicit none
integer, intent(in) :: n, nzon, nelz(nzon)
double precision, intent(in) :: x1, x2, sizez(nzon)
double precision, intent(out) :: xmesh(n)
integer :: k, i, nel
double precision :: xbeg, xend, elemsize, all_ratio, magnification, ratio
all_ratio = 0.d0
nel = 1
xmesh(nel) = 0.d0

do k = 1,nzon-1,2
    all_ratio = all_ratio + nelz(k) * sizez(k)
    do i = 1,nelz(k)
        nel = nel + 1
        xmesh(nel) = xmesh(nel-1) + sizez(k)
    end do
    ratio = (sizez(k+2)/sizez(k))**(1.0d0/dble(nelz(k+1)+1))
    do i = 1,nelz(k+1)
        nel = nel + 1
        xmesh(nel) = xmesh(nel-1) + sizez(k) * (ratio**i)
        all_ratio = all_ratio + sizez(k) * (ratio**i)
    end do
end do

do k = 1,nelz(nzon)
    nel = nel + 1
    xmesh(nel) = xmesh(nel-1) + sizez(nzon)
end do
all_ratio = all_ratio + nelz(nzon) * sizez(nzon)
magnification = (x2 - x1)/all_ratio

do nel = 1,n
    xmesh(nel) = x1 + xmesh(nel) * magnification
end do

return
end
