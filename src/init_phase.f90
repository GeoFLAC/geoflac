!--------------------------------------------------------------
! Initialization of phases in each zone 
!--------------------------------------------------------------

subroutine init_phase
use arrays
use params
implicit none
integer :: i, j, k, ii, jj, k1, k2, u
if (irphase .gt. 0) then
    !  Read distribution of the phases from the dat file
    open(newunit=u, file=phasefile)
    do i=1,nx-1
        do j=1,nz-1
            read(u,*) ii,jj,iphase(j,i)
        enddo
    enddo
    close(u)
endif

!   Put different rheologies for inclusions 
do i = 1,inhom
    ! Rectangular shape:
    if (igeom(i) .eq.0) then
        do j = ix1(i),ix2(i)
            do k = iy1(i),iy2(i)
                if( inphase(i) > 0) iphase(k,j) = inphase(i)
                aps(k,j)=xinitaps(i)
            end do
        end do
    endif

    ! weak zone at 45 degree
    if (igeom (i) .eq.3) then
        do j = ix1(i),ix2(i)
            k = nint(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            if( inphase(i) > 0) iphase(k,j) = inphase(i)
        end do
    endif
    
    ! Weak zone in accumulated plastic strain at 45 degree        
    if (igeom (i).eq.4) then
        do j =ix1(i),ix2(i)
            k1 = floor(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            k2 = ceiling(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            aps(k1:k2,j)=xinitaps(i)
            if( inphase(i) > 0) iphase(k1:k2,j) = inphase(i)
        end do
    endif
end do

return
end


subroutine check_visc_rheol
use arrays
use params
implicit none
integer :: i, j, iph

! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(j,i)
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do

return
end
