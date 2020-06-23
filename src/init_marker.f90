subroutine init_marker

USE marker_data

use arrays
use params
implicit none
double precision :: a(3,2), b(3,2), points(9,2)
double precision, parameter :: half = 0.5d0
double precision, parameter :: onesixth = 0.1666666666666666666666d0
double precision, parameter :: fivesixth = 0.8333333333333333333333d0
integer :: i, j, i1, i2, iamp, inc, itop, iwidth, ixc, k, k1, k2, kph, n, l
double precision :: ddx, ddy, dx, dy, r, rx, ry, xx, ycol1, ycol2, ycol3, ycol4, &
                    yy, yyy

! define euler coordinate of the markers
! Distribute evenly first then randomize the distribution
! to start 9 markers per elements
mark_id_elem = 0
nmark_elem = 0
nmarkers = 0

! zones with 9 markers per elements
! calculate the id (element number) of the zones of high res

!call random_seed
!write(333,*) 'Call to random_seed(), result may be stochastic'

do i = 1 , nx-1
    do j = 1 , nz-1
        ! Alog the edge of an element; a and b are the nodes
        !   a - o -- x -- v - b
        !
        ! x is located at the midpint
        !   x = a / 2 + b / 2
        !
        ! o is located at a distance of 1/6 length
        !   o = 5/6 * a + 1/6 * b
        !
        ! v is located at a distance of 5/6 length
        !   o = 1/6 * a + 5/6 * b
        !
        ! Considering two elements
        !   a - o1 -- x1 -- v1 - b - o2 -- x2 -- v2 - c
        !
        ! o1, x1, v1, o2, x2, v2 will be equi-distant
        !

        a(1,:) = cord(j,i,:)*fivesixth + cord(j+1,i,:)*onesixth
        a(2,:) = cord(j,i,:)*half + cord(j+1,i,:)*half
        a(3,:) = cord(j,i,:)*onesixth + cord(j+1,i,:)*fivesixth

        b(1,:) = cord(j,i+1,:)*fivesixth + cord(j+1,i+1,:)*onesixth
        b(2,:) = cord(j,i+1,:)*half + cord(j+1,i+1,:)*half
        b(3,:) = cord(j,i+1,:)*onesixth + cord(j+1,i+1,:)*fivesixth

        points(1,:) = a(1,:)*fivesixth + b(1,:)*onesixth
        points(2,:) = a(1,:)*half + b(1,:)*half
        points(3,:) = a(1,:)*onesixth + b(1,:)*fivesixth

        points(4,:) = a(2,:)*fivesixth + b(2,:)*onesixth
        points(5,:) = a(2,:)*half + b(2,:)*half
        points(6,:) = a(2,:)*onesixth + b(2,:)*fivesixth

        points(7,:) = a(3,:)*fivesixth + b(3,:)*onesixth
        points(8,:) = a(3,:)*half + b(3,:)*half
        points(9,:) = a(3,:)*onesixth + b(3,:)*fivesixth

        dx = cord(j,i+1,1) - cord(j,i,1)
        dy = cord(j+1,i,2) - cord(j,i,2)

! randomize the new coordinates inside the element
        l = 1
        do while (l .le. 9)
            ! initialize kph in case nzone_age = 0
            kph = iphase(j,i)

            ! position of the marker
            call random_number(rx)
            call random_number(ry)
            rx = 0.5d0 - rx
            ry = 0.5d0 - ry
            ddx = dx*rx/3
            ddy = dy*ry/3
            xx = points(l,1) + ddx
            yy = points(l,2) + ddy

            ! phase of the marker
            ! smooth transition of marker phase
            do n = 1, nzone_age
                if (i<ixtb1(n) .or. i>ixtb2(n)) cycle
                if (n /= 1) then
                    if (iph_col_trans(n-1) == 1) cycle
                endif
                ycol1 = hc1(n)
                ycol2 = hc2(n)
                ycol3 = hc3(n)
                ycol4 = hc4(n)
                if (iph_col_trans(n) == 1) then
                    i1 = ixtb1(n)
                    i2 = ixtb2(n)
                    r = (cord(1,i,1) - cord(1,i1,1)) / (cord(1,i2,1) - cord(1,i1,1))
                    ycol1 = hc1(n) + (hc1(n+1) - hc1(n)) * r
                    ycol2 = hc2(n) + (hc2(n+1) - hc2(n)) * r
                    ycol3 = hc3(n) + (hc3(n+1) - hc3(n)) * r
                    ycol4 = hc4(n) + (hc4(n+1) - hc4(n)) * r
                endif

                ! layer
                yyy = yy * (-1d-3)

                if (yyy.lt.ycol1) then
                    kph = iph_col1(n)
                else if (yyy.lt.ycol2) then
                    kph = iph_col2(n)
                else if (yyy.lt.ycol3) then
                    kph = iph_col3(n)
                else if (yyy.lt.ycol4) then
                    kph = iph_col4(n)
                else
                    kph = iph_col5(n)
                end if
                exit
            end do

            call add_marker(xx, yy, kph, 0.d0, nmarkers, j, i, inc)
            if(inc.eq.0) cycle

            l = l + 1
            !print *, xx, yy, mark_a1(kk), mark_a2(kk), mark_ntriag(kk)
        enddo
    enddo
enddo

!   Put initial heterogeneities
do i = 1,inhom
    if (inphase(i) < 0) cycle  ! skip

    ! Rectangular shape:
    if (igeom(i) .eq.0) then
        call newphase2marker(iy1(i),iy2(i),ix1(i),ix2(i),inphase(i))
    endif

    ! Gauss shape:
    if (igeom(i).eq.1.or.igeom(i).eq.2) then
        stop 1
    endif

    ! weak zone at 45 degree
    if (igeom (i) .eq.3) then
        do j = ix1(i),ix2(i)
            k = nint(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            call newphase2marker(k,k,j,j,inphase(i))
        end do
    endif

    ! Weak zone in accumulated plastic strain at 45 degree
    if (igeom (i).eq.4) then
        do j =ix1(i),ix2(i)
            k1 = floor(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            k2 = ceiling(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            if( inphase(i) .ge. 0) then
                call newphase2marker(k1,k2,j,j,inphase(i))
            endif
        end do
    endif
end do

write(333,*) '# of markers', nmarkers

!!! This is the first GPU code during initialization
!$ACC update device(nmark_elem, mark_id_elem)
call count_phase_ratio_all
!$ACC update host(phase_ratio, iphase)
!!! Part of initialization later still runs on CPU
return
end subroutine init_marker
