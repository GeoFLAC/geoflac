
!------ calculation parameters of triangle (coef of basic functions)

! THE area(n,it) is inverse of "real" DOUBLE area (=1./det) =>
! area (n,it) ( in program) = 1./(2 real_area)
! real_area = 0.5* (1./area(n,t))

subroutine init_areas
use arrays
use params
implicit none
integer :: i, j
double precision :: x1, y1, x2, y2, x3, y3, x4, y4, det, det1
!$OMP parallel do private(j,x1,y1,x2,y2,x3,y3,x4,y4,det,det1)
do i = 1,nx-1
    do j = 1,nz-1

        ! Coordinates
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)
        
	! (1) Element A:
        det=((x2*y3-y2*x3)-(x1*y3-y1*x3)+(x1*y2-y1*x2))
        det1 = 1.d0/det
        area(j,i,1) = det1

        ! (2) Element B:
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        det1 = 1.d0/det
        area(j,i,2) = det1

        ! (3) Element C:
        det=((x2*y4-y2*x4)-(x1*y4-y1*x4)+(x1*y2-y1*x2))
        det1 = 1.d0/det
        area(j,i,3) = det1

        ! (4) Element D:
        det=((x4*y3-y4*x3)-(x1*y3-y1*x3)+(x1*y4-y1*x4))
        det1 = 1.d0/det
        area(j,i,4) = det1
        if( area(j,i,1).lt.0.or.area(j,i,2).lt.0.or.area(j,i,3).lt.0.or.area(j,i,4).lt.0 ) then
            call SysMsg('INIT_AREAS: Negative area!')
            stop 41
        endif
    enddo
enddo
!$OMP end parallel do
return
end subroutine init_areas

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


function total_area( iph )
use arrays
use params
implicit none
integer, intent(in) :: iph
double precision :: total_area
integer :: i, j
double precision :: area_t, x1, y1, x2, y2, x3, y3, x4, y4, det
area_t = 0
!$OMP parallel do private(j,x1,y1,x2,y2,x3,y3,x4,y4,det) reduction(+:area_t)
do i = 1,nx-1
    do j = 1,nz-1

        if( iph .ne. 0 ) then
            if( iphase(j,i) .ne. iph ) cycle
        endif
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)

        ! (1) Element A:
        det=((x2*y3-y2*x3)-(x1*y3-y1*x3)+(x1*y2-y1*x2))
        area_t = area_t + det/2

        ! (2) Element B:
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        area_t = area_t + det/2

    end do
end do
!$OMP end parallel do

total_area = area_t

return
end
