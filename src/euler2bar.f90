subroutine euler2bar(x,y,bar1,bar2,ntr,ii,jj,inc)
!$ACC routine seq
!$ACC routine(check_inside) seq

use arrays
use params
implicit none
double precision, intent(in) :: x, y
double precision, intent(out) :: bar1, bar2
integer, intent(out) :: ntr, inc
integer, intent(inout) :: ii, jj
integer :: i, j, n1, n2, n3
!character*200 msg

! find the triangle in which the marker belongs

! check (jj,ii) elem first
i = ii
j = jj
call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
if(inc .eq. 1) return

! search the neighboring elem
n1 = 2
do j = max(1, jj-n1), min(nz-1, jj+n1)
    do i = max(1, ii-n1), min(nx-1, ii+n1)
        if (i == ii .and. j == jj) cycle  ! already checked
        call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
        if(inc .eq. 1) then
            ii = i
            jj = j
            return
        endif
    enddo
enddo

n2 = 6
do j = max(1, jj-n2), min(nz-1, jj+n2)
    do i = max(1, ii-n2), min(nx-1, ii+n2)
        if (abs(i - ii) <= n1 .and. abs(j - jj) <= n1) cycle  ! already checked
        call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
        if(inc .eq. 1) then
            ii = i
            jj = j
            return
        endif
    enddo
enddo

n3 = 12
do j = max(1, jj-n3), min(nz-1, jj+n3)
    do i = max(1, ii-n3), min(nx-1, ii+n3)
        if (abs(i - ii) <= n2 .and. abs(j - jj) <= n2) cycle  ! already checked
        call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
        if(inc .eq. 1) then
            ii = i
            jj = j
            return
        endif
    enddo
enddo

! More extensive search
! disabled since most models won't need this
if(.false.) then

! search all surface elem
do i = 1, nx-1
    call check_inside(x,y,bar1,bar2,ntr,i,1,inc)
    if(inc .eq. 1) then
        ii = i
        jj = 1
        return
    endif
enddo

! search all elem, usually it means the marker is "dead"
!write(msg,*) 'Searching all elem. for marker, original (i,j) ', ii, jj
!call SysMsg(msg)
do j = 1, nz-1
    do i = 1, nx-1
        call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
        if(inc .eq. 1) then
            ! If the marker is found, that means its current element is
            ! to far away from its original element
            !!!write(msg,*) 'Found at (i,j)', i, j, ' Might need more frequent remeshing?'
            !!!call SysMsg(msg)
            ii = i
            jj = j
            return
        endif
    enddo
enddo
endif

inc = 0
return
end subroutine euler2bar



subroutine check_inside(x,y,bar1,bar2,ntr,i,j,inc)
  !$ACC routine seq
  use arrays
  use params
  implicit none
  double precision, intent(in) :: x, y
  double precision, intent(out) :: bar1, bar2
  integer, intent(out) :: ntr, inc
  integer, intent(in) :: i, j
  integer :: k, ninters
  double precision :: x1, x2, x3, y1, y2, y3, det, xxmpt(2,3)


  inc = 0

  do k = 1 , 2

      ninters = 0
      if (k.eq.1) then
          x1 = cord(j  ,i  ,1)
          x2 = cord(j+1,i  ,1)
          x3 = cord(j  ,i+1,1)
          y1 = cord(j  ,i  ,2)
          y2 = cord(j+1,i  ,2)
          y3 = cord(j  ,i+1,2)
      else  !if (k.eq.2) then
          x1 = cord(j  ,i+1,1)
          x2 = cord(j+1,i  ,1)
          x3 = cord(j+1,i+1,1)
          y1 = cord(j  ,i+1,2)
          y2 = cord(j+1,i  ,2)
          y3 = cord(j+1,i+1,2)
      endif

      ! Calculate triangle properties
      det=( (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2) )

      !Find the parameters ONLY for 2 vertices
      xxmpt(1,1) = (x2*y3-y2*x3)/det
      xxmpt(1,2) = (y2-y3)/det
      xxmpt(1,3) = (x3-x2)/det
      xxmpt(2,1) = (x3*y1-y3*x1)/det
      xxmpt(2,2) = (y3-y1)/det
      xxmpt(2,3) = (x1-x3)/det

      ! Calculate barycentic coordinates
      bar1 = xxmpt(1,1) + xxmpt(1,2)*x + y*xxmpt(1,3)
      bar2 = xxmpt(2,1) + xxmpt(2,2)*x + y*xxmpt(2,3)

      ! found the triangle
      if(bar1 >= 0.0d0 .and. bar2 >= 0.0d0 .and. bar1+bar2 <= 1.0d0) then
          ntr = 2 * ( (nz-1)*(i-1)+j-1) + k
          inc = 1
          return
      endif

  enddo

  return
end subroutine check_inside

