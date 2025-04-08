subroutine bar2euler
use arrays
use params
USE marker_data
implicit none

double precision :: ba1, ba2, ba3, x, y
integer :: i, ii, i1, i2, i3, jj, j1, j2, j3, k, n

! calculate the new paramters for the triangles
!$ACC parallel loop async(1)
!$OMP parallel do private(ba1, ba2, ba3, x, y, i, ii, i1, i2, i3, jj, j1, j2, j3, k, n)
do i = 1 , nmarkers
  if (mark_dead(i).eq.0) cycle
  n = mark_ntriag(i)
  k = mod(n-1, 2) + 1
  jj = mod((n - k) / 2, nz-1) + 1
  ii = (n - k) / 2 / (nz - 1) + 1
  ba1 = mark_a1(i)
  ba2 = mark_a2(i)
  ba3 = 1.0d0 - ba1 - ba2

  ! assign the mel values and aps to markers in element
  mark_zpressn(i) = zpressm(jj,ii)
  mark_Emeltm(i) = Eff_melt(jj,ii)

  if (k .eq. 1) then
    i1 = ii
    i2 = ii
    i3 = ii + 1
    j1 = jj
    j2 = jj + 1
    j3 = jj
  else
    i1 = ii + 1
    i2 = ii
    i3 = ii + 1
    j1 = jj
    j2 = jj + 1
    j3 = jj + 1
  endif

  ! interpolate nodal values to the marker
  x = cord(j1,i1,1)*ba1 + cord(j2,i2,1)*ba2 + cord(j3,i3,1)*ba3
  y = cord(j1,i1,2)*ba1 + cord(j2,i2,2)*ba2 + cord(j3,i3,2)*ba3
  mark_x(i) = x
  mark_y(i) = y
enddo
!$OMP end parallel do
return
end subroutine bar2euler


subroutine shape_functions(j, i, shp2)
  !$ACC routine seq
  use arrays

  include 'precision.inc'
  double precision, intent(out) :: shp2(2,3,2)
  do k = 1 , 2
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
      shp2(1,1,k) = (x2*y3-y2*x3)/det
      shp2(1,2,k) = (y2-y3)/det
      shp2(1,3,k) = (x3-x2)/det
      shp2(2,1,k) = (x3*y1-y3*x1)/det
      shp2(2,2,k) = (y3-y1)/det
      shp2(2,3,k) = (x1-x3)/det
  enddo
end subroutine shape_functions


! For brevity shp(1,2) --> s12 etc
!
! a1 = s11 + s12*x + s13*y
! a2 = s21 + s22*x + s23*y
!
! solve for x and y
subroutine bar2xy(ba1, ba2, shp, x, y)
  !$ACC routine seq
  use arrays

  include 'precision.inc'
  double precision :: shp(2,3)

  xnum = ba2*shp(1,3) - shp(2,1)*shp(1,3) - shp(2,3)*ba1 + shp(2,3)*shp(1,1)
  xdem = shp(1,3)*shp(2,2) - shp(2,3)*shp(1,2)
  x = xnum / xdem
  y = (ba1 - shp(1,1) - shp(1,2)*(xnum/xdem)) / shp(1,3)
end subroutine bar2xy
