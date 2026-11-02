! Reconstruct markers' Euler (absolute x,y) coordinates from their
! barycentric coordinates (mark_a1, mark_a2) within their current
! triangle (mark_ntriag) and that triangle's current (deformed) node
! positions. Euler coordinates are not stored persistently; this is how
! they get (re)computed wherever needed.
!
! bar2euler (no args, called right before lpeuler2bar in par.f90's remesh
! path): results are written into mark_a1/mark_a2 in place. This is safe
! there specifically because lpeuler2bar unconditionally overwrites
! mark_a1/mark_a2 with newly-computed barycentric coordinates immediately
! after, so the old barycentric values are about to be discarded anyway --
! using them as scratch space for x,y avoids needing a separate buffer.
!
! bar2euler_xy (outmarker.f90, and anywhere else that still needs the true
! mark_a1/mark_a2 afterward): results go into the caller's own buffer
! instead, leaving mark_a1/mark_a2 untouched.
!
! These are two separate subroutines rather than one with an optional
! argument: optional and assumed-shape dummy arguments require an explicit
! interface at the call site to work correctly, but this codebase calls
! every subroutine (including this one) without module wrapping or
! interface blocks, so such arguments would be passed under undefined
! (F77-style implicit-interface) semantics. Both share the per-marker
! interpolation logic via marker2xy below.
subroutine bar2euler
!$ACC routine(marker2xy) seq
use arrays
use params
USE marker_data
implicit none

double precision :: x, y
integer :: i

! calculate the new paramters for the triangles
!$ACC parallel loop async(1) private(x, y)
!$OMP parallel do private(i, x, y)
do i = 1 , nmarkers
  if (mark_dead(i).eq.0) cycle
  call marker2xy(i, x, y)
  mark_a1(i) = x
  mark_a2(i) = y
enddo
!$OMP end parallel do
return
end subroutine bar2euler


subroutine bar2euler_xy(xout, yout)
!$ACC routine(marker2xy) seq
use arrays
use params
USE marker_data
implicit none

double precision, intent(out) :: xout(nmarkers), yout(nmarkers)
integer :: i

!$ACC parallel loop async(1)
!$OMP parallel do private(i)
do i = 1 , nmarkers
  if (mark_dead(i).eq.0) cycle
  call marker2xy(i, xout(i), yout(i))
enddo
!$OMP end parallel do
return
end subroutine bar2euler_xy


! Euler (x,y) position of marker i, interpolated from its barycentric
! coordinates (mark_a1(i), mark_a2(i)) within its current triangle
! (mark_ntriag(i)) and that triangle's current (deformed) node positions.
subroutine marker2xy(i, x, y)
  !$ACC routine seq
  use arrays
  use params
  USE marker_data
  implicit none

  integer, intent(in) :: i
  double precision, intent(out) :: x, y
  double precision :: ba1, ba2, ba3
  integer :: ii, i1, i2, i3, jj, j1, j2, j3, k, n

  n = mark_ntriag(i)
  k = mod(n-1, 2) + 1
  jj = mod((n - k) / 2, nz-1) + 1
  ii = (n - k) / 2 / (nz - 1) + 1
  ba1 = mark_a1(i)
  ba2 = mark_a2(i)
  ba3 = 1.0d0 - ba1 - ba2

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

  ! adjust y for non-tectonic topographic accumulation on surface nodes
  if (j1 .eq. 1) y = y - dhacc_correction(i1)*ba1
  if (j2 .eq. 1) y = y - dhacc_correction(i2)*ba2
  if (j3 .eq. 1) y = y - dhacc_correction(i3)*ba3
end subroutine marker2xy


subroutine shape_functions(j, i, shp2)
  !$ACC routine seq
  use arrays
  implicit none

  integer, intent(in) :: j, i
  double precision, intent(out) :: shp2(2,3,2)
  double precision :: x1, x2, x3, y1, y2, y3, det
  integer :: k

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
  implicit none

  double precision, intent(in) :: ba1, ba2, shp(2,3)
  double precision, intent(out) :: x, y
  double precision :: xnum, xdem

  xnum = ba2*shp(1,3) - shp(2,1)*shp(1,3) - shp(2,3)*ba1 + shp(2,3)*shp(1,1)
  xdem = shp(1,3)*shp(2,2) - shp(2,3)*shp(1,2)
  x = xnum / xdem
  y = (ba1 - shp(1,1) - shp(1,2)*(xnum/xdem)) / shp(1,3)
end subroutine bar2xy
