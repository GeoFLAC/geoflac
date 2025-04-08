
! 1) Testing for remeshing
! find the smallest angle in each 4 subtriangles in each qudralateral   
! and compare to angle_for_remeshing
 
! 2) Find the zone where it is needed to be remeshed:
! go from the left to the right and backwards
 
subroutine itest_mesh(need_remeshing)
  use arrays
  use params
  implicit none
  integer, intent(out) :: need_remeshing
  integer :: iv(4), jv(4), i, ii, imint, j, jmint, k
  double precision :: angle(3), testcr, shortening, topochanges, dx_accr, &
                      pi, raddeg, degrad, xa, xb, xxal, xxbl, ya, yb, &
                      anglemin1, anglemint, arc_extrusion_rate

  call check_nan
  need_remeshing = 0

  ! if remeshing with adding material on the sides then
  ! remesh at pre-given shortening
  ! dx_rem*dx - critical distance of shortnening
  if( mode_rem .eq. 11.or.mode_rem.eq.3 ) then
      !$ACC serial copyout(shortening) copy(need_remeshing) async(1)
      testcr = dx_rem * rxbo / (nx-1)
      shortening = abs(cord(1,nx,1) - cord(1,1,1) - rxbo)
      if ( shortening .gt. testcr ) then
          need_remeshing = 2
      endif
      !$ACC end serial
      !$ACC wait(1)
      if( need_remeshing .eq. 2 ) then
          if( dtout_screen .ne. 0 ) then
              print *, 'Remeshing due to shortening required: ', shortening
              write(333,*) 'Remeshing due to shortening required: ', shortening
          else
              call SysMsg('TEST_MESH: Remeshing due to shortening required')
          endif
          return
      endif
  end if

  ! remesh if the surface topo changes too much
  arc_extrusion_rate = 1.d0 - ratio_mantle_mzone
  if (topo_kappa > 0 .or. arc_extrusion_rate > 0) then
      !$ACC parallel loop copyout(topochanges) copy(need_remeshing) async(1)
      do i = 1, nx-1
          ! XXX: 1/3 thickness of the top-left element
          testcr = (cord(1,1,2) - cord(2,1,2)) * 0.3333d0
          topochanges = abs(dhacc(i)) + extr_acc(i)
          if (topochanges > testcr) then
              !$ACC atomic write
              need_remeshing = 3
          endif
      enddo
      !$ACC wait(1)
      if( need_remeshing .eq. 3 ) then
          if( dtout_screen .ne. 0 ) then
              print *, 'Remeshing due to surface topo changes required: ', topochanges
              write(333,*) 'Remeshing due to surface topo changes required: ', topochanges
          else
              call SysMsg('TEST_MESH: Remeshing due to surface topo changes required')
          endif
          return
      endif
  endif

  pi = 3.14159265358979323846d0
  degrad = pi/180.d0
  raddeg = 180.d0/pi
  anglemint = 180.d0
  imint = 0
  jmint = 0

  if (need_remeshing == 0) then
    !$OMP parallel do collapse(3) reduction(min:anglemint) private(iv,jv,angle,anglemin1,xa,xb,xxal,xxbl,ya,yb)
    !$ACC parallel loop collapse(3) reduction(min:anglemint) private(iv,jv,angle) async(1)
    do i = 1, nx-1
        do j = 1,nz-1
            ! loop for each 4 sub-triangles
            do ii = 1,4
                if (ii.eq.1) then
                    iv(1) = i ; jv(1) = j ; iv(2) = i ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j
                elseif (ii.eq.2) then
                    iv(1) = i ; jv(1) = j+1 ; iv(2) = i+1 ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j
                elseif (ii.eq.3) then
                    iv(1) = i ; jv(1) = j ; iv(2) = i ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j+1
                elseif (ii.eq.4) then
                    iv(1) = i ; jv(1) = j ; iv(2) = i+1 ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j
                endif
                iv(4) = iv(1)
                jv(4) = jv(1)

                ! Find all angles using vector dot product a*b = |a||b|cos(a)
                do k = 2,3
                    xa = cord(jv(k+1),iv(k+1),1)-cord(jv(k),iv(k),1)
                    ya = cord(jv(k+1),iv(k+1),2)-cord(jv(k),iv(k),2)
                    xxal = sqrt(xa*xa + ya*ya) 
                    xb = cord(jv(k-1),iv(k-1),1)-cord(jv(k),iv(k),1)
                    yb = cord(jv(k-1),iv(k-1),2)-cord(jv(k),iv(k),2)
                    xxbl = sqrt(xb*xb + yb*yb) 

                    angle(k) = raddeg*acos((xa*xb+ya*yb)/(xxal*xxbl))
                end do
                angle (1) = 180.d0-angle(2)-angle(3)

                ! min angle in one trianle
                anglemin1 = min(angle(1),angle(2),angle(3))

                ! min angle in the whole mesh
                anglemint = min(anglemint, anglemin1)

            end do

        end do
    end do

    if( dtout_screen .ne. 0 ) then
        !$ACC wait(1)
        write (6,'(A,F6.2,A,F10.6)') '        min.angle=',anglemint, '     dt(yr)=',dt/sec_year
        write (333,'(A,F6.2,A,F10.6)') '        min.angle=',anglemint, '     dt(yr)=',dt/sec_year
        flush (333)
    endif
    ! check if the angle is smaller than angle of remeshing  
    if (anglemint .le. angle_rem) then
        if( dtout_screen .ne. 0 ) then
            print *, 'Remeshing due to angle required.'
            write(333,*) 'Remeshing due to angle required.'
        else
            call SysMsg('TEST_MESH: Remeshing due to angle required.')
        endif
        need_remeshing = 1
    endif
  endif

  return
end subroutine itest_mesh
