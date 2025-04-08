subroutine marker2elem
  !$ACC routine gang
  use myrandom_mod
  use marker_data
  use arrays
  use params
  use phases
  implicit none

  integer :: i, j, kinc, inc, iseed, icount
  double precision :: x1, x2, y1, y2, xx, yy, r1, r2
  
  !character*200 msg

  ! Interpolate marker properties into elements
  ! Find the triangle in which each marker belongs

  !$OMP parallel do private(kinc,r1,r2,x1,y1,x2,y2,xx,yy,inc,icount)
  !$ACC loop collapse(2) gang vector
  do i = 1 , nx-1
      do j = 1 , nz-1
          iseed = nloop + i + j
          kinc = nmark_elem(j,i)
          icount = 0

          !  if there are too few markers in the element, create a new one
          !  with age 0 (similar to initial marker)
          !if(kinc.le.4) then
          !    write(msg,*) 'marker2elem: , create a new marker in the element (i,j))', i, j
          !    call SysMsg(msg)
          !endif

          do while (kinc.le.4)
              call myrandom(iseed, r1)
              call myrandom(iseed, r2)

              ! (x1, y1) and (x2, y2)
              x1 = cord(j  ,i,1)*(1-r1) + cord(j  ,i+1,1)*r1
              y1 = cord(j  ,i,2)*(1-r1) + cord(j  ,i+1,2)*r1
              x2 = cord(j+1,i,1)*(1-r1) + cord(j+1,i+1,1)*r1
              y2 = cord(j+1,i,2)*(1-r1) + cord(j+1,i+1,2)*r1

              ! connect
              ! (this point is not uniformly distributed within the element area
              ! and is biased against the thicker side of the element, but this
              ! point is almost gauranteed to be inside the element)
              xx = x1*(1-r2) + x2*r2
              yy = y1*(1-r2) + y2*r2

              call add_marker(xx, yy, iphase(j,i), zpressm(j,i), Eff_melt(j,i), 0.d0, j, i, inc)
              icount = icount + 1
              if(icount > 100) stop 133
              if(inc.le.0) cycle

              kinc = kinc + 1
              zpresscounter(j,i) = zpresscounter(j,i) + zpressm(j,i)
              Emeltcounter(j,i) = Emeltcounter(j,i) + Eff_melt(j,i)
          enddo

          Eff_melt(j,i) = Emeltcounter(j,i)/float(kinc)
          zpressm(j,i) = zpresscounter(j,i)/float(kinc)

          call count_phase_ratio(j,i)

      enddo
  enddo
  !$OMP end parallel do

! Find the Moho
if (itype_melting.eq.2) then
    !$OMP parallel do private(i,j)
    !$ACC loop auto
    do i = 1, nx-1
        jmoho(i) = nz-1
        do j = 1, nz-1
            if (sum(phase_ratio(mantle_phases,j,i)) > 0.5d0) then
                if (phase_ratio(kocean0,j,i) > 0.1.or.fmelt(j,i)>0.) then
                    jmoho(i) = j - 1
                else
                    jmoho(i) = j
                endif
                exit
            endif
        enddo
        if (jmoho(i).le.8) jmoho(i) = 8
    enddo
    !$OMP end parallel do
else
    !$OMP parallel do private(i,j)
    !$ACC loop auto
    do i = 1, nx-1
        jmoho(i) = nz-1
        do j = 1, nz-1
            if (sum(phase_ratio(mantle_phases,j,i)) > 0.5d0) then
                jmoho(i) = j
                exit
            endif
        enddo
    enddo
    !$OMP end parallel do
endif

if (nloop .ne. 0 .and. itype_melting.eq.2) call check_chamber

return
end subroutine marker2elem
