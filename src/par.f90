! -*- F90 -*-

program DREZINA
!$ACC routine(marker2elem) gang
use arrays
use params
use marker_data

character*200 inputfile
real*4 secnds,time0
integer :: narg, iargc, j, irestart, u_out
double precision :: dtacc_file, dtacc_save, dtacc_screen
logical :: file_exists

narg = iargc()
if(narg /= 1) then
    write(*,*) 'usage: flac inputfile'
    stop 1
endif
call getarg(1, inputfile)

! unit 333 is also written to directly (without an open) from
! rem_test.f90, init_marker.f90, and fl_move.f90; connecting it here
! to output.asc means those writes land there too, instead of the
! default fort.333 fallback.
u_out = 333
open( unit=u_out, file='output.asc' )

time0 = secnds(0.0)

! Read task parameters
call read_params(inputfile)
call allocate_arrays(nz, nx, nphase)
call allocate_markers(nz, nx)

! Try to read save-file contents. If file exist - restart, othewise - new start
inquire(file='_contents.rs', exist=file_exists)
if (file_exists) then
    irestart = 1
else
    irestart = 0
endif

if ( irestart .eq. 1 ) then  !file exists - restart
    call rsflac
    if( dtout_screen .ne. 0 ) then
        write(6,*) 'you CONTINUE from  ', nloop, ' step'
    else
        call SysMsg('you CONTINUE the execution')
    endif
else ! file does not exist - new start
    if( dtout_screen .ne. 0 ) then
        write(6,*) 'you have NEW start conditions'
        write(u_out,*) 'you have NEW start conditions'
    else
        call SysMsg('you have NEW start conditions')
    endif
    call setflac
    ! Output of initial configuration
    call outflac
    call outmarker
end if


!      ****************** running ********************************
dtacc_screen = 0
dtacc_file = 0
dtacc_save = 0

do while( time .le. time_max )
  if( dtout_screen .ne. 0 ) then
    if( dtacc_screen .gt. dtout_screen ) then
       write(*,'(I10,A,F7.3,A,F8.1,A)') nloop,'''s step. Time[My]=', time/sec_year/1.d+6, &
                ',  elapsed sec-', secnds(time0)
       write(u_out,'(I10,A,F7.3,A,F8.1,A)') nloop,'''s step. Time[My]=', time/sec_year/1.d+6, &
                ',  elapsed sec-', secnds(time0)
       dtacc_screen = 0
    endif
  endif

  do j = 1, ntest_rem
    ! FLAC
    call flac

    nloop = nloop + 1
    time = time + dt
    !$ACC update device(time,nloop) async(1)
    dtacc_screen = dtacc_screen + dt
    dtacc_file   = dtacc_file   + dt
    dtacc_save   = dtacc_save   + dt
  end do

  ! Remeshing
  if( ny_rem.eq.1 .and. itherm.ne.2 ) then
    call itest_mesh(need_remeshing)
    if( need_remeshing .ne. 0 ) then
      ! If there are markers recalculate their x,y global coordinate and assign them aps, eII, press, temp
      call bar2euler
      call remesh
      ! If markers are present recalculate a1,a2 local coordinates and assign elements with phase ratio vector
      call lpeuler2bar
      !$ACC parallel async(1)
      call marker2elem
      !$ACC end parallel
      !$ACC update self(nmarkers) async(1)
      !$ACC kernels async(1)
      dhacc_correction(:) = 0.d0
      !$ACC end kernels
    endif
  endif

    ! OUTPUT  
  !$ACC wait
  if( dtout_file .ne. 0 ) then 
    if( dtacc_file .gt. dtout_file ) then
      call outflac
      call outmarker
      dtacc_file = 0
    endif
  endif
  ! SAVING the restart information
  if( dtsave_file .ne. 0 ) then 
    if( dtacc_save .gt. dtsave_file ) then
      call saveflac
      dtacc_save = 0
    endif
  endif

end do

! SAVING the restart information of last step
call saveflac

close(u_out)

call SysMsg('Congratulations !')
end program DREZINA
