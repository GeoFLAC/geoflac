! Reads problem parameters from input file

subroutine read_params(inputfile)
use arrays
use params
implicit none
character*200, intent(in) :: inputfile
integer :: i, j, iu, line, u11

iu = 4
open( iu, file=inputfile )
line = 1

! MESH
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) nx, nz
line = line + 1
open(newunit=u11, file='nxnz.0')
write(u11,*) nx, nz
close(u11)

nx = nx+1
nz = nz+1
nmarkers = 9*(nz-1)*(nx-1)
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) x0, z0
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) rxbo, rzbo
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) ircoord, coordfile
line = line + 1
if (ircoord .gt. 0) then
    ! ignore next two lines
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) nzonx
    line = line + 1
    if (nzonx .ne. 0) stop 'ircoord is set, but nzonx is not 0!'
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) nzony
    line = line + 1
    if (nzony .ne. 0) stop 'ircoord is set, but nzony is not 0!'
else
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) nzonx
    line = line + 1
    if (nzonx .eq. 0) then
        nzonx = 1
        nelz_x(1) = nx - 1
        sizez_x(1) = 1.d0
    else
        do i = 1, nzonx
            call AdvanceToNextInputLine(4, line)
            read(4,*,err=1000) nelz_x(i), sizez_x(i)
            line = line + 1
        end do
    endif

    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) nzony
    line = line + 1
    if (nzony .eq. 0) then
        nzony = 1
        nelz_y(1)    = nz - 1
        sizez_y(1)   = 1.d0
    else
        do i = 1,nzony
            sizez_y(i) = 1.d0
            call AdvanceToNextInputLine(4, line)
            read(4,*,err=1000) nelz_y(i), sizez_y(i)
            line = line + 1
        end do
    endif
endif

! MECHANICAL CONDITIONS
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) nystressbc,nydrsides
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) nofbc
line = line + 1
call AdvanceToNextInputLine(4, line)
        do i = 1,nofbc
        read(iu,*,err=1000) nofside(i),nbc1(i),nbc2(i),nbc(i),  &
        bca(i),bcb(i),bcc(i),  &
        bcd(i),bce(i),bcf(i),bcg(i),bch(i),bci(i)
        line = line + 1
        enddo
call AdvanceToNextInputLine(4, line)
! hydrostatic pressure applied at the bottom
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) nyhydro,pisos,iphsub,drosub,damp_vis
line = line + 1
! gravity
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) g
line = line + 1



! THERMAL CONDITIONS
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) i_prestress
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) extra_pres
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) itherm
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) istress_therm
line = line + 1
call AdvanceToNextInputLine(4, line)             ! thermal stresses
read (4,*,err=1000) ishearh                           ! shear heating
line = line + 1
call AdvanceToNextInputLine(4, line)
read (4,*,err=1000) t_top
line = line + 1
call AdvanceToNextInputLine(4, line)
read (4,*,err=1000) t_bot
line = line + 1
call AdvanceToNextInputLine(4, line)
read (4,*,err=1000) hs, hr
line = line + 1
! boundary conditions at the bottom (1-T,2-Flux) 
call AdvanceToNextInputLine(4, line)
read (4,*,err=1000) itemp_bc, bot_bc
line = line + 1
if( itemp_bc.eq.2 ) bot_bc = bot_bc/1000  ! convert in W/m3
! Predefined distributions
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) irtemp
line = line + 1
if ( irtemp .gt. 0 ) then
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) tempfile
    line = line + 1
else
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000)
    line = line + 1
endif
! temp structure
call AdvanceToNextInputLine(4, line)
read (4,*,err=1000) nzone_age
line = line + 1
call AdvanceToNextInputLine(4, line)
do i = 1, nzone_age
    read (4,*,err=1000) ictherm(i),age_1(i),tp1(i),tp2(i),ixtb1(i),ixtb2(i)
    line = line + 1
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) nph_layer(i), (hc(i,j), j=1,nph_layer(i)-1)
    line = line + 1
    !print *, (hc(i,j), j=1,nph_layer(i)-1)
    read(4,*,err=1000) (iph_col(i,j), j=1,nph_layer(i))
    line = line + 1
    !print *, (iph_col(i,j), j=1,nph_layer(i))
enddo
! check smooth nzone_age
do i = 1, nzone_age-1
    iph_col_trans(i) = 0
    if (ixtb2(i) == -1) then
        if (ixtb1(i+1) /= -1) then
            print *, 'Error: ixtb1 is not -1 at', i+1, 'column!'
        endif
        if (ictherm(i) /= ictherm(i+1)) then
            print *, 'Error: ictherm at', i, i+1, 'columns are not equal!'
        endif
        if (nph_layer(i) /= nph_layer(i+1)) then
            print *, 'Error: nph_layer at', i, i+1, 'columns are not equal!'
        endif
        do j = 1, nph_layer(i)
            ! the phases in the columns must match
            if (iph_col(i,j) /= iph_col(i+1,j)) then
                print *, 'Error: iph_col at', i, i+1, 'columns are not equal!'
            endif
        enddo
        ! flag indicates smooth transition
        iph_col_trans(i) = 1
        ! remove '-1'
        ixtb2(i) = ixtb2(i+1)
        ixtb1(i+1) = ixtb1(i)
    endif
enddo

! RHEOLOGY
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) nphase
line = line + 1
do i = 1, nphase 
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) irheol(i),visc(i),den(i),alfa(i),beta(i),pln(i),acoef(i),eactiv(i),vactiv(i), &
         rl(i),rm(i),plstrain1(i),plstrain2(i),fric1(i),fric2(i),cohesion1(i),cohesion2(i), &
         dilat1(i),dilat2(i), &
         conduct(i),cp(i),ts(i),tl(i),tk(i),erodibility(i)
    line = line + 1
end do
print *, nphase, 'phases'
! Flag to take initial phase distribution from a file
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) irphase
line = line + 1
!write(*,*) irphase

if ( irphase .gt. 0 ) then
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) phasefile
    line = line + 1
else
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000)
    line = line + 1
endif

! inclusions
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) inhom
line = line + 1
if( inhom .gt. maxinh ) then
    call SysMsg('Read_params: Increase arrays for inhomogenities')
    stop 26
endif
do i = 1, inhom
    call AdvanceToNextInputLine(4, line)
    read(4,*,err=1000) ix1(i), ix2(i), iy1(i), iy2(i), inphase(i), igeom(i), xinitaps(i)
    line = line + 1
end do
! Tension cut-off
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) ten_off
line = line + 1
!linear healing parameter
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) tau_heal
line = line + 1
! viscosity limits
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) v_min, v_max, ivis_shape,efoldc
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) itype_melting, nelem_serp, nelem_dike, prod_magma, prod_magma2, rho_magma
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) angle_mzone, fmagma_max, ratio_mantle_mzone, ratio_crust_mzone
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) latent_heat_magma, lambda_freeze, lambda_freeze_tdep
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) weaken_ratio_plastic, weaken_ratio_viscous
line = line + 1
weaken_ratio_viscous = log(weaken_ratio_viscous)

! REMESHING
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  ny_rem, mode_rem, ntest_rem, angle_rem
line = line + 1
if ( mode_rem.ne.1 .and. mode_rem.ne.11 .and. mode_rem.ne.3 .and. mode_rem.ne.4) then
    call SysMsg('Illegal remeshing mode! Allowable - 1, 3, 4or 11')
    stop
endif
! dx_rem - remeshing criteria for mode_rem=11
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  dx_rem
line = line + 1
! diffusion of topography
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  topo_kappa, fac_kappa
line = line + 1


! PROCESS CONTROL
! inertial mass scaling
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  idt_scale
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  dt_scale, strain_inert
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  i_rey,xReyn
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  ifreq_rmasses
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  ifreq_imasses
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  ifreq_visc
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000)  ifreq_avgsr
line = line + 1

if (mod(ntest_rem, ifreq_rmasses) + mod(ntest_rem, ifreq_imasses) &
    + mod(ntest_rem, ifreq_visc) + mod(ntest_rem, ifreq_avgsr) .ne. 0) then
        call SysMsg('ntest_rem must be multiples of process frequency')
        stop
endif

call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) frac, fracm
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) movegrid,ndim
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) demf, mix_strain, mix_stress
line = line + 1

! OUTPUT PARAMETERS
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) time_max
line = line + 1
time_max = time_max * 1000 * sec_year
call AdvanceToNextInputLine(4, line)
read (4,*,err=1000) dtout_screen
line = line + 1
dtout_screen = dtout_screen * 1000 * sec_year
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) dtout_file
line = line + 1
dtout_file = dtout_file * 1000 * sec_year
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz,io_sxz,io_pres, &
    io_temp,io_phase,io_visc,io_unused,io_density,io_src,io_diss,io_forc,io_hfl,io_topo
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) lastout
line = line + 1
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) dtsave_file
line = line + 1
dtsave_file = dtsave_file * 1000 * sec_year
call AdvanceToNextInputLine(4, line)
read(4,*,err=1000) lastsave
line = line + 1

close (iu)

return
1000  print *, 'Error reading file "', trim(inputfile), '" at line', line
stop 11
end


subroutine AdvanceToNextInputLine(iu, line)
implicit none
integer, intent(in) :: iu
integer, intent(inout) :: line
character(len=1) :: buf
integer :: ios

do
    line = line + 1
    read(iu, '(A1)', iostat=ios) buf
    if (ios /= 0) then
        print *, 'AdvanceToNextInputLine: EOF or error reached!'
        stop 11
    endif
    if (buf(1:1) /= ';') then
        ! this line is not a comment, rewind
        backspace(iu)
        line = line - 1
        return
    endif
end do

end subroutine AdvanceToNextInputLine
