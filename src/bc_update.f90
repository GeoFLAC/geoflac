!--------------------------------------------------------------------
!
! Update and apply the boundary conditins in STRESSES 
!
!--------------------------------------------------------------------
subroutine bc_update
  !----------------------- determ. boundary conditions -----------
  !$ACC routine(Eff_dens) seq
  use arrays
  use params
  implicit none
  integer :: i, j, iph, ii1, ii2, jj1, jj2
  double precision :: rogh, densT, dh1, dh2, dh, dPT, dP, press, dlx, dly, press_norm, s_normal, s_shear, Eff_dens
  ! -------------------------------------------------------------
  !      LEFT BOUNDARY
  !   Update for hydrostatic force (normal to the surface)
  !   Shear component = 0 (perfect fluid)
  ! -------------------------------------------------------------
  !$ACC kernels async(1)
  force = 0.0d0
  !$ACC end kernels

  if (nydrsides.eq. 1) then
      !$ACC serial async(1)
      rogh = 0.d0
      do j=1,nz-1
          iph = iphase(j,1)
          densT = Eff_dens(j,1)
          dh1 = cord (j,1,2) - cord (j+1,1,2)
          dh2 = cord (j,2,2) - cord (j+1,2,2)
          dh  = 0.5d0 * (dh1+dh2)
          dPT = densT * g * dh
          dP = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )
          press = rogh +0.5d0*dP

          dlx = cord(j+1,1,1)-cord(j,1,1)
          dly = cord(j+1,1,2)-cord(j,1,2)
          press_norm = -press
          ! Normal component
          !           projection on x
          force(j,1,1) = force(j,1,1)+0.5d0*press_norm*dly
          force(j+1,1,1) = force(j+1,1,1)+0.5d0*press_norm*dly
          !        write(*,*) j,force(j,1,1),force(j+1,1,1)

          !           projection on y
          force(j,1,2) = force(j,1,2)-0.5d0*press_norm*dlx
          force(j+1,1,2) = force(j+1,1,2)-0.5d0*press_norm*dlx

          rogh = rogh +dP
      enddo

      ! -------------------------------------------------------------
      !      RIGHT BOUNDARY
      !   Update for hydrostatic force (normal to the surface)
      !   Shear component = 0 (perfect fluid)
      ! -------------------------------------------------------------
      rogh = 0.
      do j=1,nz-1
          iph = iphase(j,nx-1)
          densT = Eff_dens(j,nx-1)
          dh1 = cord (j,nx-1,2) - cord (j+1,nx-1,2)
          dh2 = cord (j,nx,2) - cord (j+1,nx,2)
          dh  = 0.5d0 * (dh1+dh2)
          dPT = densT * g * dh
          dP = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )
          press = rogh +0.5d0*dP

          !        write(*,*) j,stress0(j,nx-1,1,1),press_norm,densT
          dlx = cord(j+1,nx,1)-cord(j,nx,1)
          dly = cord(j+1,nx,2)-cord(j,nx,2)
          press_norm = -press
          ! Normal component
          !
          !           projection on x
          force(j,nx,1) = force(j,nx,1)-0.5d0*press_norm*dly
          force(j+1,nx,1) = force(j+1,nx,1)-0.5d0*press_norm*dly

          !           projection on y
          force(j,nx,2) = force(j,nx,2)+0.5d0*press_norm*dlx
          force(j+1,nx,2) = force(j+1,nx,2)+0.5d0*press_norm*dlx

          rogh = rogh +dP
      enddo
      !$ACC end serial

  endif

  !----------------------------------
  ! APPLY  STATIC STRESSES:
  !----------------------------------
  !$ACC parallel loop async(1)
  do i=1,nopbmax

      ii1 = nopbou(i,1)
      ii2 = nopbou(i,2)
      jj1 = nopbou(i,3)
      jj2 = nopbou(i,4)
      dlx = cord(jj2,ii2,1)-cord(jj1,ii1,1)
      dly = cord(jj2,ii2,2)-cord(jj1,ii1,2)

      !
      ! Normal component
      !

      s_normal = bcstress(i,1)
      !           projection on x
      !$ACC atomic update
      force(jj1,ii1,1) = force(jj1,ii1,1)-0.5d0*s_normal*dly
      !$ACC atomic update
      force(jj2,ii2,1) = force(jj2,ii2,1)-0.5d0*s_normal*dly
      !           projection on z
      !$ACC atomic update
      force(jj1,ii1,2) = force(jj1,ii1,2)+0.5d0*s_normal*dlx
      !$ACC atomic update
      force(jj2,ii2,2) = force(jj2,ii2,2)+0.5d0*s_normal*dlx
      !       write(*,*) jj1,jj2,ii1,ii2,s_normal,force(jj1,ii1,1), force(jj2,ii2,1)
      !
      ! Shear component
      !

      s_shear  = bcstress(i,2)
      !           projection on x
      !$ACC atomic update
      force(jj1,ii1,1) = force(jj1,ii1,1)+0.5d0*s_shear*dlx
      !$ACC atomic update
      force(jj2,ii2,1) = force(jj2,ii2,1)+0.5d0*s_shear*dlx
      !           projection on z
      !$ACC atomic update
      force(jj1,ii1,2) = force(jj1,ii1,2)+0.5d0*s_shear*dly
      !$ACC atomic update
      force(jj2,ii2,2) = force(jj2,ii2,2)+0.5d0*s_shear*dly


      !       if (i_sl.eq.1) then
      !       s_shearo  = bcstress(i,3)
      !           projection on x
      !       force(nn1,3) = force(nn1,3)+0.5d0*s_shearo*dlx
      !       force(nn2,3) = force(nn2,3)+0.5d0*s_shearo*dlx
      !       write(*,*) bcstress(i,3)
      !       endif
  enddo
  return
end subroutine bc_update
