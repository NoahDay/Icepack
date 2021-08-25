!=======================================================================
!
!BOP
!
! !MODULE: icepack_floe - Floe size distribution - from the cicewithwaves
!                                                  project
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors Noah Day, using the work from
!         Luke Bennetts and Siobhan O'Farrell
! date: 25-08-21
!
! !INTERFACE:
!
      module icepack_floe
!
!
! !USES:
!
use icepack_kinds
use icepack_parameters, only: c0, c1, c2, c4, p01, p1, p5, puny
use icepack_parameters, only: pi, floeshape, wave_spec, bignum, gravit, rhoi
use icepack_tracers, only: nt_fsd, tr_fsd
use icepack_warnings, only: warnstr, icepack_warnings_add
use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted
!
!EOP
!
  implicit none

  !logical (kind=log_kind) :: &
     !tr_fsd,       & ! if .true., use floe tracer
     !restart_fsd,   & ! if .true., read floe tracer restart file
     !tr_wspec,       & ! if .true., use wave spec tracer
     !restart_spec      ! if .true., read wave spec tracer restart file

!=======================================================================

    contains

!=======================================================================
!BOP
!
! !ROUTINE: init_floe
!
! !DESCRIPTION:
!
!  Initialize ice floe tracer (call prior to reading restart data)
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  subroutine init_floe
if (tr_fsd) then
  write(warnstr,*) 'INIT FLOE', nt_fsd
  call icepack_warnings_add(warnstr)
end if
!
! !USES:
!
!  use ice_state, only: nt_fsd, trcrn
!  use ice_flux,  only: floe_rad_c
!
!EOP
!
!  integer (kind=int_kind)            :: i, j

!  if (trim(runtype) == 'continue') restart_fsd = .true.

!  if (restart_fsd) then
!     call read_restart_fsd
!  else
!     call init_floe_0
!         trcrn(i,j,nt_fsd,:,:) = c0
!         do j = 1, wavemask
!          do i = 1, nx_block
!           trcrn(i,j,nt_fsd,:,:) = c300
!          enddo
!         enddo
!         floe_rad_c(i,j,:) = c0
!         do j = 1, wavemask
!          do i = 1, nx_block
!           floe_rad_c(i,j,:) = c300
!          enddo
!         enddo
!  endif

  end subroutine init_floe

!=======================================================================
!BOP
!
! !ROUTINE: init_floe_0
!
! !DESCRIPTION:
!
!  Initialize ice floe tracer
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  subroutine init_floe_0
!
! !USES:
!
!  use ice_state, only: nt_fsd, trcrn
!  use ice_flux,  only: floe_rad_c
!
!EOP
!
!  integer (kind=int_kind)            :: i, j

!trcrn(:,:,nt_fsd,:,:) = c0
!floe_rad_c(:,:,:) = c0
!do j = 1, ny_block
!do i = 1, nx_block
! trcrn(i,j,nt_fsd,:,:) = c300
! floe_rad_c(i,j,:)              = c300
!enddo
!enddo

end subroutine init_floe_0






!=======================================================================

end module icepack_floe

!=======================================================================
