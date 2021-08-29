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
use icepack_parameters, only: p01, p5, c0, c1, c2, c3, c4, c10
use icepack_parameters, only: bignum, puny, gravit, pi
use icepack_tracers, only: nt_fsd
use icepack_warnings, only: warnstr, icepack_warnings_add,  icepack_warnings_aborted
use icepack_fsd
!use ice_domain_size, only: ncat, nfsd, nfreq
!use ice_state ! noah day wim 004
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
subroutine init_floe(trcrn)
!
! !USES:
!
!  use ice_state, only: nt_fsd, trcrn
!  use ice_flux,  only: floe_rad_c
!
!EOP
!
!  integer (kind=int_kind)            :: i, j


real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
   trcrn           ! tracer array
! noah day wim 005 -------

write(warnstr,*) 'INIT FLOE'
call icepack_warnings_add(warnstr)

! --------------------

!  if (trim(runtype) == 'continue') restart_fsd = .true.



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
  subroutine init_floe_0(trcrn)
!
! !USES:
!
!  use ice_state, only: nt_fsd, trcrn
!  use ice_kinds_mod,  only: floe_rad_c
!
!EOP
!
!  integer (kind=int_kind)            :: i, j
real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
   trcrn           ! tracer array


!trcrn(:,:,nt_fsd,:,:) = c0
!floe_rad_c(:,:,:) = c0
!do j = 1, ny_block
!do i = 1, nx_block
! trcrn(i,j,nt_fsd,:,:) = c300
! floe_rad_c(i,j,:)              = c300
!enddo
!enddo

! noah day wim 006 -------

write(warnstr,*) 'INIT FLOE', trcrn(nt_fsd,:)
call icepack_warnings_add(warnstr)

! --------------------

end subroutine init_floe_0

!=======================================================================

!BOP
!
! !ROUTINE: increment_floe
!
! !DESCRIPTION:
!
!  Increase ice floe tracer by scaled timestep length.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine increment_floe (nx_block, ny_block,                   &
                                dt, tmask, Lcell,                      &
                                loc_swh, loc_ppd, loc_mwd,             &
                                ifloe,  afice,  vfice, dum_wavemask )
!
! !USES:
!
      !use m_prams_waveice
      !use m_waveice
      !use m_waveattn

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block       ! block dimensions

!    ! time step
      real (kind=dbl_kind), intent(in) :: &
         dt

!     ! ice concentration & thickness (1st category)
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         afice, vfice

!     ! ice floe size param, significant wave height, wave peak period
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         ifloe , loc_swh, loc_ppd, loc_mwd

!     ! ice floe size param, significant wave height, wave peak period
!      real (kind=dbl_kind), dimension(nx_block,ny_block) ::          loc_mwd

!     ! land/boundary mask, thickness (T-cell)
      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          Lcell

	  integer (kind=int_kind), intent(in)   :: dum_wavemask
!
!  local variables
!
!     ! Counters
      integer (kind=int_kind)            :: i, j, ij

!	  ! For the wave-ice code:
!     ! counter
      integer :: lp_i, lp_j
!     ! points in freq & angular domains
      integer, parameter          :: nw_in=31, nth_in=1
!     ! max floe size param
      real(kind=dbl_kind)                :: D1
!     ! cell length
!      real(kind=dbl_kind), parameter     :: Lcell = 25*1000d0
!     ! open ocean wave params
      real(kind=dbl_kind), parameter                  :: fmin = 1d0/16d0, fmax = 1d0/6d0
                                                  ! freq min/max
      !real(kind=dbl_kind), parameter                  :: om1=2*pi*fmin, om2=2*pi*fmax
                                                  ! ang freqs
      !real(kind=dbl_kind), parameter                  :: om_0 = (om2 - om1)/(nw_in-1)
      !real(kind=dbl_kind), dimension(nw_in)           :: om_in, T
                                                  ! freqency
      real(kind=dbl_kind), dimension(nth_in)          :: th_in
                                                  ! direction (not used yet!!!!)
 	  real(kind=dbl_kind), dimension(nw_in)           :: S_init_in
 	                                              ! initial energy spec

 	  real(kind=dbl_kind), dimension(nw_in)           :: lam_wtr_in, k_wtr_in
 	                                              ! wavelen/num

 	  real(kind=dbl_kind)                             :: dum_sm0, dum_sm2
 	  											  ! dummy spectral moments

 	  real(kind=dbl_kind), parameter                  :: ws_tol = 1.0e-1_dbl_kind

!     ! spectrum passed through ice-covered ocean
      real (kind=dbl_kind), dimension(nx_block,nw_in) :: &
         wspec_row

!     ! dummy version of wspec_row to help with directionality
      real (kind=dbl_kind), dimension(nx_block,nw_in) :: &
         wspec_row_hld

!     ! sine mean wave directions for row
      real (kind=dbl_kind), dimension(nx_block) :: &
         mwd_row

!     ! array to hold values for calculation of mean wave directions
      real (kind=dbl_kind), dimension(2,nx_block) :: &
         mwd_hld

      integer, dimension(nx_block)                    :: tmt, tmt_hld ! WIM termination flag

!     ! attenuation
	  integer              :: idd_alp
      real(kind=dbl_kind)  :: dum_alp
!      character(30) :: fname_alp
!
!EOP
!


!!! Calculate wavenumbers & wavelengths

	 ! do lp_i=1,nw_in
    !   om_in(lp_i)        = om1 + (lp_i-1)*om_0
    !   T(lp_i)            = 2d0*pi/om_in(lp_i)
    !   lam_wtr_in(lp_i)   = gravity*(T(lp_i)**2d0)/2d0/pi
    !   k_wtr_in(lp_i)     = 2d0*pi/lam_wtr_in(lp_i)
    !end do


    ! noah day wim 007 -------

    !write(warnstr,*) 'INIT FLOE', lam_wtr_in(lp_i)
    !call icepack_warnings_add(warnstr)

    ! --------------------


      end subroutine increment_floe

!=======================================================================

function fn_Attn_MBK(dum_om)
!
! !USES:
!
! !INPUT PARAMETERS:
!

!real(kind=8), intent (in) :: dum_om        ! ang freq
real(kind=dbl_kind), intent (in), dimension(25) :: dum_om        ! ang freq

!
! !OUTPUT PARAMETERS
!

real(kind=dbl_kind), dimension(25) :: fn_Attn_MBK

!
!EOP
!

real(kind=8), parameter :: beta0 = 5.376168295200780E-005, &
   beta1 = 2.947870279251530E-005

fn_Attn_MBK = beta0*(dum_om**2) + beta1*(dum_om**4)

!fn_Attn_MBK = attn_fac*fn_Attn_MBK

end function fn_Attn_MBK


!=======================================================================

end module icepack_floe

!=======================================================================
