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
        use ice_domain_size, only: ncat, nfsd, nfreq
!use ice_state ! noah day wim 004
        use ice_fileunits, only: nu_diag
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
!
! !USES:
 use ice_arrays_column, only: floe_rad_c
 use ice_state, only: trcrn
 use icepack_tracers, only: nt_fsd, tr_fsd
!
!EOP
!
!  integer (kind=int_kind)            :: i, j


!real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
!   trcrn           ! tracer array
! noah day wim 005 -------

    call init_floe_0

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
  subroutine init_floe_0
!
! !USES:
!
    use ice_state, only: trcrn!, nt_ifloe
    use icepack_tracers, only: nt_fsd, nfsd
    use ice_flux,  only: ifd
    use ice_blocks, only: nx_block, ny_block
    use ice_constants, only: c0, c1, c2, c3, c4, c10, c30, c110, c300, c1000, p25, p01
!
!EOP
!
  integer (kind=int_kind)            :: i, j

trcrn(:,:,nt_fsd,:,:) = c0
ifd(:,:,:) = c0
do j = 1, ny_block
  do i = 1, nx_block
    trcrn(:,:,nt_fsd,1,:) = c300 ! Noah Day WIM, instead of :,: at the end
    ifd(i,j,:)            = c300
  enddo
enddo


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

! !USES:
!
      use m_prams_waveice
      use m_waveice
      use m_waveattn

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
          loc_swh, loc_ppd, loc_mwd, ifloe

!     ! ice floe size param, significant wave height, wave peak period

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
      real(kind=dbl_kind), parameter                  :: om1=2*pi*fmin, om2=2*pi*fmax
                                                  ! ang freqs
      real(kind=dbl_kind), parameter                  :: om_0 = (om2 - om1)/(nw_in-1)
      real(kind=dbl_kind), dimension(nw_in)           :: om_in, T
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
!
!EOP
!

!!! Define wave mask edge and propagate waves into ice:
     if (cmt.ne.0) then
      write(nu_diag,*) '    in increment_floe -> wavemask=', dum_wavemask
     endif

!!! Define attenuation coefficient coefficients

!	 if (ATTEN_MODEL.eq.1) then
!	  write(nu_diag,*) 'opening: ', waveicedatadir , fname_alp, '...'

!	  open(newunit=idd_alp,file=waveicedatadir // fname_alp)
!	  do lp_j=1,Ncheb_h+1
!	   do lp_i=1,Ncheb_f+1
!		read(idd_alp,*) dum_alp
	!	alp_coeffs(lp_i,lp_j)=dum_alp
	 !  end do
	 ! end do
	 ! close(idd_alp)
	 ! write(nu_diag,*) '... done'
	! endif



   !!! Calculate wavenumbers & wavelengths

   	  do lp_i=1,nw_in
          om_in(lp_i)        = om1 + (lp_i-1)*om_0
          T(lp_i)            = 2d0*pi/om_in(lp_i)
          lam_wtr_in(lp_i)   = gravity*(T(lp_i)**2d0)/2d0/pi
          k_wtr_in(lp_i)     = 2d0*pi/lam_wtr_in(lp_i)
         end do

         if (nth_in.ne.1) then
          do lp_i=1,nth_in
           th_in(lp_i)        = -pi/2 + (lp_i-1)*pi/(nth_in-1)
          end do
         else
          th_in = 0d0
         end if

    	 tmt(:)             = 0      ! flag to terminate routine
    	 tmt_hld(:)         = 1      ! flag to terminate routine
       wspec_row(:,:)     = c0     ! a dummy vector
       wspec_row_hld(:,:) = c0     ! a dummy dummy vector
       mwd_hld(:,:)       = c0     ! another dummy dummy vector
       !loc_mwd(:,:)       = c0     ! another dummy dummy vector
   	   S_init_in(:)       = c0

       !!! Begin at wavemask (only difference is initialisation)
       if (cmt.ne.0) then
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
        write(nu_diag,*) '                      -> starting wavemask loop'
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
       endif

       j = dum_wavemask

       ! A0. initialise with: Bretschneider
       do i=1,nx_block
        if (tmask(i,j)) then
        ! initialise direction
        mwd_row(i) = loc_mwd(i,j)
        if (loc_swh(i,j).lt.ws_tol) then
         if (cmt.ne.0) then
          write(nu_diag,*) '>>>--------------------------------------------->>>'
          write(nu_diag,*) '                         i, j = ', i, j
          write(nu_diag,*) '                         NOT running wave-ice routine'
          write(nu_diag,*) '                      -> sig wave ht <', ws_tol
          !write(nu_diag,*) '<<<---------------------------------------------<<<'
         endif
         tmt(i) = 1
         if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
          ifloe(i,j) = floe_sz_pancake
          if (cmt.ne.0) then
           if (vfice(i,j).lt.tolh) then
            write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
           endif
           if (afice(i,j).lt.tolice) then
            write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
           endif
           write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
          endif ! END IF COMMENT
         endif  ! END IF h<tolh OR c<tolc
        else ! IF there are waves:
         do lp_i=1,nw_in
          S_init_in(lp_i) = SDF_Bretschneider(om_in(lp_i),0,loc_swh(i,j),loc_ppd(i,j))
         end do

         ! want consistency between SWH definitions:
         !if (i.eq.1.and.cmt.ne.0) then

         if (cmt.ne.0) then
            write(nu_diag,*) '                      -> check: swh ', loc_swh(i,j)
            write(nu_diag,*) '                      ->        ppd ', loc_ppd(i,j)
            dum_sm0        = fn_SpecMoment(S_init_in,nw_in,nth_in,om_in,th_in,0,nu_diag)
   	     dum_sm2        = fn_SpecMoment(S_init_in,nw_in,nth_in,om_in,th_in,2,nu_diag)
            write(nu_diag,*) '                      ->        swh ', 4d0*(dum_sm0**0.5d0)
            write(nu_diag,*) '                      ->        ppd ', &
                 											2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
        endif ! cmt
         !endif ! END COMMENT

          ! A1. Use WIM to update wave spectrum and floe sizes

          if (cmt.ne.0) then
           write(nu_diag,*) '>>>--------------------------------------------->>>'
           write(nu_diag,*) '        INPUT         -> i,j   =', i, j
 		  write(nu_diag,*) '                      -> ifloe =', ifloe(i,j)
           write(nu_diag,*) '                      -> swh   =', loc_swh(i,j)
           write(nu_diag,*) '                      -> ppd   =', loc_ppd(i,j)
           write(nu_diag,*) '                      -> mwd   =', 180d0*mwd_row(i)/pi
    	     endif
    	     if (do_coupled.ne.0) then
    	      !call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j), &
 		  !nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag)
 		 else
 		  call sub_Uncoupled(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j), &
 		  nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag)
      !call sub_Uncoupled(Length_cell, hice_init, aice_in, &
      !  nfreq,nthh,dwavefreq,wave_spectrum_in,wave_spectrum_out,tmtt,nu_diag)
 		 endif
 		 ifloe(i,j) = D1
 		endif ! ENDIF ws_tol
        else
         if (cmt.ne.0) then
          write(nu_diag,*) '>>>--------------------------------------------->>>'
          write(nu_diag,*) '                         i, j = ', i, j
          write(nu_diag,*) '                      -> LAND'
          write(nu_diag,*) '<<<---------------------------------------------<<<'
         endif
         tmt(i) = 1
        endif ! ENDIF tmask(i,j)
       enddo ! ENDDO i=1,nx_block

       if (cmt.ne.0) then
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
        write(nu_diag,*) '                      -> do wave directions', dum_wavemask
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
       endif

       !print*, 'doing wave directions on wavemask loop'

       ! A2. Redistribute energy according to mean directions
        ! LH Boundary:
         i=1
         ! if the cell isn't land ...
       if (tmask(i,j).and.tmt(i).ne.1) then
        ! Calculate a dummy swh to weight directions:
        dum_sm0        = &
           fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
 	   dum_sm0        = 4d0*(dum_sm0**0.5d0)
        if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
         !write(nu_diag,*) '                      -> going south', 1
         !print*, '                      -> going south', 1
         wspec_row(1,:)  = wspec_row(1,:) + &
        						(1d0-sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
         mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(1)*dum_sm0
         mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
         tmt_hld(1)      = 0
        endif ! ENDIF -pi/4<mwd<pi/4
        ! if wave energy needs to be advected to the east ...
        if (tmask(2,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
         !write(nu_diag,*) '                      -> going east', 1
         wspec_row(2,:) = wspec_row(2,:) + &
         				(sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
         mwd_hld(1,2)   = mwd_hld(1,2) + mwd_row(1)*dum_sm0
         mwd_hld(2,2)   = mwd_hld(2,2) + dum_sm0
         tmt_hld(2)     = 0
         ! if wave energy needs to be advected to the west (wrap around vector) ...
        elseif (tmask(nx_block,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
         !write(nu_diag,*) '                      -> going west', 1
         wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                                (sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
         mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + mwd_row(1)*dum_sm0
         mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
         tmt_hld(nx_block)     = 0
        endif
 	  endif ! END IF TMASK

! ======

  ! RH Boundary:
    i=nx_block
      ! if the cell isn't land ...
    if (tmask(i,j).and.tmt(i).ne.1) then
     dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
   dum_sm0        = 4d0*(dum_sm0**0.5d0)
     if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
      !write(nu_diag,*) '                      -> going south', nx_block
      wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                            (1d0-sin(2d0*mwd_row(nx_block))**2d0)*wspec_row_hld(nx_block,:)
      mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + &
                             mwd_row(nx_block)*dum_sm0
      mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
      tmt_hld(nx_block)     = 0
     endif   ! ENDIF -pi/4<mwd<pi/4
      ! if wave energy needs to be advected to the WEST ...
     if (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
      !write(nu_diag,*) '                         going west', nx_block
      wspec_row(i-1,:) = wspec_row(i-1,:) + &
     	                         (sin(2d0*mwd_row(i))**2)*wspec_row_hld(i,:)
      mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
      mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
      tmt_hld(i-1)     = 0
      ! if wave energy needs to be advected to the EAST ...
     elseif (tmask(1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
      !write(nu_diag,*) '                         going east', nx_block
      wspec_row(1,:)  = wspec_row(1,:) + &
                             (sin(2d0*mwd_row(nx_block))**2)*wspec_row(nx_block,:)
      mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(nx_block)*dum_sm0
      mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
      tmt_hld(1)      = 0
     endif ! ENDIF (tmask(nx_block-1,j).and.sinmwd_row(nx_block).gt.c0)
  endif ! END IF TMASK

    ! loop the inner cells...
    do i=2,nx_block-1
     if (tmask(i,j).and.tmt(i).ne.1) then
      dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
    dum_sm0        = 4d0*(dum_sm0**0.5d0)
      if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       ! South
       !write(nu_diag,*) '                      -> going south', i
       wspec_row(i,:) = wspec_row(i,:) + &
       			(1d0-sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
       mwd_hld(1,i)   = mwd_hld(1,i) + mwd_row(i)*dum_sm0
       mwd_hld(2,i)   = mwd_hld(2,i) + dum_sm0
       tmt_hld(i)     = 0
      endif ! ENDIF -pi/4<mwd<pi/4
      ! East
      if (tmask(i+1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
       !write(nu_diag,*) '                      -> going east', i
       wspec_row(i+1,:) = wspec_row(i+1,:) + &
       					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
       mwd_hld(1,i+1)   = mwd_hld(1,i+1) + mwd_row(i)*dum_sm0
       mwd_hld(2,i+1)   = mwd_hld(2,i+1) + dum_sm0
       tmt_hld(i+1)     = 0
      elseif (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
       !write(nu_diag,*) '                      -> going west', i
       wspec_row(i-1,:) = wspec_row(i-1,:) + &
       				(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
       mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
       mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
       tmt_hld(i-1)     = 0
      end if
     endif ! END IF TMASK
  enddo

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> update mean values', dum_wavemask
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
    endif

    ! A3. Calculate mean parameters
    do i=1,nx_block

     if (tmask(i,j)) then
      if (mwd_hld(2,i).gt.c0) then
        mwd_row(i) = mwd_hld(1,i)/mwd_hld(2,i)
      else
        mwd_row(i) = c0
      endif
      tmt(i)         = tmt_hld(i)
      loc_mwd(i,j)   = mwd_row(i)
      if (tmt(i).eq.1) then
       loc_swh(i,j)   = 0d0
       loc_ppd(i,j)   = 0d0
      else
	 dum_sm0        = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
	 dum_sm2        = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,2,nu_diag)
	 loc_swh(i,j)   = 4d0*(dum_sm0**0.5d0)
	 loc_ppd(i,j)   = 2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
	endif
	if (cmt.ne.0) then
       write(nu_diag,*) '           OUTPUT     -> i,j   =', i, j
	 write(nu_diag,*) '                      -> tmt      =', tmt(i)
     write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
       write(nu_diag,*) '                      -> swh      =', loc_swh(i,j)
       write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j)
       write(nu_diag,*) '                      -> mwd      =', 180d0*loc_mwd(i,j)/pi
       write(nu_diag,*) '<<<---------------------------------------------<<<'
      endif
     endif ! ENDIF tmask(i,j)
    enddo ! ENDDO i=1,nblock

    wspec_row_hld(:,:) = c0     ! reset the dummy dummy vector
    mwd_hld(:,:)       = c0     ! ditto
    tmt_hld(:)         = 1

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> starting wave propagation'
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
    endif

    do j=dum_wavemask-1,2,-1
     if (cmt.ne.0) then
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
      write(nu_diag,*) '                      -> j=', j
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
     endif

     ! B1. Use WIM to update wave spectrum and floe sizes
     do i=1,nx_block
      !print*, 'Hi: i,j=', i, j
      if (tmt(i).eq.0) then
      !print*, '... tmt=0'
      if (tmask(i,j)) then
       if (loc_swh(i,j+1).lt.ws_tol) then
        if (cmt.ne.0) then
         write(nu_diag,*) '>>>--------------------------------------------->>>'
         write(nu_diag,*) '                         i, j = ', i, j
         write(nu_diag,*) '                         NOT running wave-ice routine'
         write(nu_diag,*) '                      -> sig wave ht <', ws_tol
        endif
        tmt(i) = 1
        if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
         ifloe(i,j) = floe_sz_pancake
         if (cmt.ne.0) then
          if (vfice(i,j).lt.tolh) then
           write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
          endif
          if (afice(i,j).lt.tolice) then
           write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
          endif
          write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
         endif ! END IF COMMENT
        endif  ! END IF h<tolh OR c<tolc
       else ! IF there are waves:
       !print*, '... tmask=true'
       if (cmt.ne.0) then
        write(nu_diag,*) '>>>--------------------------------------------->>>'
        write(nu_diag,*) '       INPUT          -> i,j      =', i, j
        !print*, '       INPUT          -> i,j      =', i, j
	  write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
	  !print*, '                      -> ifloe    =', ifloe(i,j)
        write(nu_diag,*) '                      -> swh      =', loc_swh(i,j+1)
        !print*, '                      -> swh      =', loc_swh(i,j+1)
        write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j+1)
        !print*, '                      -> ppd      =', loc_ppd(i,j+1)
        write(nu_diag,*) '                      -> mwd      =', 180d0*mwd_row(i)/pi
        !print*, '                      -> mwd      =', 180d0*mwd_row(i)/pi
       endif  ! ENDIF if (cmt.ne.0) then
       if (do_coupled.ne.0) then
         !call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j),nw_in,&
          !nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag)
       else
         call sub_Uncoupled(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j),nw_in,&
          nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag)
          !call sub_Uncoupled(Length_cell, hice_init, aice_in, &
          !  nfreq,nthh,dwavefreq,wave_spectrum_in,wave_spectrum_out,tmtt,nu_diag)
       endif ! ENDIF (do_coupled.ne.0)
       !print*, '... done wave-ice routine'
       ifloe(i,j) = D1
        !print*, '... set ifloe(i,j)=', ifloe(i,j), 'D1=', D1
       endif ! ENDIF ws_tol
       else
        if (cmt.ne.0) then
         write(nu_diag,*) '>>>--------------------------------------------->>>'
         write(nu_diag,*) '                         i, j = ', i, j
         write(nu_diag,*) '                      -> LAND'
         write(nu_diag,*) '<<<---------------------------------------------<<<'
        endif
       tmt(i)         = 1
      endif ! ENDIF tmask(i,j)
      else
       if (cmt.ne.0) then
        write(nu_diag,*) '>>>--------------------------------------------->>>'
        write(nu_diag,*) '                         i, j = ', i, j
        write(nu_diag,*) '                      -> tmt(i)=1'
       endif
       if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
        ifloe(i,j) = floe_sz_pancake
        if (cmt.ne.0) then
         if (vfice(i,j).lt.tolh) then
          write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
         endif
         if (afice(i,j).lt.tolice) then
          write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
         endif
         write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
        endif ! END IF COMMENT
       endif  ! END IF h<tolh OR c<tolc
      endif ! ENDIF tmt(i)
     enddo ! ENDDO i=1,nblock

     wspec_row(:,:) = c0     ! reset dummy vector

     if (cmt.ne.0) then
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
      write(nu_diag,*) '                      -> do wave directions', j
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
     endif

     ! B2. Redistribute energy according to mean directions
      ! LH Boundary:
      i=1
       ! if the cell isn't land ...
     if (tmask(1,j).and.tmt(1).ne.1) then
      dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
    dum_sm0        = 4d0*(dum_sm0**0.5d0)
      if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       wspec_row(1,:)  = wspec_row(1,:) + &
      				(1d0-sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
       mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(1)*dum_sm0
       mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
       tmt_hld(1)      = 0
      endif ! ENDIF -pi/4<mwd<pi/4
      ! if wave energy needs to be advected to the east ...
      if (tmask(2,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
       wspec_row(2,:) = wspec_row(2,:) + &
       				(sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
       mwd_hld(1,2)   = mwd_hld(1,2) + mwd_row(1)*dum_sm0
       mwd_hld(2,2)   = mwd_hld(2,2) + dum_sm0
       tmt_hld(2)     = 0
       ! if wave energy needs to be advected to the west (wrap around vector) ...
      elseif (tmask(nx_block,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
       wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                             (sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(1,:)
       mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + mwd_row(1)*dum_sm0
       mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
       tmt_hld(nx_block)     = 0
      endif ! ENDIF (tmask(2,j).and.sinmwd_row(1).lt.c0)
     endif ! END IF TMASK
     ! RH Boundary:
     i=nx_block
      ! if the cell isn't land ...
     if (tmask(nx_block,j).and.tmt(nx_block).ne.1) then
      dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
    dum_sm0        = 4d0*(dum_sm0**0.5d0)
      if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                            (1d0-sin(2d0*mwd_row(nx_block))**2d0)*wspec_row_hld(nx_block,:)
       mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + &
                             mwd_row(nx_block)*dum_sm0
       mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
       tmt_hld(nx_block)     = 0
      endif ! ENDIF -pi/4<mwd<pi/4
      ! if wave energy needs to be advected to the west ...
      if (tmask(nx_block-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
       wspec_row(nx_block-1,:) = wspec_row(nx_block-1,:) + &
     	                         (sin(2d0*mwd_row(nx_block))**2)*wspec_row_hld(nx_block,:)
       mwd_hld(1,nx_block-1)   = mwd_hld(1,nx_block-1) + &
     							  mwd_row(nx_block)*dum_sm0
       mwd_hld(2,nx_block-1)   = mwd_hld(1,nx_block-1) + dum_sm0
       tmt_hld(nx_block-1)     = 0
      ! if wave energy needs to be advected to the east ...
      elseif (tmask(1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
       wspec_row(1,:) = wspec_row(1,:) + &
       				(sin(2d0*mwd_row(nx_block))**2)*wspec_row_hld(nx_block,:)
       mwd_hld(1,1)   = mwd_hld(1,1) + mwd_row(nx_block)*dum_sm0
       mwd_hld(2,1)   = mwd_hld(2,1) + dum_sm0
       tmt_hld(1)     = 0
      endif ! ENDIF (tmask(nx_block-1,j).and.sinmwd_hld(nx_block).gt.c0)
     endif ! ENDIF TMASK
     ! loop the inner cells...
     do i=2,nx_block-1
      if (tmask(i,j).and.tmt(i).ne.1) then
       dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
     dum_sm0        = 4d0*(dum_sm0**0.5d0)
       if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       ! South
        wspec_row(i,:) = wspec_row(i,:) + &
       				(1d0-sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
        mwd_hld(1,i)   = mwd_hld(1,i) + mwd_row(i)*dum_sm0
        mwd_hld(2,i)   = mwd_hld(2,i) + dum_sm0
        tmt_hld(i)     = 0
       endif ! ENDIF -pi/4<mwd<pi/4
       ! East
       if (tmask(i+1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
        wspec_row(i+1,:) = wspec_row(i+1,:) + &
        					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
        mwd_hld(1,i+1)   = mwd_hld(1,i+1) + mwd_row(i)*dum_sm0
        mwd_hld(2,i+1)   = mwd_hld(2,i+1) + dum_sm0
        tmt_hld(i+1)     = 0
       elseif (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
        wspec_row(i-1,:) = wspec_row(i-1,:) + &
        					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
        mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
        mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
        tmt_hld(i-1)     = 0
       endif ! ENDIF (tmask(i+1,j).and.sinmwd_row(i).lt.c0)
      endif ! ENDIF TMASK
     end do ! ENDDO i=2,nx_block-1

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> update mean values', j
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
    endif

    ! B3. Calculate mean parameters
     do i=1,nx_block
      if (tmask(i,j)) then
       if (mwd_hld(2,i).gt.c0) then
        mwd_row(i) = mwd_hld(1,i)/mwd_hld(2,i)
       else
        mwd_row(i) = c0
       endif
       loc_mwd(i,j) = mwd_row(i)
       tmt(i)         = tmt_hld(i)
       if (tmt(i).eq.1) then
        loc_swh(i,j)   = 0d0
        loc_ppd(i,j)   = 0d0
       else
	  dum_sm0        = &
         fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
	  dum_sm2        = &
	   fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,2,nu_diag)
	  loc_swh(i,j)   = 4d0*(dum_sm0**0.5d0)
	  loc_ppd(i,j)   = 2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
       endif

       if (cmt.ne.0) then
        write(nu_diag,*) '         OUTPUT       -> i,j      =', i, j
        write(nu_diag,*) '                      -> tmt      =', tmt(i)
        write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
        write(nu_diag,*) '                      -> swh      =', loc_swh(i,j)
        write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j)
        write(nu_diag,*) '                      -> mwd      =', 180d0*loc_mwd(i,j)/pi
        write(nu_diag,*) '<<<---------------------------------------------<<<'
       endif
       endif ! END IF tmask(i,j)
      enddo ! END do i=1,nx_block

     wspec_row_hld(:,:) = c0     ! reset the dummy dummy vector
     mwd_hld(:,:)       = c0     ! ditto
     tmt_hld(:)         = 1

    enddo ! END do j=dum_wavemask-1,2,-1

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> finished increment_floe'
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
    endif





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

!!!!!!!!!!!!!!!!!!!!!
!!! fn_SpecMoment !!!
!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_SpecMoment
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!

function fn_SpecMoment(dum_S, nw_in, nth_in, om_in, mom, idl)

!
! !USES:
!
!
! !INPUT PARAMETERS:
!
integer, intent(in)                                :: nw_in, nth_in, mom, idl
real (kind=8), dimension(nw_in), intent(in)        :: om_in
!real (kind=8), dimension(nth_in), intent(in)       :: th_in
real (kind=8), dimension(nw_in*nth_in), intent(in) :: dum_S
!
! !OUTPUT PARAMETERS
!
real (kind=8)                            :: fn_SpecMoment
!
!EOP
!

integer                                  :: loop_w, loop_th
real (kind=8), dimension(nw_in*nth_in)   :: wt_simp, wt_int
real (kind=8), dimension(nw_in)          :: dum_simp
real (kind=8), dimension(nth_in)         :: dum_simp_th
real (kind=8), dimension(nw_in*nth_in)   :: dum_v
real (kind=8)							  :: dom_local, dth_local

dom_local = om_in(2)-om_in(1)
! Noah Day WIM,..
!if (nth_in.eq.1) then
! dth_local = 3d0         ! for Simpson's rule
!else
! dth_local = th_in(2)-th_in(1)
!end if
!-----
! if (idl.ne.0) then
!  write(idl,*) '                     --> into fn_SpecMoment...'
!  write(idl,*) '                         nw=', nw_in
!  write(idl,*) '                         th=', nth_in
!  write(idl,*) '                         mom=', mom
!  write(idl,*) '                         om=', om_in(1), '->', om_in(nw_in)
!  write(idl,*) '                         th=', th_in(1), '->', th_in(nth_in)
!  write(idl,*) '                         S=',  dum_S(1), '->', dum_S(nw_in*nth_in)
!  write(idl,*) '                         dom=', dom_local
!  write(idl,*) '                         dth=', dth_local
! end if

wt_int = c1
dum_v = c1

fn_SpecMoment   = dot_product(wt_int,dum_v)

! wt_simp(1)     = 1d0
! wt_simp(nw_in) = 1d0
!
! do loop_w=2,nw_in-1,2
!  wt_simp(loop_w) = 4d0
! end do
!
! do loop_w=3,nw_in-1,2
!  wt_simp(loop_w) = 2d0
! end do

!! if (idl.ne.0) then
!!  write(idl,*) '--> wt_simp', wt_simp(1), '->', wt_simp(nw_in)
!! end if

! do loop_w=1,nw_in
!  wt_int(loop_w) = (dom_local/3d0)*wt_simp(loop_w)
! end do

!! if (idl.ne.0) then
!!  write(idl,*) '--> wt_int', wt_int(1), '->', wt_int(nw_in)
!! end if

! fn_SpecMoment   = dot_product(wt_int,dum_S*(om_in**mom))

! if (idl.ne.0) then
!  write(idl,*) '... exiting fn_SpecMoment ', fn_SpecMoment
! end if

end function fn_SpecMoment

!=======================================================================

!=======================================================================

end module icepack_floe

!=======================================================================
