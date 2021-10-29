!=======================================================================
!BOP
!
! !MODULE: m_prams_waveice - define parameters for wave-ice interaction code
!
! !DESCRIPTION:
!
! Defines variable precision for all common data types \\
!
! !REVISION HISTORY:
!  SVN:$Id: m_waveice_prams.F90 2013-08-22 13:48:20 lgb566 $
!
! author: Luke Bennetts, Uni Adelaide
!
! !INTERFACE:
!
 module m_prams_waveice
!
! !USES:
!
 implicit none
!
!EOP
!=======================================================================

 ! BASIC PARAMS

 real(kind=8), parameter  :: pi = 3.141592653589793d0

 real(kind=8), parameter  :: gravity = 9.81d0               ! gravity
 real(kind=8), parameter  :: water_density = 1025d0         ! water density

 real(kind=8), parameter  :: reldens = 0.9d0                ! ice density/fluid density

 ! NUMERICAL PARAMS
 integer, parameter       :: WIM=1 		                ! waves in ice, (0 = off, 1 = on)

 integer, parameter       :: cmt=0 		                ! `write' parameter

 integer, parameter       :: do_coupled=0                   ! attn/floe size coupled?

 integer, parameter       :: Wsq_METH=0                     ! strain method
 															! 0 => W=1
 															! 1 => W=k_ice/k_wtr
 integer, parameter       :: ATTEN_METH=0                   ! attenuation method (0=Av,1=int)
 integer, parameter       :: ATTEN_MODEL=0                  ! attenuation model
  															! 0=fn_apxAttn
  															! 1=fn_Attn_WIM_v1
 integer, parameter       :: res=500                        ! resolution of floe size interval

 real(kind=8), parameter  :: tolice = 1.0d-2               ! conc<tolice treated as zero ice
 real(kind=8), parameter  :: tolh   = 1.0d-1               ! h<tolh      treated as zero ice
 real(kind=8), parameter  :: toli   = 1.0d-16

 ! ATTENUATION

 integer, parameter                           :: Ncheb_f=26 ! deg of poly in period/freq
 integer, parameter                           :: Ncheb_h=26 ! deg of poly in thickness

 real(kind=8), dimension(Ncheb_f+1,Ncheb_h+1) :: alp_coeffs

 real(kind=8), parameter                      :: attn_fac=1.0d0 ! adjust attenuation rate

 ! ICE PARAMS

 real(kind=8), parameter  :: Y0=10d0**10                    ! Young's modulus parameter
 real(kind=8), parameter  :: nub=0.1d0                      ! brine volume content
 real(kind=8), parameter  :: Y=Y0*(1-3.51*nub) - 10d0**9    ! eff Young's mod
 real(kind=8), parameter  :: poisson=0.3d0                  ! Poisson's ratio

 real(kind=8), parameter  :: bk0=1.76d0*10**6, bk1=5.88d0   ! breaking strain params
 real(kind=8), parameter  :: Pc = exp(-1d0)                 ! critical strain
 real(kind=8), parameter  :: epsc=bk0*exp(-bk1*sqrt(nub))/Y ! breaking strain

 ! FSD params
 real(kind=8), parameter  :: floe_sz_min = 20d0, floe_sz_crit=25d0
 real(kind=8), parameter  :: i0 = floe_sz_min+0.1d0, i1 = 500d0
 real(kind=8), parameter  :: gamma0=1.15d0, gamma1=2.1d0
 real(kind=8), parameter  :: gamma =2d0+log(0.9d0)/log(2d0)
 real(kind=8), parameter  :: floe_sz_brk_emp = 20d0
 real(kind=8), parameter  :: floe_sz_pancake = 5d0

 ! WAVES

 integer, parameter                    :: WAVE_METH = 1   ! inc waves (0=user,1=ww3)
 real(kind=8), allocatable             :: ww3_lat(:,:), ww3_lon(:,:), ww3_tm(:,:)
 real(kind=8), allocatable             :: ww3_swh(:,:), ww3_fp(:,:), ww3_dir(:,:), &
 				ww3_swh_full(:,:,:), ww3_fp_full(:,:,:), ww3_dir_full(:,:,:)
 integer, parameter                    :: nww3_dt = 1 ! ww3 time step relative to cice

 ! DATA

 character(100)            :: waveicedatadir='/home/cawcr_data/gx1/2005'!19 '../../waveice_data/'
 character(10)            :: fname_alp='alp_coeffs'
 !character(24)            :: fname_ww3='waves/ww3.197803_full.nc'
 character(14)            :: fname_ww3='ww3_2005'!'waves/ww3.1978'
 integer, parameter       :: OVERWRITE_DIRS = 1   ! overwrite wave directions with usr set ones (0=no,1=yes)
 integer, parameter       :: WIM_BREAKUP = 0   ! use the WIM breakup method (0=no,1=yes)

!=======================================================================

 end module m_prams_waveice

!=======================================================================
