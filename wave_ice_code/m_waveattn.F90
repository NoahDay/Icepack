!=======================================================================
!BOP
!
! !MODULE: m_waveattn - attenuation models
!
! !DESCRIPTION:
!
! blah blah blah \\
!
! !REVISION HISTORY:
!  SVN:$Id: m_waveattn.F90 2013-09-02 13:48:20 lgb566 $
!
! author: Luke Bennetts, Uni Adelaide
!
! !INTERFACE:
!
 module m_waveattn
!
! !USES:

 use m_prams_waveice

 implicit none

! !PUBLIC MEMBER FUNCTIONS:
!
!EOP
!=======================================================================

 ! PARAMETERS:

 contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!! Subfunctions & subroutines !!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!
 !!! fn_Attn_WIM_v1 !!!
 !!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_Attn_MBK
!
! !DESCRIPTION:
!
!  The attenuation model derived by Meylan et al. (2014), Geophys Res Lett
!
!  alpha = beta0 + beta1*om^2
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_Attn_MBK(dum_om)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real(kind=8), intent (in) :: dum_om        ! ang freq

!
! !OUTPUT PARAMETERS
!

 real(kind=8) :: fn_Attn_MBK

!
!EOP
!

 real(kind=8), parameter :: beta0 = 5.376168295200780E-005, &
     beta1 = 2.947870279251530E-005

 fn_Attn_MBK = beta0*(dum_om**2) + beta1*(dum_om**4)

 fn_Attn_MBK = attn_fac*fn_Attn_MBK

 end function fn_Attn_MBK

 !!!!!!!!!!!!!!!!!!!!!!
 !!! fn_Attn_WIM_v1 !!!
 !!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_Attn_WIM_v1
!
! !DESCRIPTION:
!
!  The attenuation model used for Williams et al (2013a,b) Ocean Model.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_Attn_WIM_v1(dum_om,dum_h,dum_idl)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real(kind=8), intent (in) :: dum_om, dum_h  ! ang freq, thickness

 integer, intent(in)       :: dum_idl

!
! !OUTPUT PARAMETERS
!

 real(kind=8) :: fn_Attn_WIM_v1

!
!EOP
!

 integer                    :: lp_i, lp_j !, idd_alp
 real(kind=8)               :: new_h      ! dum_alp,
 real(kind=8), dimension(4) :: Limits
! integer, parameter         :: Ncheb_f=26 ! deg of poly in period/freq
! integer, parameter         :: Ncheb_h=26 ! deg of poly in thickness
 real(kind=8), parameter    :: hmax=8.25d0,hmin=0.25d0
 real(kind=8), parameter    :: xmin=2*pi/35d0,xmax=2.692793703076966d0

! real(kind=8), dimension(Ncheb_f+1,Ncheb_h+1) :: alp_coeffs
!
! character(30) :: fname_alp
! character(len=255) :: cwd

! if (dum_idl.ne.0) then
!  call getcwd(cwd)
!  write(dum_idl,*) '                      --> Into fn_Attn_WIM_v1'
!  write(dum_idl,*) '                      --> ', trim(cwd)
!  write(dum_idl,*) '                      --> omega = ', dum_om
!  write(dum_idl,*) '                      --> h     = ', dum_h
! endif

! fname_alp='data/alp_coeffs'
!
! open(newunit=idd_alp,file=fname_alp)
!
! do lp_j=1,Ncheb_h+1
!  do lp_i=1,Ncheb_f+1
!   read(idd_alp,*) dum_alp
!   alp_coeffs(lp_i,lp_j)=dum_alp
!  end do
! end do
!
! close(idd_alp)

 !print*, 'alp_coeffs:', alp_coeffs(1,1), alp_coeffs(Ncheb_f+1,1), &
 ! alp_coeffs(1,Ncheb_h+1), alp_coeffs(Ncheb_f+1,Ncheb_h+1)

 if (dum_h.lt.hmin) then
  new_h=hmin
 elseif (dum_h.gt.hmax) then
  new_h=hmax
 else
  new_h=dum_h
 end if

 !print*, 'hnew=',new_h

 Limits=(/ xmin,xmax,hmin,hmax /)

 !print*, 'Limits',Limits

 fn_Attn_WIM_v1=&
  fn_OP_chebinterp2d(dum_om,new_h,Ncheb_f,Ncheb_h,alp_coeffs,Limits)

 fn_Attn_WIM_v1 = attn_fac*fn_Attn_WIM_v1

 !print*, 'alpha=',fn_Attn_WIM_v1

 end function fn_Attn_WIM_v1

!=======================================================================

 !!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! fn_OP_chebinterp2d !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_OP_chebinterp2d
!
! !DESCRIPTION:
!
!  Chebyshev interpolation of a function of 2 variables
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_OP_chebinterp2d(xx,yy,Ncheb1,Ncheb2,F_coeffs,Lims)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 integer, intent(in)                                     :: Ncheb1, Ncheb2
 real(kind=8), intent (in)                               :: xx, yy
 real(kind=8), dimension(Ncheb1+1,Ncheb2+1), intent (in) :: F_coeffs
 real(kind=8), dimension(4), intent (in)                 :: Lims

!
! !OUTPUT PARAMETERS
!

 real(kind=8) :: fn_OP_chebinterp2d

!
!EOP
!

 real(kind=8) :: x_min,x_max,y_min,y_max,dx,dy,tt1,tt2

 real(kind=8), dimension(:), allocatable :: TnVals1, TnVals2

 !print*, 'into cheby2'

 x_min    = Lims(1)
 x_max    = Lims(2)
 y_min    = Lims(3)
 y_max    = Lims(4)

 allocate(TnVals1(1:Ncheb1+1))
 allocate(TnVals2(1:Ncheb2+1))

 dx       = x_max - x_min
 tt1      = -1d0 + 2d0*(xx-x_min)/dx
 TnVals1  = fn_OP_chebinterp1d(tt1,Ncheb1)

 !print*, 'tt1=',tt1
 !print*,'TnVals1=',TnVals1(1),'->',TnVals1(Ncheb1+1)

 dy       = y_max - y_min
 tt2      = -1d0 + 2d0*(yy-y_min)/dy
 TnVals2  = fn_OP_chebinterp1d(tt2,Ncheb2)

 fn_OP_chebinterp2d = dot_product(TnVals1,matmul(F_coeffs,TnVals2))

 deallocate(TnVals1)
 deallocate(TnVals2)

 end function fn_OP_chebinterp2d

!=======================================================================

 !!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! fn_OP_chebinterp1d !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_OP_chebinterp1d
!
! !DESCRIPTION:
!
!  Chebyshev interpolation of a function of 1 variable
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_OP_chebinterp1d(tt,An)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real(kind=8), intent (in) :: tt
 integer, intent (in)      :: An

!
! !OUTPUT PARAMETERS
!

 real(kind=8), dimension(An+1) :: fn_OP_chebinterp1d

!
!EOP
!

 integer                        :: lp
 real(kind=8), dimension(An+1)  :: hn, f
 real(kind=8)                   :: C0, C1, Cn

 !print*, 'into cheby1'

 hn(1) = pi
 do lp=1,An
  hn(lp+1)=pi/2
 end do

 !print*,'hn=',hn(1),hn(2),'->',hn(An+1)

 C0=1d0
 C1=tt
 f(1)=C0
 f(2)=C1

 do lp=2,An
  Cn      = 2*tt*C1-C0
  f(lp+1) = Cn
  C0      = C1
  C1      = Cn
 end do

 !print*,'f=',f(1),f(2),'->',f(An+1)

 fn_OP_chebinterp1d = f

 end function fn_OP_chebinterp1d

!=======================================================================

!=======================================================================

 end module m_waveattn

!=======================================================================
