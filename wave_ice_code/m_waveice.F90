!=======================================================================
!BOP
!
! !MODULE: m_waveice - subfunctions required by m_balance_waveice
!
! !DESCRIPTION:
!
! blah blah blah \\
!
! !REVISION HISTORY:
!  SVN:$Id: m_waveice.F90 2013-09-02 13:48:20 lgb566 $
!
! author: Luke Bennetts, Uni Adelaide
!
! !INTERFACE:
!
 module m_waveice
!
! !USES:

 use m_fzero
 use m_prams_waveice
 use m_waveattn

 implicit none

! !PUBLIC MEMBER FUNCTIONS:
!
!EOP
!=======================================================================

 ! PARAMETERS:

 integer                                  :: nw, nth
 real(kind=8), dimension(:), allocatable  :: om, th, k_wtr, S_init

 real(kind=8), dimension(:), allocatable  :: k_ice, lam_ice     ! ice-coupled waveno/length

 real(kind=8)                             :: hice, conc

 real(kind=8)                             :: dom, dth

 real(kind=8)                             :: kappa              ! frequency parameter

 real(kind=8)                             :: L                  ! MIZ length

 integer                                  :: dum_idl

 real(kind=8)                             :: mass, fr           ! scaled mass & flexural rigidity

 contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!! Subfunctions & subroutines !!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!
 !!! sub_Uncoupled !!!
 !!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: sub_Uncoupled
!
! !DESCRIPTION:
!
!  attenuation does not depend on floe sizes
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!

 subroutine sub_Uncoupled(floe_sz_init,floe_sz_max,Lcell,dum_hice,dum_conc, &
  dum_nw,dum_nth,dum_om,dum_th,dum_k_wtr,dum_S_init,S_attn_out,tmt,idl)

!
! !USES:
!
! !INPUT: Log file param
!

 integer, intent(in)                             :: idl
 real(kind=8), intent(in)            	 :: floe_sz_init          ! initial max floe size
 real(kind=8), intent(in)                 :: Lcell
 real(kind=8), intent(in)                 :: dum_hice,dum_conc     ! cell length, ice thick & con

 integer, intent(in)                                 :: dum_nw, dum_nth
 real(kind=8), dimension(dum_nw), intent(in)  :: dum_om, dum_k_wtr, dum_S_init
 real(kind=8), dimension(dum_nth), intent(in) :: dum_th

!
! !OUTPUT:
!
 real(kind=8), intent(out)          		   :: floe_sz_max     ! max floe size
 real(kind=8), dimension(dum_nw), intent(out)  :: S_attn_out
 integer, intent(out)                          :: tmt             ! flag to terminate routine
!
!EOP
!

 ! Numerical params:

 integer                      :: lp_i, lp_j

 ! Ice params:

 real(kind=8)                 :: Es, ES_init              ! strain

 ! MIZ length parameter

 real(kind=8)                 :: Lc

 ! Wave params:

 real(kind=8), dimension(:), allocatable  :: S_attn, alpha   ! attn spectrum & attn coeff

 real(kind=8)                 :: lam_init

 if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '>>>--------------------------------------------->>>'
  write(idl,*) '        sub_Uncoupled  -->  h  = ',dum_hice
  write(idl,*) '                       -->  c  = ',dum_conc
  write(idl,*) '                       -->  Di = ',floe_sz_init
  write(idl,*) '                       -->  L  = ',Lcell/1000d0,'km'
 end if

 dum_idl = idl

 tmt     = 0

 nw   = dum_nw
 nth  = dum_nth
 conc = dum_conc
 hice = dum_hice

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!! No ice = nothing to do
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (conc.lt.tolice.or.hice.lt.tolh) then
  !Noah Day 25/10 floe_sz_max = floe_sz_pancake
  do lp_i=1,nw
   S_attn_out(lp_i) = dum_S_init(lp_i)
  end do

  if (idl.ne.0.and.cmt.ne.0) then
   write(idl,*) '>>>>>>> RESULT:'
   if (hice.eq.0d0) then
    write(idl,*) '                       --> no ice in this cell: h<', tolh
   else
    write(idl,*) '                       --> no ice in this cell: c<', tolice
   end if
  write(idl,*) '                           tmt =', tmt
  write(idl,*) '                   floe_sz_max =', floe_sz_max
!   write(idl,*) '<<<---------------------------------------------<<<'
  end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else !!! there is some ice !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 mass  = reldens*hice
 fr    = Y*(hice**3)/12/water_density/gravity/(1-poisson**2)

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>>>> PARAMS:'
  write(idl,*) 'mass, fr = ', mass, fr
  write(idl,*) 'E_crit=epsc*sqrt(-2/log(Pc))      = ', epsc*sqrt(-2/log(Pc))
  write(idl,*) 'nw, nth', nw, nth
  write(idl,*) 'om=',dum_om(1),'->',dum_om(nw)
  write(idl,*) 'th=',dum_th(1),'->',dum_th(nth)
 endif

! if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '>>>>>>> ALLOCATING VARIABLES...'
! endif

 allocate(om(1:nw))
 allocate(th(1:nth))

 allocate(k_wtr(1:nw))

 allocate(k_ice(1:nw))
 allocate(lam_ice(1:nw))

 allocate(S_init(1:nth*nw))
 allocate(S_attn(1:nw*nth))

 allocate(alpha(1:nw))

! if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '         ... DONE'
! endif

 do lp_i=1,nw
  om(lp_i)     = dum_om(lp_i)
  k_wtr(lp_i)  = dum_k_wtr(lp_i)
 end do

! if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '         ... OM & K_WTR CALCULATED'
! endif

 do lp_i=1,nth*nw
  S_init(lp_i) = dum_S_init(lp_i)
 end do

! if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '         ... S_INIT CALCULATED'
! endif

 do lp_i=1,nth
  th(lp_i)     = dum_th(lp_i)
 end do

! if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '         ... TH CALCULATED'
! endif

 dom = om(2)-om(1)
 if (nth.eq.1) then
  dth = 3d0         ! for Simpson's rule
 else
  dth = th(2)-th(1)
 end if

! if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '         ... DTH CALCULATED'
! endif

 do lp_i=1,nw
  kappa           = om(lp_i)**2d0/gravity
  k_ice(lp_i)     = zero(0d0,max(kappa/tanh(kappa),sqrt(sqrt(kappa*mass/fr))), &
                         toli,toli,fn_DispRel_ice_inf)
  lam_ice(lp_i)   = 2d0*pi/k_ice(lp_i)
 end do

! if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '         ... KAPPA, K_ICE & LAM_ICE CALCULATED'
! endif

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>> WAVE NUMBER SPECTRA:'
  write(idl,*) 'k_wtr   = ', k_wtr(1), '->', k_wtr(nw)
  write(idl,*) 'lam_wtr = ', 2*pi/k_wtr(1), '->', 2*pi/k_wtr(nw)
  write(idl,*) 'k_ice   = ', k_ice(1), '->', k_ice(nw)
  write(idl,*) 'lam_ice = ', lam_ice(1), '->', lam_ice(nw)
 endif

 L = 0  ! for test of incident spectrum

if (WIM_BREAKUP.eq.1) then
 call sub_StrainSpec(S_init, Es)
else
  Es = 0d0
end if

 call sub_WavelenSpec(S_init, lam_init)

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>> INCIDENT SPECTRUM:'
  write(idl,*) 'S_init   = ', S_init(1), '->', S_init(nw*nth)
  call sub_StrainSpec(S_init, Es_init)
  write(idl,*) 'lam_init = ', lam_init
  !write(idl,*) 'Es_init  = ', Es_init
 end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!! Inc strain<fracture threshold =
 !!!!! nothing to do but advect waves
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (Es.lt.epsc*sqrt(-2/log(Pc))) then

 L = Lcell ! initialise propagation length

 !Noah Day 25/10 floe_sz_max = floe_sz_init

 do lp_i=1,nw
  alpha(lp_i)  = conc*fn_Attn_MBK(om(lp_i))/0.7d0
  do lp_j=1,nth
   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
     exp(-alpha(lp_i)*cos(th(lp_j))*L)
  end do
 end do

 do lp_i=1,nw*nth
  S_attn_out(lp_i) = S_attn(lp_i)
 end do

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>> RESULT:'
  write(idl,*) '                      --> S_init isnt strong enough to break ice: ', &
   Es, epsc*sqrt(-2/log(Pc))
  write(idl,*) '                          tmt =', tmt
  write(idl,*) '                  floe_sz_max =', floe_sz_max
!   write(idl,*) '<<<---------------------------------------------<<<'
 end if

 !tmt = 0      ! Es proportional to hice -> so can't just switch off (LB Apr 14)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else !!! potential for fracture
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Noah Day 25/10  L = Lcell ! initialise propagation length

!Noah Day 25/10  if (idl.ne.0.and.cmt.ne.0) then
!Noah Day 25/10   write(idl,*) '>>>>>>> AT END OF CELL:'
  !write(idl,*) 'D1 (init) = ', floe_sz_init, ': empirical breaking size = ', floe_sz_brk_emp
 !Noah Day 25/10 end if

! Propagation of wave spectrum
 do lp_i=1,nw
  alpha(lp_i)  = conc*fn_Attn_MBK(om(lp_i))/0.75d0
  do lp_j=1,nth
   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
     exp(-alpha(lp_i)*cos(th(lp_j))*L)
  end do
 end do

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) 'alpha = ', alpha(1), '->', alpha(nw)
  write(idl,*) 'S_attn = ', S_attn(1), '->', S_attn(nw)
 end if

 if (WIM_BREAKUP.eq.1) then
  call sub_StrainSpec(S_init, Es)
 else
   Es = 0d0
 end if

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) 'Es = ', Es
 end if

 ! Set attenuated spectrum (as not floe size dependent)
 do lp_i=1,nw*nth
  S_attn_out(lp_i) = S_attn(lp_i)
 end do

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Noah Day 25/10 if (floe_sz_brk_emp.gt.floe_sz_init) then  !!! floes already
 !!!                             smaller than breakage size
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Noah Day 25/10 floe_sz_max = floe_sz_init
  !tmt = 0
!Noah Day 25/10  if (idl.ne.0.and.cmt.ne.0) then
!Noah Day 25/10   write(idl,*) '>>>>>>> RESULT:'
!Noah Day 25/10   write(idl,*) '                       --> Breakage occurs: location inconsequential'
!Noah Day 25/10   write(idl,*) '                           initial floe size', floe_sz_init, '<', &
!Noah Day 25/10      											floe_sz_brk_emp
!Noah Day 25/10   write(idl,*) '                           tmt =', tmt
!Noah Day 25/10   write(idl,*) '                   floe_sz_max =', floe_sz_max
!   write(idl,*) '<<<---------------------------------------------<<<'
!Noah Day 25/10  end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Noah Day 25/10 elseif (Es.ge.epsc*sqrt(-2/log(Pc))) then   !!! waves strong
 !!!                   enough to break ice at end of cell
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Noah Day 25/10   floe_sz_max = floe_sz_brk_emp
  !tmt = 0
!Noah Day 25/10   if (idl.ne.0.and.cmt.ne.0) then
!Noah Day 25/10    write(idl,*) '>>>>>>> RESULT:'
!Noah Day 25/10    write(idl,*) '                       --> Breakage occurs at end of cell',L/1000d0,'km'
!Noah Day 25/10    write(idl,*) '                           tmt =', tmt
!Noah Day 25/10    write(idl,*) '                   floe_sz_max =', floe_sz_max
!   write(idl,*) '<<<---------------------------------------------<<<'
!Noah Day 25/10   end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Noah Day 25/10 else        !!! search for point at which waves can break ice
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Noah Day 25/10 Lc = zero(i0,Lcell,toli,toli,fn_proplng_uncoupled)
  !Noah Day 25/10 L           = Lc
  !Noah Day 25/10 floe_sz_max = (Lc*floe_sz_brk_emp + (Lcell-Lc)*floe_sz_init)/Lcell
  !Noah Day 25/10 tmt         = 1   ! waves attenuated (in all liklihood)
!Noah Day 25/10  if (idl.ne.0.and.cmt.ne.0) then
!Noah Day 25/10   write(idl,*) '>>>>>>> RESULT:'
!Noah Day 25/10   write(idl,*) '                      --> Breakage occurs within cell ', &
!Noah Day 25/10   Lc/1000d0,'km' !, '(zero chk: ', fn_proplng_uncoupled(Lc),')'
!Noah Day 25/10   write(idl,*) '                          tmt =', tmt
!Noah Day 25/10   write(idl,*) '                  floe_sz_max =', floe_sz_max
!   write(idl,*) '<<<---------------------------------------------<<<'
!Noah Day 25/10  end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Noah Day 25/10end if ! (Es.ge.epsc*sqrt(-2/log(Pc)))
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end if ! (Esinit.lt.Ec)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 deallocate(om)
 deallocate(th)

 deallocate(k_wtr)

 deallocate(k_ice)
 deallocate(lam_ice)

 deallocate(S_init)
 deallocate(S_attn)

 deallocate(alpha)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end if  ! (conc.lt.tolice.or.hice.eq.0d0)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 end subroutine sub_Uncoupled

 !!!!!!!!!!!!!!!!!!!
 !!! sub_Balance !!!
 !!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: sub_Balance
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!

 subroutine sub_Balance(floe_sz_init,floe_sz_max,Lcell,dum_hice,dum_conc, &
  dum_nw,dum_nth,dum_om,dum_th,dum_k_wtr,dum_S_init,S_attn_out,tmt,idl)

!
! !USES:
!
! !INPUT: Log file param
!

 integer, intent(in)                      :: idl
 real(kind=8), intent(in)            	  :: floe_sz_init          ! initial max floe size
 real(kind=8), intent(in)                 :: Lcell
 real(kind=8), intent(in)                 :: dum_hice,dum_conc     ! cell length, ice thick & con

 integer, intent(in)                          :: dum_nw, dum_nth
 real(kind=8), dimension(dum_nw), intent(in)  :: dum_om, dum_k_wtr, dum_S_init
 real(kind=8), dimension(dum_nth), intent(in) :: dum_th

!
! !OUTPUT:
!
 real(kind=8), intent(out)          		   :: floe_sz_max     ! max floe size
 real(kind=8), dimension(dum_nw), intent(out)  :: S_attn_out
 integer, intent(out)                          :: tmt             ! flag to terminate routine
!
!EOP
!

 ! Numerical params:

 integer                      :: lp_i, lp_j

 ! Ice params:

 real(kind=8)                 :: Es, ES_init              ! strain

 ! MIZ length parameter

 real(kind=8)                 :: Lc

 ! Wave params:

 real(kind=8), dimension(:), allocatable  :: S_attn, alpha   ! attn spectrum & attn coeff

 real(kind=8)                 :: lam_init

 if (idl.ne.0.and.cmt.ne.0) then
!  write(idl,*) '>>>--------------------------------------------->>>'
  write(idl,*) '        sub_Balanace   -->  h  = ',dum_hice
  write(idl,*) '                       -->  c  = ',dum_conc
  write(idl,*) '                       -->  Di = ',floe_sz_init
  write(idl,*) '                       -->  L  = ',Lcell/1000d0,'km'
 end if

 dum_idl = idl

 tmt     = 0

 nw   = dum_nw
 nth  = dum_nth
 conc = dum_conc
 hice = dum_hice

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!! No ice = nothing to do
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (conc.lt.tolice.or.hice.lt.tolh) then
  floe_sz_max = 5d0
  do lp_i=1,nw
   S_attn_out(lp_i) = dum_S_init(lp_i)
  end do

  if (idl.ne.0.and.cmt.ne.0) then
   write(idl,*) '>>>>>>> RESULT:'
   if (hice.eq.0d0) then
    write(idl,*) '                       --> no ice in this cell: h<', tolh
   else
    write(idl,*) '                       --> no ice in this cell: c<', tolice
   end if
!   write(idl,*) '<<<---------------------------------------------<<<'
  end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else !!! there is some ice !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 mass  = reldens*hice
 fr    = Y*(hice**3)/12/water_density/gravity/(1-poisson**2)

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>>>> PARAMS:'
  write(idl,*) 'mass, fr = ', mass, fr
  write(idl,*) 'E_crit=epsc*sqrt(-2/log(Pc))      = ', epsc*sqrt(-2/log(Pc))
  write(idl,*) 'nw, nth', nw, nth
  write(idl,*) 'om=',dum_om(1),'->',dum_om(nw)
  write(idl,*) 'th=',dum_th(1),'->',dum_th(nth)
 endif

 allocate(om(1:nw))
 allocate(th(1:nth))

 allocate(k_wtr(1:nw))

 allocate(k_ice(1:nw))
 allocate(lam_ice(1:nw))

 allocate(S_init(1:nth*nw))
 allocate(S_attn(1:nw*nth))

 allocate(alpha(1:nw))

 do lp_i=1,nw
  om(lp_i)     = dum_om(lp_i)
  k_wtr(lp_i)  = dum_k_wtr(lp_i)
 end do

 do lp_i=1,nth*nw
  S_init(lp_i) = dum_S_init(lp_i)
 end do

 do lp_i=1,nth
  th(lp_i)     = dum_th(lp_i)
 end do

 dom = om(2)-om(1)
 if (nth.eq.1) then
  dth = 3d0         ! for Simpson's rule
 else
  dth = th(2)-th(1)
 end if

 do lp_i=1,nw
  kappa           = om(lp_i)**2d0/gravity
  k_ice(lp_i)     = zero(0d0,max(kappa/tanh(kappa),sqrt(sqrt(kappa*mass/fr))), &
                         toli,toli,fn_DispRel_ice_inf)
  lam_ice(lp_i)   = 2d0*pi/k_ice(lp_i)
 end do

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>> WAVE NUMBER SPECTRA:'
  write(idl,*) 'k_wtr   = ', k_wtr(1), '->', k_wtr(nw)
  write(idl,*) 'lam_wtr = ', 2*pi/k_wtr(1), '->', 2*pi/k_wtr(nw)
  write(idl,*) 'k_ice   = ', k_ice(1), '->', k_ice(nw)
  write(idl,*) 'lam_ice = ', lam_ice(1), '->', lam_ice(nw)
 endif

 L = 0  ! for test of incident spectrum

 call sub_StrainSpec(S_init, Es)

 call sub_WavelenSpec(S_init, lam_init)

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>> INCIDENT SPECTRUM:'
  write(idl,*) 'S_init   = ', S_init(1), '->', S_init(nw*nth)
  call sub_StrainSpec(S_init, Es_init)
  write(idl,*) 'lam_init = ', lam_init
  write(idl,*) 'Es_init  = ', Es_init
 end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!! Inc strain<fracture threshold =
 !!!!! nothing to do but advect waves
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (Es.lt.epsc*sqrt(-2/log(Pc))) then

!!!!!!!!!!!! >CAN BE COMMENTED OUT< !!!!!!!!!!!!!!!!!!!!
!
! floe_sz_max = zero(i0,i1,toli,toli,fn_wlng)
!
! if (idl.ne.0) then
!  write(idl,*) 'lam_init = ', lam_init
!  write(idl,*) 'D1 = ', floe_sz_max, ': f(z) = ', fn_wlng(floe_sz_max)
! end if
!
!!!!!!!!!!!! <CAN BE COMMENTED OUT> !!!!!!!!!!!!!!!!!!!!

 L = Lcell ! initialise propagation length

 floe_sz_max = floe_sz_init

 do lp_i=1,nw
  if (ATTEN_METH.eq.1) then
   alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
  else
   alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
  end if
!  S_attn(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
  do lp_j=1,nth
   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
     exp(-alpha(lp_i)*cos(th(lp_j))*L)
  end do
 end do

 do lp_i=1,nw*nth
  S_attn_out(lp_i) = S_attn(lp_i)
 end do

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>> RESULT:'
  write(idl,*) '                      --> S_init isnt strong enough to break ice: ', &
   Es, epsc*sqrt(-2/log(Pc))
!   write(idl,*) '<<<---------------------------------------------<<<'
 end if

 tmt = 0      ! Es proportional to hice -> so can't just switch off (LB Apr 14)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else !!! potential for fracture
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 L = Lcell ! initialise propagation length

 call sub_WavelenSpec(S_init, lam_init)

 floe_sz_max = zero(i0,i1,toli,toli,fn_wlng)

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) '>>>>>>> AT END OF CELL:'
  write(idl,*) 'D1 = ', floe_sz_max, ': f(z) = ', fn_wlng(floe_sz_max)
 end if

 do lp_i=1,nw
  if (ATTEN_METH.eq.1) then
   alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
  else
   alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
  end if
!  S_attn(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
  do lp_j=1,nth
   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
     exp(-alpha(lp_i)*cos(th(lp_j))*L)
  end do
 end do

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) 'alpha = ', alpha(1), '->', alpha(nw)
  write(idl,*) 'S_attn = ', S_attn(1), '->', S_attn(nw)
 end if

 call sub_StrainSpec(S_attn, Es)

 if (idl.ne.0.and.cmt.ne.0) then
  write(idl,*) 'Es = ', Es
 end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (Es.ge.epsc*sqrt(-2/log(Pc))) then   !!! waves strong
 !!!                   enough to break ice at end of cell
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (idl.ne.0.and.cmt.ne.0) then
   write(idl,*) '>>>>>>> RESULT:'
   write(idl,*) '                       --> Breakage occurs at end of cell',L/1000d0,'km'
!   write(idl,*) '<<<---------------------------------------------<<<'
  end if
  do lp_i=1,nw*nth
   S_attn_out(lp_i) = S_attn(lp_i)
  end do
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else        !!! search for point at which waves can break ice
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Lc = zero(i0,Lcell,toli,toli,fn_proplng)
  if (idl.ne.0.and.cmt.ne.0) then
   write(idl,*) '>>>>>>> RESULT:'
   write(idl,*) '                      --> Breakage occurs within cell ', &
   Lc/1000d0,'km' !, '(zero chk: ', fn_proplng(Lc),')'
!   write(idl,*) '<<<---------------------------------------------<<<'
  end if
  L = Lc
  floe_sz_max = zero(i0,i1,toli,toli,fn_wlng)
  do lp_i=1,nw
   if (ATTEN_METH.eq.1) then
    alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
   else
    alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
   end if
!   S_attn_out(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
   do lp_j=1,nth
    S_attn_out(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
     exp(-alpha(lp_i)*cos(th(lp_j))*L)
   end do
  end do
  tmt = 1

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end if ! (Es.ge.epsc*sqrt(-2/log(Pc)))
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end if ! (Esinit.lt.Ec)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 deallocate(om)
 deallocate(th)

 deallocate(k_wtr)

 deallocate(k_ice)
 deallocate(lam_ice)

 deallocate(S_init)
 deallocate(S_attn)

 deallocate(alpha)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end if  ! (conc.lt.tolice.or.hice.eq.0d0)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 end subroutine sub_Balance

!=======================================================================

 !!!!!!!!!!!!!!!!!!
 !!! fn_apxAttn !!!
 !!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_apxAttn  - approximate the attenuation coefficient
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_apxAttn(fl_len,lambda,gamma)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real (kind=8), intent (in) :: fl_len, lambda, gamma

!
! !OUTPUT PARAMETERS
!

 real (kind=8) :: fn_apxAttn

!
!EOP
!

 real (kind=8)              :: Rsq, A, B

 Rsq = 0.002       ! Reflected energy (10s period, 1m thickness)
 A = -2*log(1-Rsq) ! long-floe limit

 B = 5/(lambda/2)**gamma ! tanh(5) = 1

 fn_apxAttn = A*tanh(B*(fl_len**gamma))

 end function fn_apxAttn

!=======================================================================

 !!!!!!!!!!!!!!!!!!!
 !!! fn_splitPDF !!!
 !!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_splitPDF  - the split PDF of Toyota et al (2011)
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_splitPDF(D,beta0,beta1,P0)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real (kind=8), intent(in) :: D, beta0, beta1, P0

!
! !OUTPUT PARAMETERS
!

 real (kind=8) :: fn_splitPDF

!
!EOP
!

 if (D.ge.floe_sz_crit) then
    fn_splitPDF = (1-P0)*beta1*gamma1/D**(gamma1+1d0)
 elseif (D.lt.floe_sz_crit.and.D.ge.floe_sz_min) then
    fn_splitPDF = P0*beta0*gamma0/D**(gamma0+1d0)
 else
    fn_splitPDF=0
 end if

 end function fn_splitPDF

!=======================================================================

 !!!!!!!!!!!!!!!!!
 !!! fn_AvAttn !!!
 !!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_AvAttn
!
! !DESCRIPTION:
!
!  the averaged attenuation coefficient wrt floe diameter
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_AvAttn(floe_sz_max,lam,dum_om)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real(kind=8), intent(in) :: floe_sz_max, lam, dum_om

!
! !OUTPUT PARAMETERS
!

 real(kind=8) :: fn_AvAttn

!
!EOP
!

 real(kind=8)             :: Dav, alpha_nd

 ! D average

 Dav = gamma*(floe_sz_max*((floe_sz_min/floe_sz_max)**gamma)-floe_sz_min)/(1d0-gamma)

 ! non-dim attenuation

 if (ATTEN_MODEL.eq.0) then
  alpha_nd = fn_apxAttn(Dav,lam,1d0)
 elseif (ATTEN_MODEL.eq.1) then
  alpha_nd = fn_Attn_WIM_v1(dum_om,hice,dum_idl)
 end if

 ! Dimensionalise

 fn_AvAttn = conc*alpha_nd/Dav

 end function fn_AvAttn

!=======================================================================

 !!!!!!!!!!!!!!!!!!
 !!! fn_IntAttn !!!
 !!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_IntAttn
!
! !DESCRIPTION:
!
!  the integrated attenuation coefficient wrt floe diameter
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_IntAttn(floe_sz_max, lam, dum_om)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real (kind=8), intent(in)       :: floe_sz_max, lam, dum_om

!
! !OUTPUT PARAMETERS
!

 real (kind=8) :: fn_IntAttn

!
!EOP
!

 real (kind=8), dimension(res+1) :: dum_D_vec
 real (kind=8), dimension(res)   :: D_vec, alpha_vec, p_vec
 real (kind=8)                   :: dD, beta0, beta1, P0
 integer                         :: loop_w

 ! mid-points
 dD = (floe_sz_max-floe_sz_min)/(res+1)
 do loop_w=1,res+1
  dum_D_vec(loop_w) = floe_sz_min + (loop_w-1)*dD
 end do
 do loop_w=1,res
  D_vec(loop_w) = (dum_D_vec(loop_w+1)+dum_D_vec(loop_w))/2d0
 end do

 ! FSD params

 beta0 = 1/(floe_sz_min**(-gamma0) - floe_sz_max**(-gamma1))
 beta1 = floe_sz_crit**gamma1

 if (floe_sz_max.ge.floe_sz_crit) then
  P0 = 1d0-0.05*((floe_sz_max/floe_sz_crit)**gamma1)
 else
  P0=0d0
 end  if

 if (P0.ge.1.0) then
  P0=1d0
 elseif (P0.le.0.0) then
  P0=0d0
 end if

 do loop_w=1,res
  if (ATTEN_MODEL.eq.0) then
   alpha_vec(loop_w) = fn_apxAttn(D_vec(loop_w),lam,1d0)
  elseif (ATTEN_MODEL.eq.1) then
   alpha_vec(loop_w) = fn_Attn_WIM_v1(dum_om, hice, dum_idl)
  end if
  p_vec(loop_w) = fn_splitPDF(D_vec(loop_w),beta0,beta1,P0)
 end do

 fn_IntAttn = dD*conc*dot_product(alpha_vec,p_vec/D_vec)

 end function fn_IntAttn

!=======================================================================

 !!!!!!!!!!!!!!!!!!
 !!! fn_proplng !!!
 !!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_proplng  -
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_proplng(dum_L)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
 real(kind=8), intent(in)     :: dum_L
!
! !OUTPUT PARAMETERS
!
 real(kind=8)                 :: fn_proplng
!
!EOP
!

 real(kind=8)                 :: dum_D1, dum_Es
 real(kind=8), dimension(nw)  :: alpha, S_attn
 integer                      :: loop_w

 L = dum_L

 dum_D1 = zero(i0,i1,toli,toli,fn_wlng)

 do loop_w=1,nw
  if (ATTEN_METH.eq.1) then
   alpha(loop_w) = fn_IntAttn(dum_D1, lam_ice(loop_w), om(loop_w))
  else
   alpha(loop_w) = fn_AvAttn(dum_D1, lam_ice(loop_w), om(loop_w))
  end if
  S_attn(loop_w) = S_init(loop_w)*exp(-alpha(loop_w)*L)
 end do

 call sub_StrainSpec(S_attn, dum_Es)

 fn_proplng = dum_Es - epsc*sqrt(-2/log(Pc))

 end function fn_proplng

 !=======================================================================

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! fn_proplng_uncoupled !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_proplng_uncoupled
!
! !DESCRIPTION:
!
!  fn_proplng for problem when attn/floe length are uncoupled
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_proplng_uncoupled(dum_L)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
 real(kind=8), intent(in)     :: dum_L
!
! !OUTPUT PARAMETERS
!
 real(kind=8)                 :: fn_proplng_uncoupled
!
!EOP
!

 real(kind=8)                 :: dum_Es
 real(kind=8), dimension(nw)  :: alpha, S_attn
 integer                      :: loop_w

 L = dum_L

 do loop_w=1,nw
  alpha(loop_w)  = conc*fn_Attn_MBK(om(loop_w))/0.7d0
  S_attn(loop_w) = S_init(loop_w)*exp(-alpha(loop_w)*L)
 end do

 call sub_StrainSpec(S_attn, dum_Es)

 fn_proplng_uncoupled = dum_Es - epsc*sqrt(-2/log(Pc))

 end function fn_proplng_uncoupled

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

 function fn_SpecMoment(dum_S, nw_in, nth_in, om_in, th_in, mom, idl)

!
! !USES:
!
!
! !INPUT PARAMETERS:
!
 integer, intent(in)                                :: nw_in, nth_in, mom, idl
 real (kind=8), dimension(nw_in), intent(in)        :: om_in
 real (kind=8), dimension(nth_in), intent(in)       :: th_in
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

 if (nth_in.eq.1) then
  dth_local = 3d0         ! for Simpson's rule
 else
  dth_local = th_in(2)-th_in(1)
 end if

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

 dum_simp(1) = 1d0
 dum_simp(nw_in) = 1d0

 do loop_w=2,nw_in-1,2
  dum_simp(loop_w) = 4d0
 end do

 do loop_w=3,nw_in-1,2
  dum_simp(loop_w) = 2d0
 end do

 dum_simp_th(1) = 1d0
 dum_simp_th(nth_in) = 1d0

 do loop_th=2,nth_in-1,2
  dum_simp_th(loop_th) = 4d0
 end do

 do loop_th=3,nth_in-1,2
  dum_simp_th(loop_th) = 2d0
 end do

 do loop_w=1,nw_in
  do loop_th=1,nth_in
   wt_simp(loop_w+nw_in*(loop_th-1)) = dum_simp(loop_w)*dum_simp_th(loop_th)
  end do
 end do

 do loop_w=1,nw_in
  do loop_th=1,nth_in
   wt_int(loop_w+nw_in*(loop_th-1)) = (dom_local/3d0)*(dth_local/3d0)* &
                              wt_simp(loop_w+nw_in*(loop_th-1))
   dum_v(loop_w+nw_in*(loop_th-1))  = dum_S(loop_w+nw_in*(loop_th-1))* &
                              (om_in(loop_w)**mom)
  end do
 end do

! if (idl.ne.0) then
!  write(idl,*) '                         wt_int=', wt_int(1), '->', wt_int(nw_in*nth_in)
!  write(idl,*) '                         v=', dum_v(1), '->', dum_v(nw_in*nth_in)
! end if


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

 !!!!!!!!!!!!!!!!!!!!!!!
 !!! sub_WavelenSpec !!!
 !!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: sub_WavelenSpec
!
! !DESCRIPTION:
!
!  calculates a representative wavelength for the spectrum
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!

 subroutine sub_WavelenSpec(dum_S, wlng_crest)

!
! !USES:
!
!
! !INPUT PARAMETERS:
!

 real (kind=8), dimension(nw*nth), intent(in) :: dum_S

!
! !OUTPUT PARAMETERS
!

 real (kind=8), intent(out)                   :: wlng_crest

!
!EOP
!

 integer                                      :: loop_w, loop_th
 real (kind=8), dimension(nw)                 :: dum_simp
 real (kind=8), dimension(nth)                :: dum_simp_th
 real (kind=8), dimension(nw*nth)             :: dum_v, dum_u !F, om_lng
 real (kind=8), dimension(nw*nth)             :: wt_simp, wt_int
 real (kind=8)                                :: mom0, mom2, T_crit, om_crit

 dum_simp(1) = 1d0
 dum_simp(nw) = 1d0

 do loop_w=2,nw-1,2
  dum_simp(loop_w) = 4d0
 end do

 do loop_w=3,nw-1,2
  dum_simp(loop_w) = 2d0
 end do

 dum_simp_th(1) = 1d0
 dum_simp_th(nth) = 1d0

 do loop_th=2,nth-1,2
  dum_simp_th(loop_th) = 4d0
 end do

 do loop_th=3,nth-1,2
  dum_simp_th(loop_th) = 2d0
 end do

 do loop_w=1,nw
  do loop_th=1,nth
   wt_simp(loop_w+nw*(loop_th-1)) = dum_simp(loop_w)*dum_simp_th(loop_th)
  end do
 end do

 do loop_w=1,nw
  do loop_th=1,nth
   wt_int(loop_w+nw*(loop_th-1)) = (dom/3d0)*(dth/3d0)*wt_simp(loop_w+nw*(loop_th-1))
   !om_lng(loop_w+nw*(loop_th-1)) = om(loop_w)
   if (hice.eq.0d0) then
    !F(loop_w+nw*(loop_th-1)) = 1d0
    dum_u(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))
    dum_v(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))*(om(loop_w)**2)
   else
    if (Wsq_METH.eq.0) then
     !F(loop_w+nw*(loop_th-1)) = 1d0
     dum_u(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))
     dum_v(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))*(om(loop_w)**2)
    elseif (Wsq_METH.eq.1) then
     !F(loop_w+nw*(loop_th-1)) = k_ice(loop_w)/k_wtr(loop_w)
     dum_u(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))* &
      ((k_ice(loop_w)/k_wtr(loop_w))**2)
     dum_v(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))*(om(loop_w)**2)* &
      ((k_ice(loop_w)/k_wtr(loop_w))**2)
    end if
   end if
  end do
 end do

 mom0     = dot_product(wt_int,dum_u)
 mom2     = dot_product(wt_int,dum_v)

 om_crit = sqrt(mom2/mom0)

 !if (out_ct.eq.1) then
 ! print*, 'wt_int=', wt_int
 ! print*, 'S     =', dum_S
 ! print*, 'om    =', om
 ! print*, 'mom0,mom2,om_crit=', mom0,mom2,om_crit
 ! out_ct = out_ct+1
 !end if

 kappa = om_crit**2d0/gravity

 wlng_crest = 2d0*pi/zero(0d0,max(kappa/tanh(kappa),sqrt(sqrt(kappa*mass/fr))), &
                          toli,toli,fn_DispRel_ice_inf)

 end subroutine sub_WavelenSpec

!=======================================================================

 !!!!!!!!!!!!!!!!!!!!!!
 !!! sub_StrainSpec !!!
 !!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: sub_StrainSpec
!
! !DESCRIPTION:
!
!  calculates the strain imposed on the ice by the wave spectrum
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!

 subroutine sub_StrainSpec(dum_S, Es)

!
! !USES:
!
!
! !INPUT PARAMETERS:
!

 real (kind=8), dimension(nw), intent(in) :: dum_S

!
! !OUTPUT PARAMETERS
!

 real (kind=8), intent(out)               :: Es

!
!EOP
!

 integer                                  :: loop_w, loop_th
 real (kind=8), dimension(nw)             :: dum_simp
 real (kind=8), dimension(nth)            :: dum_simp_th
 real (kind=8), dimension(nw*nth)         :: dum_vec !, F
 real (kind=8), dimension(nw*nth)         :: wt_simp, wt_int
 real (kind=8)                            :: mom0_eps, T_crit

 dum_simp(1) = 1d0
 dum_simp(nw) = 1d0

 do loop_w=2,nw-1,2
  dum_simp(loop_w) = 4d0
 end do

 do loop_w=3,nw-1,2
  dum_simp(loop_w) = 2d0
 end do

 dum_simp_th(1) = 1d0
 dum_simp_th(nth) = 1d0

 do loop_th=2,nth-1,2
  dum_simp_th(loop_th) = 4d0
 end do

 do loop_th=3,nth-1,2
  dum_simp_th(loop_th) = 2d0
 end do

 do loop_w=1,nw
  do loop_th=1,nth
   wt_simp(loop_w+nw*(loop_th-1)) = dum_simp(loop_w)*dum_simp_th(loop_th)
  end do
 end do

! if (out_ct.eq.1) then
!  print*, 'dum_simp   ', dum_simp
!  print*, 'dum_simp_th', dum_simp_th
!  out_ct=out_ct+1
! end if

 do loop_w=1,nw
  do loop_th=1,nth
   wt_int(loop_w+nw*(loop_th-1))   = (dom/3d0)*(dth/3d0)*wt_simp(loop_w+nw*(loop_th-1))
   if (hice.eq.0d0) then
    !F(loop_w+nw*(loop_th-1)) = 1d0
    dum_vec(loop_w+nw*(loop_th-1)) = 0.25d0*(hice**2d0)*dum_S(loop_w+nw*(loop_th-1)) &
         *(k_ice(loop_w)**4d0)
   else
    if (Wsq_METH.eq.0) then
     !F(loop_w+nw*(loop_th-1)) = 1d0
     dum_vec(loop_w+nw*(loop_th-1)) = 0.25d0*(hice**2d0)*dum_S(loop_w+nw*(loop_th-1)) &
         *(k_ice(loop_w)**4d0)
    elseif (Wsq_METH.eq.1) then
     !F(loop_w+nw*(loop_th-1)) = k_ice(loop_w)/k_wtr(loop_w)
     dum_vec(loop_w+nw*(loop_th-1)) = 0.25d0*(hice**2d0)*dum_S(loop_w+nw*(loop_th-1)) &
         *(k_ice(loop_w)**4d0)*((k_ice(loop_w)/k_wtr(loop_w))**2d0)
    end if
   end if
  end do
 end do

 mom0_eps = dot_product(wt_int,dum_vec)

 Es = 2d0*sqrt(mom0_eps)

 end subroutine sub_StrainSpec

!=======================================================================

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!! fn_DispRel_ice_inf !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_DispRel_ice_inf  - infinite depth ice-coupled dispersion relation
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_DispRel_ice_inf(dum_k)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real (kind=8), intent(in) :: dum_k

!
! !OUTPUT PARAMETERS
!

 real (kind=8) :: fn_DispRel_ice_inf

!
!EOP
!

 fn_DispRel_ice_inf = (1d0-(mass*kappa)+(fr*(dum_k**4d0)))*dum_k - kappa

 end function fn_DispRel_ice_inf

!=======================================================================

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! fn_wlng = lam/2 - D1 !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: fn_wlng  - calculate wavelength
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function fn_wlng(dum_D1)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 real(kind=8), intent(in)          :: dum_D1

!
! !OUTPUT PARAMETERS
!

 real(kind=8) :: fn_wlng

!
!EOP
!
 real(kind=8), dimension(nw)       :: alpha
 real(kind=8), dimension(nw*nth)   :: S_attn
 integer                           :: loop_w, loop_th
 real(kind=8)                      :: lam_attn

 do loop_w=1,nw
  if (ATTEN_METH.eq.1) then
   alpha(loop_w) = fn_IntAttn(dum_D1, lam_ice(loop_w), om(loop_w))
  else
   alpha(loop_w) = fn_AvAttn(dum_D1, lam_ice(loop_w), om(loop_w))
  end if
  do loop_th=1,nth
   S_attn(loop_w+nw*(loop_th-1)) = S_init(loop_w+nw*(loop_th-1))* &
     exp(-alpha(loop_w)*cos(th(loop_th))*L)
  end do
 end do

 call sub_WavelenSpec(S_attn, lam_attn)

 fn_wlng = lam_attn/2d0 - dum_D1

 end function fn_wlng

 !=======================================================================

 !!!!!!!!!!!!!!!!!!!!!!!!!
 !!! SDF_Bretschneider !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!BOP
!
! !ROUTINE: SDF_Bretschneider  - Bretschneider wave spectrum
!
! !DESCRIPTION:
!
!  blah blah blah
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  function SDF_Bretschneider(omega,moment_no,Hs,Tm)
!
! !USES:
!
! !INPUT PARAMETERS:
!

 integer, intent(in)       :: moment_no
 real (kind=8), intent(in) :: omega
 real (kind=8), intent(in) :: Hs       ! sig wave height
 real (kind=8), intent(in) :: Tm       ! peak period


!
! !OUTPUT PARAMETERS
!

 real (kind=8) :: SDF_Bretschneider

!
!EOP
!

 real (kind=8)             :: om_m, tau
 real (kind=8)             :: f1, f2, f3

 om_m = 2d0*pi/Tm
 tau  = 2d0*pi/omega

 f1 = (5d0/16d0)*(Hs**2)*(om_m**4)
 f2 = omega**(moment_no-5)
 f3 = exp(-1.25d0*((tau/Tm)**4))

 SDF_Bretschneider = f1*f2*f3

 end function SDF_Bretschneider

!=======================================================================

 end module m_waveice

!=======================================================================
