module HAL_N_GAL_BARYONS

  ! module which deals with exchanges between halos and galaxies

  use BOOKKEEP_BARYONS

  public

  type halo_return
    real(kind=8)    :: mhot  ! mass of gas to be thrown back to the hot halo phase
    real(kind=8)    :: mhotz ! metals ...
    real(kind=8)    :: mout  ! mass of gas to be thrown out the halo
    real(kind=8)    :: moutz ! metals ... 
  end type halo_return


contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!*****************************************************************************************************************  
  subroutine evolve_halo(hbm,ham,i_t,st)

    ! subroutine which deals with galaxy-galaxy mergers during time interval [i_t(1),i_t(2)].
    ! Computes orbital radii as well.    
    ! make a temporary list of ghost galaxies (descendents of hbm) and evolve it. Then make remaining gals the progs of ham gals.

    implicit none

    type(halo)                    :: hbm          ! halo before galaxy mergers
    type(halo)                    :: ham          ! halo after galaxy mergers
    real(kind=8),intent(in)       :: i_t(2)       ! i_t(1) and i_t(2) are cosmic times in between which we evolve the galaxies
    real(kind=8)                  :: t1           
    integer(kind=4),intent(in)    :: st           ! timestep to which we evolve the galaxies (i.e. ham%my_timestep)
    type(galaxy),allocatable      :: list_gals(:) ! list of galaxies in the halo
    integer(kind=4)               :: nbgals       ! number of galaxies in halo
    integer(kind=4)               :: i
    integer(kind=4)               :: nbgal_in
    integer(kind=4)               :: renum

    nbgals    = hbm%datas%nbgal           ! necessarily > 0 (routine not called otherwise)
    nbgal_in  = nbgals

   ! define a list of ghost galaxies which are progenitors of galaxies of ham 
    allocate(list_gals(nbgals))
    if (i_t(1) == tsno(st-1)%age_univ) then
       renum   = hbm%my_number
#ifdef RENUMBERING
       call renumbering(renum,st-1)
#endif
       ! make each galaxy in hbm point to itself (as a prog) before copy to list_gal
       ! -> gals of list_gal, when ending up in ham, will point to proper progs.... 
       if (tsno(st-1)%liste_halos(renum)%datas%nbgal > 0) then
          do i = 1,nbgals
             if (hbm%datas%liste_galaxies(i)%my_progs%nb_prog > 0) then
                deallocate(hbm%datas%liste_galaxies(i)%my_progs%hno,hbm%datas%liste_galaxies(i)%my_progs%gno)
             end if
             hbm%datas%liste_galaxies(i)%my_progs%nb_prog    = 1
             allocate(hbm%datas%liste_galaxies(i)%my_progs%hno(1),hbm%datas%liste_galaxies(i)%my_progs%gno(1))
             hbm%datas%liste_galaxies(i)%my_progs%hno(1)     = hbm%my_number
             hbm%datas%liste_galaxies(i)%my_progs%gno(1)     = hbm%datas%liste_galaxies(i)%my_number
             ! also reset SFR counters ...
             call reset_sfr(hbm%datas%liste_galaxies(i))
             call clear_gal_bookkeep(hbm%datas%liste_galaxies(i)%gbk)
          end do
       end if
    end if
    do i = 1,nbgals
       call init_galaxy(list_gals(i))
       call copy_gal_props(list_gals(i),hbm%datas%liste_galaxies(i))
    end do
    
    

    ! deal with mergers and passive evolution
    t1 = i_t(1)
#ifndef NO_MERGERS
    call compute_mergers(hbm,list_gals,t1,i_t,nbgals,st) 
#else
    call passive_evolution(hbm,list_gals,t1,i_t(2),nbgals,st)
#endif

    ! now update final halo and its galaxy list
    call alloc_hlist(ham,nbgals)
    do i=1,nbgals
       call copy_gal_props(ham%datas%liste_galaxies(i),list_gals(i))
    end do

    ! deallocate temporary galaxy list
    do i = 1,nbgal_in
       call clear_galaxy(list_gals(i))
    end do
    deallocate(list_gals)

  contains 

!--------------------------------------------------------------------------------------------------
    recursive subroutine compute_mergers(h,list_gals,t1,i_t,nbgals,st)

   ! routine which performs the galaxy mergers and builds the galaxy tree from i_t(1) to i_t(2)

      implicit none

      type(halo)          :: h                 ! halo to be evolved
      integer(kind=4)     :: nbgals            ! number of galaxies in the halo at time t1
      type (galaxy)       :: list_gals(nbgals) ! list of galaxies to evolve
      real(kind=8)        :: t1                ! time at which the evolution starts for current call to compute_mergers
      real(kind=8)        :: i_t(2)            ! bounds to time interval covered by compute_mergers (with i_t(1)<i_t(2))
      integer(kind=4)     :: st                ! timestep. Used to update merger stats in tsno
      !integer(kind=4)     :: i
      real(kind=8)        :: t2
      logical(kind=4)     :: there_is_a_sat_merger
      integer(kind=4)     :: i_central
      real(kind=8)        :: df_time, tmerge

      if (nbgals > 2) then  ! possible sat-sat mergers

         ! Compute times needed for each galaxy to fall to the center of halo due to dynamical friction.
         ! df_time and i_central are the shortest time to reach center and the index of the galaxy to which it
         ! corresponds.
         call time_to_fall(h,list_gals,df_time,i_central,nbgals)
       
         ! Check if any galaxy has time to fall to the center during timestep
         if (df_time .le. i_t(2) - t1) then

            ! define time interval between t1 and time of DF merger as if there was no sat-sat mergers
            t2 = min(t1 + df_time,i_t(2))            
            ! look for (and perform) sat-sat mergers during this time interval [t1,t2]
            call sat_sat_merger(h,list_gals,t1,t2,tmerge,there_is_a_sat_merger,nbgals,st)
            if(there_is_a_sat_merger) then 
               ! there has been a sat-sat merger before first galaxy reached center so change t1 to tmerge and 
               ! update number of gals in list. Then, call compute_mergers for the new list.
               t1     = tmerge
               tsno(st)%n_ss = tsno(st)%n_ss + 1  ! update timestep stats
               if (t1 < i_t(2)) then 
                  call compute_mergers(h,list_gals,t1,i_t,nbgals,st)
               end if
            else
               ! there is no sat-sat merger before t2, so evolve all galaxies up to t2 and merge galaxy which reaches center at t2
               ! Then, update t1 to t2 and nbgals, and re-call compute_mergers if t2 < i_t(2)
               call central_merger(h,list_gals,t1,t2,i_central,nbgals,st)
               t1     = t2
               if (t1 < i_t(2)) then 
                  call compute_mergers(h,list_gals,t1,i_t,nbgals,st)
               end if
            end if

         else   ! no DF merger --> check for sat-sat mergers between t1 and i_t(2)
            
            call sat_sat_merger(h,list_gals,t1,i_t(2),tmerge,there_is_a_sat_merger,nbgals,st)
            if (there_is_a_sat_merger) then
               ! there is a sat-sat merger --> list of galaxies has been updated in sat_sat_merger, so 
               ! redefine t1 and call compute_merger again for evolved halo
               t1     = tmerge
               tsno(st)%n_ss = tsno(st)%n_ss + 1
               if (t1 < i_t(2)) then 
                  call compute_mergers(h,list_gals,t1,i_t,nbgals,st)
               end if
            else
               ! passive evolution of all galaxies from t1 to i_t(2)
               if (t1 < i_t(2)) then 
                  call passive_evolution(h,list_gals,t1,i_t(2),nbgals,st)
               end if
            end if

         end if

      else if (nbgals == 2) then ! only dynamical friction merger possible
         
         ! Compute times needed for each galaxy to fall to the center of halo due to dynamical friction.
         ! df_time and i_central are the shortest time to reach center and the index of the galaxy to which it
         ! corresponds.

         call time_to_fall(h,list_gals,df_time,i_central,nbgals)

         ! Check if any galaxy has time to fall to the center during timestep
         if (df_time .le. i_t(2) - t1) then
            ! define time interval between t1 and time of DF merger as if there was no sat-sat mergers
            t2 = min(t1 + df_time,i_t(2))            
            call central_merger(h,list_gals,t1,t2,i_central,nbgals,st)
            t1     = t2
            if (t1 < i_t(2)) then 
               call compute_mergers(h,list_gals,t1,i_t,nbgals,st)
            end if
         else   ! no DF merger --> passive evolution between t1 and i_t(2)
            call passive_evolution(h,list_gals,t1,i_t(2),nbgals,st)
         end if
         
      else  ! only one galaxy left in the halo so evolve it passively up to i_t(2)

         call passive_evolution(h,list_gals,t1,i_t(2),nbgals,st)
 
      end if

      return

    end subroutine compute_mergers

!--------------------------------------------------------------------------------------------------
    subroutine time_to_fall(h,list_gals,time,i_central,nbgals)

    ! compute time to fall, i.e. the time after which delta_r due to dyn fric is equal to r_orbit
    ! NB : also redefines core mass of satellite

      implicit none

      type(halo),intent(inout)    :: h                 ! halo in which we merge the galaxies
      real(kind=8),intent(out)    :: time              ! time needed for galaxy to fall to center due to DF
      integer(kind=4),intent(out) :: i_central         ! index of galaxy which reaches center first
      integer(kind=4)             :: nbgals            ! number of gals in halo
      type(galaxy)                :: list_gals(nbgals) ! list of galaxies to evolve

      integer(kind=4),parameter   :: nr   = 100        ! number of bins needed to integrate dynamical friction diff. eq.
      real(kind=8),parameter      :: fact = 2.655d0    ! this is 1/(2*0.428*1.023e-3*Gravconst) (see galaxy_orbits)
      real(kind=8)                :: r                 ! initial orbital radius 
      real(kind=8)                :: msat              ! mass of the satellite (baryons + DM)
      real(kind=8)                :: a                 ! Coulomb logarithm
      integer(kind=4)             :: i,j               ! loop variable
      real(kind=8)                :: t                 ! local 'time'
      real(kind=8)                :: rgal,mgal,cgal    ! radius, mass and core_mass of gal
      real(kind=8)                :: delta_r,vc 
      real(kind=8)                :: time2  
      real(kind=8),parameter      :: gravconst_galics = 0.00044965535 ! Grav. const. in Galics units = Mpc^3.(1.d11Msun)^-1.Gyr^2

      ! initialize time
      time      = 100.d0 ! just a number > age of the universe
      i_central = 0

      do i = 1, nbgals
         ! define initial orbital radius
         r       = list_gals(i)%r 
         if (r > 0.0d0) then
            rgal    = max(disc_params%s_to_m * list_gals(i)%disc%rgal, bulge_params%s_to_m * list_gals(i)%bulge%rgal, &
                 & bulge_params%s_to_m * list_gals(i)%burst%rgal)
            mgal    = total_mass(list_gals(i))
            cgal    = list_gals(i)%core%mass
            t       = 0.0d0
            delta_r = (r - rgal) / dble(nr)

!tibo comments integration of df_time calculation cuz tdyn replaces df_time

!!$            do j = 1,nr
!!$               ! define mass of satellite halo after core stripping (only DM)
!!$               call compute_mass_satellite(r,h,msat,list_gals(i))
!!$               call compute_circular_velocity(r,vc,h)
!!$
!!$               cgal = min(cgal,msat)
!!$               ! add baryon contribution to sat mass 
!!$               msat = cgal + mgal
!!$               
!!$               if (msat == 0.0d0) then ! of course, this never happens ... 
!!$                  write(errunit,*) 'crash ahead ... ',list_gals(i)%core%r_0, list_gals(i)%core%rho_0
!!$                  write(errunit,*) cgal,list_gals(i)%core%mass
!!$               end if
!!$
!!$               a    = log(1.0d0 + h%mfof / msat)   ! Coulomb logarithm
!!$               ! simple trapezoidal integration 
!!$               t    = t + delta_r * (r - rgal) * vc / msat / a 
!!$               r    = r - delta_r
!!$            end do
!!$            t    = t * fact

!tibo : change merging time: dyn_friction_time => dyn_halo_time (faux: virer le 8)
            t = sqrt(h%datas%rvir**3. / (8.0d0 * gravconst_galics * h%datas%mvir))

            !tibo
            !t = 0.001d0 ! Gyr = 1 Myr = very short
            
            if (time > t) then  ! time is time taken for fastest galaxy to fall to the center
               time      = t
               !tibo
               i_central = i !i
               
            end if
         end if
      end do

      return

    end subroutine time_to_fall

!--------------------------------------------------------------------------------------------------
    subroutine passive_evolution(h,list_gals,t1,t2,nbgals,st)

    ! Evolve passively (i.e. no mergers) all the nbgals galaxies of list_gals between t1 and t2, inside halo h.
    ! This includes dynamical friction decrement of orbital radii.

      implicit none

      type(halo)       :: h                 ! halo harbouring galaxies of list_gals
      integer(kind=4)  :: nbgals            ! number of gals in halo at time t1
      type(galaxy)     :: list_gals(nbgals) ! list of galaxies to evolve
      real(kind=8)     :: t1, t2            ! evolve gals from t1 to t2
      integer(kind=4)  :: st                ! timestep
      integer(kind=4)  :: i
      real(kind=8)     :: delta_t

    
      delta_t = t2 - t1
      do i = 1, nbgals
         ! compute star formation and feedback + associated bookkeeping from time t1 to time t2
         call star_feed(list_gals(i),t1,t2,h,st)
         list_gals(i)%tgal = t2 ! record that gal has been evolved up to t2
         ! compute disturbance on morphology due to a possible recent merger 
         call calculate_disturbance(list_gals(i),delta_t)
      end do      

#ifndef NO_MERGERS
      ! compute dynamical friction 
      call galaxy_orbits(list_gals,t1,t2,h,nbgals)
#endif

      return

    end subroutine passive_evolution

!--------------------------------------------------------------------------------------------------
    subroutine sat_sat_merger(h,list_gals,t1,t2,tmerge,there_is_a_sat_merger,nbgals,st)

    ! subroutine which checks for sat-sat mergers in time interval [t1,t2]. In case there are mergers, 
    ! one only is performed, and all non-merging galaxies are evolved from t1 to tmerge. 
    ! NB : collision with central galaxy is not allowed.

      implicit none      

      type(halo)                  :: h
      integer(kind=4)             :: nbgals                ! nb of gals in halo at t1 (returns nbgals-1 if there is a merger)
      type(galaxy)                :: list_gals(nbgals)     ! list of galaxies to evolve
      real(kind=8),intent(in)     :: t1, t2                ! time interval during which we look for sat-sat mergers
      real(kind=8)                :: tmerge                ! time at which merger occurs
      logical(kind=4)             :: there_is_a_sat_merger ! returns true if a merger was found
      integer(kind=4)             :: st                    ! timestep
      !integer(kind=4)             :: i,j
      real(kind=8)                :: delta_t               ! t2 - t1
      real(kind=8)                :: m, r, v, p_merge, res
      integer(kind=4)             :: nbmerge
      integer(kind=4),allocatable :: merge_index(:)        ! index of gals willing to merge
      integer(kind=4)             :: gals_to_merge(2)      ! index of gals to merge (allways binary mergers)
      integer(kind=4)             :: maxiteration          ! stop infinite while loops ... 

      ! initialize a few things
      there_is_a_sat_merger = .false.
      delta_t               = (t2 - t1) * psi * dble(nbgals - 1)  ! this is more useful than delta_t later...
      nbmerge               = 0
      gals_to_merge(:)      = 0
      allocate(merge_index(nbgals))

      ! spot galaxies that may merge during delta_t
      do i = 1, nbgals
         if (list_gals(i)%r > 0.0d0) then
            m = total_mass(list_gals(i))
            if (m > 0.0d0) then ! make v and r mass-weighted averages.  
               v = (list_gals(i)%disc%mgal    * list_gals(i)%disc%speed  + & 
                    & list_gals(i)%bulge%mgal * list_gals(i)%bulge%speed + & 
                    & list_gals(i)%burst%mgal * list_gals(i)%burst%speed) / m
               r = (list_gals(i)%disc%mgal * list_gals(i)%disc%rgal  * disc_params%s_to_m  + & 
                    & list_gals(i)%bulge%mgal * list_gals(i)%bulge%rgal * bulge_params%s_to_m + & 
                    & list_gals(i)%burst%mgal * list_gals(i)%burst%rgal * bulge_params%s_to_m) / m
               p_merge = delta_t * (r / h%rfof)**3 * (v / sqrt(gravconst*h%mfof/h%rfof))**3 * (v / r)
               p_merge = 1.d0 - exp(-p_merge)
               call ran1(iseed,res)
               if (res < p_merge) then  ! there is a merger for this galaxy
                  nbmerge              = nbmerge + 1
                  merge_index(nbmerge) = i
               end if
            else 
               if (accretion_remaining_mass(list_gals(i)%wind,0.0d0) &
                    & + accretion_remaining_mass(list_gals(i)%acc,0.0d0) <= 0.0d0) then 
                  write(errunit,*) '> fatal error: no mass for gal in ss_merging '
                  stop
               end if
            end if
         end if
      end do

      if (nbmerge > 1) then  ! there are more than two galaxies willing to merge

         ! --> select a random pair of galaxies willing to merge and perform their merger
         there_is_a_sat_merger = .true.
         call ran1(iseed,res)
         i = min(nint(res * nbmerge + 0.5d0),nbmerge)
         gals_to_merge(1) = merge_index(i)
         gals_to_merge(2) = gals_to_merge(1)
         maxiteration = 0
         do while (gals_to_merge(1) == gals_to_merge(2)) 
            maxiteration = maxiteration + 1
            call ran1(iseed,res)
            i = min(nint(res * nbmerge + 0.5d0),nbmerge)
            gals_to_merge(2) = merge_index(i)   ! cannot be central
            if (maxiteration > 100) then 
               ! cut the crap... 
               gals_to_merge(1) = merge_index(1)
               gals_to_merge(2) = merge_index(2)
               exit
            end if
         end do
         ! draw a random time for merger between t1 and t2 
         call ran1(iseed,res)
         tmerge           = t1 + res * (t2 - t1)

         ! jeje : i suppress the following as it is really really useless ...
         ! (and contains an infinite while loop which is really really annoying !)
!!$      else if (nbmerge == 1) then ! only one galaxy willing to merge
!!$
!!$         ! --> half a chance to merge with any other galaxy. Select other galaxy randomly if merger
!!$         call ran1(iseed,res)
!!$         if (res < 0.5) then ! merger
!!$            there_is_a_sat_merger = .true.
!!$            gals_to_merge(1) = merge_index(1)
!!$            gals_to_merge(2) = gals_to_merge(1)
!!$            do while (gals_to_merge(1) == gals_to_merge(2)) 
!!$               if (st == 65) print*,'in the while sat-sat 2'
!!$               call ran1(iseed,res)
!!$               gals_to_merge(2) = min(nint(res * nbgals + 0.5),nbgals)
!!$               if (list_gals(gals_to_merge(2))%r == 0.0) then ! prevent mergers with central galaxy
!!$                  gals_to_merge(2) = gals_to_merge(1)
!!$               end if
!!$            end do
!!$            ! draw a random time for merger between t1 and t2 
!!$            call ran1(iseed,res)
!!$            tmerge           = t1 + res * (t2 - t1)
!!$         end if

      end if                 ! otherwise, no merger

      deallocate(merge_index)

      ! if there is a merger, evolve all galaxies passively to tmerge and then merge the two chosen ones
      if (there_is_a_sat_merger) then 

         call passive_evolution(h,list_gals,t1,tmerge,nbgals,st)

         call perform_galaxy_merger(h,list_gals,gals_to_merge(1),gals_to_merge(2),nbgals,there_is_a_sat_merger,st)

      end if ! otherwise, return

      return

    end subroutine sat_sat_merger

!--------------------------------------------------------------------------------------------------
    subroutine central_merger(h,list_gals,t1,t2,i_central,nbgals,st)

    ! subroutine which performs the merger of a galaxy with the central galaxy (if there is a central galaxy)
    ! and evolves all other galaxies from t1 to t2 (t2 = time of merger)

      implicit none

      type(halo)                  :: h
      integer(kind=4)             :: nbgals                ! nb of gals in halo at t1 (returns nbgals-1 if there is a merger)
      type(galaxy)                :: list_gals(nbgals)     ! list of galaxies to evolve
      real(kind=8),intent(in)     :: t1, t2                ! time interval during which we look for sat-sat mergers
      integer(kind=4)             :: i_central             ! index of galaxy to merge with central galaxy 
      integer(kind=4)             :: st
      integer(kind=4)             :: ind_center            ! index of central galaxy (at t1)
      real(kind=8)                :: r 
      logical(kind=4)             :: ss            

      ! find out if there is a central galaxy
      r          = list_gals(1)%r
      ind_center = 1
      do while (r > 0 .and. ind_center < nbgals)
         ind_center = ind_center + 1
         if (list_gals(ind_center)%r == 0.0d0) then
            r = 0.0d0
         end if
      end do

      if (r == 0.0d0) then  ! there is a central galaxy : list_gals(ind_center)

         ! --> have to perform a merger between gals i_central and ind_center
         ! but first, evolve all galaxies passively to t2 
         if (t1 < t2) then 
            call passive_evolution(h,list_gals,t1,t2,nbgals,st)
         end if
         ss = .false. ! this is NOT a sat-sat merger
        
         call perform_galaxy_merger(h,list_gals,ind_center,i_central,nbgals,ss,st)
         tsno(st)%n_dynfric = tsno(st)%n_dynfric + 1

      else                ! there is no central galaxy

         ! --> just place galaxy which fell to the center at time t2 to the center. Evolve all gals to t2, and return.
         if (t1 < t2) then 
            call passive_evolution(h,list_gals,t1,t2,nbgals,st)
         end if
         ! make sure that galaxy i_central is at the center 
         list_gals(i_central)%r = 0.0d0

      end if

      return

    end subroutine central_merger

!--------------------------------------------------------------------------------------------------
    subroutine galaxy_orbits(list_gals,t1,t2,h,nbgals)

    ! using the dynamical friction formula of Primack and Somerville (eq 9), this routine
    ! computes the change in satellite galaxies of h1 distances from centre, during delta_t.
    ! (see GALICS I sections 5.1,5.2,5.3 for prescriptions)
    ! NB: 1.023e-3 factor in expression is for units conversion (see GALICS I, secs 5.1, 5.2, 5.3)
    !     factor 2 is to take into account non-circular orbits 

      implicit none

      type(halo)                 :: h                 ! halo in which the galaxies are
      real(kind=8),intent(in)    :: t1, t2            ! compute DF from t1 to t2
      integer(kind=4)            :: nbgals            ! number of gals in the halo at t1
      type(galaxy)               :: list_gals(nbgals) ! list of galaxies in the halo
      real(kind=8)               :: r                 ! local orbital radius
      real(kind=8)               :: msat              ! mass of satellite
      real(kind=8)               :: a                 ! Coulomb log
      real(kind=8)               :: vc                ! circular vel of halo at radius r
      integer(kind=4)            :: i,j
      real(kind=8)               :: delta_t,dt        ! t2 - t1 = time during which we compute DF
      integer(kind=4),parameter  :: nt   = 100        ! number of bins for integration of DF diff. eq.
      real(kind=8),parameter     :: fact = 0.377d0    ! = 2. * 0.428 * 1.023e-3 * gravconst
      real(kind=8)               :: rgal,cgal,mgal    ! galaxy properties             

      delta_t = t2 - t1

      do i = 1, nbgals               ! loop over all the galaxies in halo at time t1

         r      = list_gals(i)%r     ! define local short-name for orbital radius.

         if (r > 0.0d0) then           ! galaxy is not central
            ! galaxy properties
            rgal = max(disc_params%s_to_m * list_gals(i)%disc%rgal, bulge_params%s_to_m * list_gals(i)%bulge%rgal, &
                 & bulge_params%s_to_m * list_gals(i)%burst%rgal)
            mgal = total_mass(list_gals(i))
            cgal = list_gals(i)%core%mass
            dt   = delta_t / dble(nt)
            do j = 1,nt
               ! define mass of satellite halo after core stripping
               call compute_mass_satellite(r,h,msat,list_gals(i))
               ! msat is the stripped core mass
               ! %core%mass is DM core associated with the galaxy and is whatever it was previously  
               ! the core mass we want to use is therefore min(msat,core%mass)
               cgal = min(cgal,msat)
               ! add baryons to DM core
               msat = cgal + mgal 
               a    = h%mfof / msat ! dimensionless                
               ! now express variation of orbital radius due to dynamical friction (equs. 5.2 and 5.3 of 
               ! GALICS I)
               call compute_circular_velocity(r,vc,h)
               r = r - fact * msat * log(1.0d0+a) / r / vc * dt
               if (r < rgal) then
                  if (nbgals == 1) then
                     r    = 0.0d0
                     rgal = 0.0d0
                  end if
                  exit
               end if
            end do
            ! update global stats
            tsno(st)%tau_count       = tsno(st)%tau_count+1                
            tsno(st)%average_tau     = tsno(st)%average_tau+(((list_gals(i)%r-max(r,rgal))/delta_t)/ &
                 & sqrt(gravconst*h%mfof/h%rfof)) * 977.9d0
            ! update galaxy field
            list_gals(i)%r           = max(r,rgal) 
            list_gals(i)%core%mass   = cgal
         end if

      end do

      return

    end subroutine galaxy_orbits

!--------------------------------------------------------------------------------------------------
    subroutine compute_mass_satellite(r,h,msat,gal)

    ! computes the mass msat of a satellite falling in a host halo h after its DM core has been stripped
    ! by solving rho_sat(r_t) = rho_halo(r) (r_t is the truncation radius and r the orbit radii of the sat)
    ! and integrating from 0 to r_t (R1_det).
    ! assumes the same density profile for both host and satellite halo.

      implicit none 

      type (halo)                :: h
      real(kind=8),intent(in)    :: r
      type(galaxy)               :: gal
      real(kind=8)               :: msat
      real(kind=8)               :: rho_halo,rho_0,r_0
      real(kind=8)               :: rho_0_sat,r_0_sat,A_det,B_det,C_det,R1_det,R2_det

      rho_0      = h%halo_profile%rho_0
      r_0        = h%halo_profile%r_c
      rho_0_sat  = gal%core%rho_0
      r_0_sat    = gal%core%r_0

      ! have to check first if profile is defined: it is not if halo is not at least partially virialized
      if (h%datas%rvir > 0.0d0) then 

         if (profile == 'TSIS') then 

            ! rho_0 is actually rho_vir and r_0 is r_vir as the density profile is singular in r=0
            ! we have to use the same "trick" as for the potential i.e. we truncate the density profile @ 0.1 kpc
            ! idem for the satellite: rho_0_sat is rho_vir_sat and r_0_sat is r_vir_sat in the TSIS case
            if (r > min_size) then 
               R1_det = sqrt(rho_0_sat/rho_0)*r_0_sat/r_0*r
            else
               R1_det = sqrt(rho_0_sat/rho_0)*r_0_sat/r_0*min_size
            endif
            ! satellite is truncated at its virial radius in any case ...
            R1_det = min(R1_det,r_0_sat)
            msat   = 4.d0*pi*rho_0_sat*r_0_sat**2*R1_det

         else if (profile == 'NSTIS') then

            rho_halo = rho_0*(ais/(a2is+(r/r_0)**2) - bis/(b2is+(r/r_0)**2))
            ! compute core stripping :  solve the second degree equation to obtain r_t (R1_det) if it exists ...
            ! NB: following if statement handles numerical precision ...
            if ((rho_0_sat/rho_halo-1.d0) > 1e-3) then ! satellite is dense enough: strip ext regions ONLY  
               A_det = rho_halo / rho_0_sat
               B_det = A_det*(a2is+b2is) - ais + bis
               C_det = A_det*a2is*b2is - ais*b2is + bis*a2is
               if (B_det*B_det - 4.d0*A_det*C_det < 0.0d0) then
                  write(errunit,*) 'first det problem in galaxy_orbits',B_det*B_det-4*A_det*C_det
                  stop
               else
                  R1_det      = (-B_det + sqrt(B_det*B_det - 4.d0*A_det*C_det))/(2.d0*A_det) 
                  R2_det      = (-B_det - sqrt(B_det*B_det - 4.d0*A_det*C_det))/(2.d0*A_det)
                  if (R1_det < 0.0d0) then 
                     if (R2_det < 0.0d0) then 
                        write(errunit,*) 'second det problem in galaxy_orbits',R1_det,R2_det
                        write(errunit,*) A_det,B_det,C_det
                        write(errunit,*) h%my_number,r,rho_halo,rho_0_sat
                        stop 
                     else
                        R1_det = R2_det
                     endif
                  else
                     if (R2_det > 0.0d0) then
                        if (R1_det == 0.0d0) then
                           R1_det = R2_det
                        else
                           R1_det = min(R1_det,R2_det)
                        endif
                     endif
                  endif
                  R1_det = sqrt(R1_det)
               endif
               msat   = 4*pi*rho_0_sat*r_0_sat**3*((ais-bis)*R1_det - sqrt(a2is)*ais*atan(R1_det/sqrt(a2is)) + &
                    & sqrt(b2is)*bis*atan(R1_det/sqrt(b2is)))
            else ! satellite is not dense enough: must be stripped ENTIRELY
               msat   = 0.0d0
            endif

         else if (profile == 'NFW') then

            write(errunit,*) '> error in compute_halo_density: profile ',profile,' not implemented yet'
            stop

         endif
         
      else ! halo is not virialized at all --> assume DM profile is homogeneous sphere for host halo

         rho_halo = h%mfof/(4.d0/3.d0*pi*h%rfof**3)
         ! check what is the profile of the satellite
         if (profile == 'TSIS') then 

            ! rho_0_sat is rho_vir_sat and r_0_sat is r_vir_sat in the TSIS case
            R1_det = sqrt(rho_0_sat/rho_halo)*r_0_sat
            ! satellite is truncated at its virial radius in any case ...
            R1_det = min(R1_det,r_0_sat)
            msat   = 4.d0*pi*rho_0_sat*r_0_sat**2*R1_det
            
         else

            write(errunit,*) '> error in compute_halo_density: profile not implemented yet'
            stop

         endif
            
      endif

      if (msat < 0.0d0) then
         write(errunit,*) '> Problem in compute_mass_satellite : msat < 0'
         write(errunit,*) '> ',msat,i,h%my_number,h%my_timestep
         stop
      end if

      return

    end subroutine compute_mass_satellite

!--------------------------------------------------------------------------------------------------
    subroutine perform_galaxy_merger(h,list_gals,i1,i2,nbgals,ss,st)

    ! routine which indeed performs a merger between two galaxies among a list, and updates all family links.

      implicit none

      integer(kind=4)              :: nbgals            ! nb of gals in halo (or list), before merger
      type(halo)                   :: h                 ! halo in which merger occurs
      type(galaxy)                 :: list_gals(nbgals) ! list of galaxies in the halo
      integer(kind=4)              :: i1, i2            ! indices of the two gals to merge
      logical(kind=4)              :: ss                ! if true, sat-sat merger
      integer(kind=4)              :: st                ! timestep 
      integer(kind=4)              :: i, i3
      type(galaxy)                 :: newgal            ! product of the merger


      ! define i1 as min(i1,i2) and i2 as max(i1,i2) 
      if (i1 > i2) then
         i3 = i1
         i1 = i2
         i2 = i3
      end if
  
      ! define properties of 'newgal', the product of the merger
      call modif_gal_simple(h,list_gals(i1),list_gals(i2),newgal,ss,st)

      ! update list_gals so that all galaxies in list_gals become progs of new list_gals
      ! (this is because the gals in list_gals have already been evolved up to the time of the merger
      ! and so after the merger, we consider their descendents as the new list of galaxies to evolve untill
      ! the end of the timestep or the next merger).
      call copy_gal_props(list_gals(i1),newgal)

      ! decrement nbgals 
      nbgals = nbgals - 1

      ! copy descendents into temporary list 
      do i = i2, nbgals
         call copy_gal_props(list_gals(i),list_gals(i+1))
      end do

      return      

    end subroutine perform_galaxy_merger

!--------------------------------------------------------------------------------------------------
    subroutine modif_gal_simple(h,gal1,gal2,newgal,ss,st)

      ! gal1 and gal2 merger within halo h, and the result of the merger is newgal.
      ! ss is true if sat-sat merger, false if central merger.
      
      ! define here the properties of the merged galaxy
      ! NB: it may happen that a galaxy merges while it has no mass -> in that case, just transfer its accretion budget...

      implicit none

      type(halo)                   :: h          ! halo in which the merger takes place
      type(galaxy)                 :: gal1, gal2 ! the two galaxies that merge (evolved up to time of merging)
      type(galaxy)                 :: newgal     ! the result of the merger
      logical(kind=4)              :: ss         ! if true, the merger is between two satellites
      integer(kind=4)              :: st         ! timestep
      type(galaxy)                 :: g
      real(kind=8)                 :: frac_tran, Agas(3,3), Astar(3,3)
      real(kind=8)                 :: bar_mass1, bar_mass2, m1, m2, mnewgal, mg
      real(kind=8)                 :: r1_half, r2_half, rgal, rhalfmgal
      real(kind=8)                 :: core_mass, core_mass1, core_mass2, sum_coremass, tot_mass, mass_rat
      real(kind=8)                 :: disturbance
      real(kind=8)                 :: e_orbit, r_buxxx 
      integer(kind=4)              :: n1,n2

#ifdef SHELL
      type(shell)                  :: shell
#endif

      call init_galaxy(newgal)
      call init_galaxy(g)

      ! copy most massive gal into newgal and least massive into g
      m1 = total_mass(gal1)
      m2 = total_mass(gal2)
      n1 = gal1%my_progs%nb_prog
      n2 = gal2%my_progs%nb_prog
      if (m1 > m2) then
         ! update major/minor mergers counts : 
         if (m1 > chi * m2) then    ! minor merger
            tsno(st)%nb_min_mergs = tsno(st)%nb_min_mergs + 1
         else                       ! major merger
            tsno(st)%nb_maj_mergs = tsno(st)%nb_maj_mergs + 1
         end if
         call copy_gal_props(newgal,gal1)         
         if (n1 > 0) deallocate(newgal%my_progs%hno,newgal%my_progs%gno)
         call copy_gal_props(g,gal2)
      else
         ! update major/minor mergers counts : 
         if (m2 > chi * m1) then    ! minor merger
            tsno(st)%nb_min_mergs = tsno(st)%nb_min_mergs + 1
         else                       ! major merger
            tsno(st)%nb_maj_mergs = tsno(st)%nb_maj_mergs + 1
         end if
         call copy_gal_props(newgal,gal2)
         if (n2 > 0) deallocate(newgal%my_progs%hno,newgal%my_progs%gno)
         call copy_gal_props(g,gal1)
      end if
      newgal%my_progs%nb_prog         = n1 + n2

      ! jeje : this is not really a problem ... is it ? 
      if (n1+n2 == 0) then 
         write(errunit,*) '> Pb : no progs at all for merging galaxies ... '
         stop
      end if

      if (n1 + n2 > 0) allocate(newgal%my_progs%hno(n1+n2),newgal%my_progs%gno(n1+n2))
      ! have to be careful: if a galaxy is created during a timestep and merges with another 
      ! during this very same timestep, then it has no progenitor at a previous timestep i.e. n1 or n2 = 0
      if (n1 > 0) then 
         newgal%my_progs%hno(1:n1)       = gal1%my_progs%hno
         newgal%my_progs%gno(1:n1)       = gal1%my_progs%gno
         if (n2 > 0) then 
            newgal%my_progs%hno(n1+1:n1+n2) = gal2%my_progs%hno
            newgal%my_progs%gno(n1+1:n1+n2) = gal2%my_progs%gno
         endif
      else
         if (n2 > 0) then 
            newgal%my_progs%hno(1:n2) = gal2%my_progs%hno
            newgal%my_progs%gno(1:n2) = gal2%my_progs%gno
         endif
      endif
      
      ! add up accretion events now (why later ?)
      call merge_accretion_events(newgal%acc,g%acc)       ! add g%acc to newgal%acc
      call merge_accretion_events(newgal%wind,g%wind)     ! add g%wind to newgal%wind
      
#ifdef SHELL
      !call shell_when_merger(newgal%shell,g%shell)
#endif

      ! add bookkeeped stuff
      call merge_gal_bookkeep(newgal%gbk,g%gbk) 

      ! it may happen that g has no mass (accretion hasn't arrived yet...)
      if (total_mass(g) <= 0.0d0) then 
         ! do some bookkeeping
         newgal%nb_merg   = newgal%nb_merg + g%nb_merg + 1 
         newgal%tbirth    = min(newgal%tbirth,g%tbirth)
         ! and leave
         write(errunit,*) 'weird case ... handled smoothly, though... '
         return
      end if

      disturbance                     = newgal%disturb

      ! define what goes where during merger between newgal and g 
      ! (frac_tran is the mass fraction of material (star + gas) transfered from the disc to the bulge-burst combo)
      mnewgal = total_mass(newgal)
      mg = total_mass(g) 
      call mass_recipe_burst(mnewgal,mg,frac_tran,Agas,Astar)

      ! compute props of newgal
      bar_mass1    = (1.d0-frac_tran)*newgal%disc%mgal + newgal%bulge%mgal + newgal%burst%mgal
      if (bar_mass1 > 0.0d0) then
         r_buxxx    = max(newgal%bulge%rgal,newgal%burst%rgal)   ! bulge may have zero size if it has zero mass.
         r1_half    = ((1.d0-frac_tran)*newgal%disc%mgal*newgal%disc%rgal*disc_params%s_to_m +  &
              &  (newgal%bulge%mgal+newgal%burst%mgal)*r_buxxx*bulge_params%s_to_m) / bar_mass1
      else
         ! useless because it will not be used afterwards but prevents e_orbit from being = to infinity
         r1_half    =  newgal%disc%mgal*disc_params%s_to_m
      endif
      call compute_core_mass(newgal,r1_half,core_mass1,'nope')
      ! factor 2 is bc mass is computed in half mass radius so have to correct for that
      m1            = bar_mass1 + newgal%core%fdm*min(2.d0*core_mass1,newgal%core%mass)

      ! compute props of g
      bar_mass2    = (1.d0-frac_tran) * g%disc%mgal + g%bulge%mgal + g%burst%mgal
      if (bar_mass2 > 0.0d0) then
         r_buxxx    = max(g%bulge%rgal,g%burst%rgal)   ! bulge may have zero size if it has zero mass.
         r2_half    = ((1.d0-frac_tran)*g%disc%mgal*g%disc%rgal*disc_params%s_to_m +  &
              &  (g%bulge%mgal+g%burst%mgal)*r_buxxx*bulge_params%s_to_m) / bar_mass2
      else
         ! useless because it will not be used afterwards but prevents e_orbit from being = to infinity
         r2_half    =  g%disc%mgal*disc_params%s_to_m
      endif
      call compute_core_mass(g,r2_half,core_mass2,'nope')
      m2            = bar_mass2 + g%core%fdm*min(2.d0*core_mass2,g%core%mass)


      ! sum of core masses of gal1 and gal2 
      sum_coremass  = newgal%core%fdm*min(2.d0*core_mass1,newgal%core%mass) + &
           & g%core%fdm*min(2.d0*core_mass2,g%core%mass)


      ! now calculate the energy interaction which is going to be deposited in final bulge  
      if (m1 > 0.0d0) then 
         ! disturbance of morphology due to recent merger. assume the amount of disturbed mass follows the mass.  Hence:
         disturbance = (disturbance * m1 + g%disturb * m2) / (m1 + m2)
         ! E_orbit = -0.5 G M_1 M_2 / (r1+r2) per galaxy --> no G because it simplifies out later
         e_orbit     = - m2*m1 / (r1_half+r2_half)
      else
         write(errunit,*) '> Pb in modif_gal: m1 = 0',m1
         stop
      end if

      ! do the actual transfer of stars and gas between components
      call increment_mcold(newgal,g)
      call increment_mcoldz(newgal,g)
      call increment_minstar(newgal,g)     ! also adds SFRs ... 
      call add_gas_props(newgal%disc%mcold,newgal%bulge%mcold,newgal%burst%mcold,Agas)
      call add_gas_props(newgal%disc%mcoldz,newgal%bulge%mcoldz,newgal%burst%mcoldz,Agas)
      call add_star_props(newgal%disc,newgal%bulge,newgal%burst,Astar)

      ! now, the total number of mergers undergone by newgal is the sum of all mergers plus this one
      newgal%nb_merg   = newgal%nb_merg + g%nb_merg + 1 
      ! update time of galaxy birth as the oldest of the two progenitors
      newgal%tbirth    = min(newgal%tbirth,g%tbirth)

      ! how much additional disturbance?
      mass_rat         = min(m2/m1,m1/m2)
      disturbance      = disturbance  + mass_rat
      newgal%disturb   = disturbance       
      
      ! update the mass fields of each components
      call get_mgal(newgal,'in modif_gal')

      ! now define the half mass radius of the final bulge by assuming that both components are virialized
      ! before and after the galaxy merger and galaxy conservation
      if (newgal%bulge%mgal+newgal%burst%mgal > 0.d0) then

         ! pot ener    --> 0.4 factor comes from Hernquist profile for bulges (ApJ 1990)
         newgal%bulge%rgal    = (m1+m2)*(m1+m2) / (-e_orbit/0.4d0+m1*m1/r1_half+m2*m2/r2_half) / bulge_params%s_to_m
         ! in fact do not authorize the bulge radius to be greater than the radius of the host halo
         ! (corresponds to enforcing a minimal orbital energy threshold in absolute value) 
         newgal%bulge%rgal    = min(newgal%bulge%rgal,h%datas%rvir/bulge_params%s_to_m)
         ! limit the minimum bulge size to min_size = 100 pc
         newgal%bulge%rgal    = max(newgal%bulge%rgal,min_size)
         ! assume that the disc which is the more extended survives but in the same time do not allow 
         ! the new disc to be smaller than the bulge 
         newgal%disc%rgal     = max(newgal%disc%rgal,g%disc%rgal,newgal%bulge%rgal*bulge_params%s_to_m/disc_params%s_to_m)
         if (.not. ss) then
            ! compute new core properties before recomputing disc properties if the galaxy is central
            rhalfmgal       = (newgal%disc%mgal*newgal%disc%rgal*disc_params%s_to_m                                           &
                  &         + (newgal%bulge%mgal*newgal%bulge%rgal+newgal%burst%mgal*newgal%burst%rgal)*bulge_params%s_to_m)  &
                  &         / total_mass(newgal)
            call compute_core_mass(newgal,rhalfmgal,core_mass,'calc')
            newgal%core%mass  = min(newgal%core%mass,sum_coremass)
         endif

         call disc_dynamics(newgal,0.0d0,h)
         call compute_core_mass(newgal,newgal%bulge%rgal*bulge_params%s_to_m,core_mass,'nope')
         ! 3D vel disp @ half mass radius --> 0.177 factor comes also from Hernquist profile for bulges (ApJ 1990)
         tot_mass             = newgal%core%fdm*2.d0*core_mass + newgal%bulge%mgal + newgal%burst%mgal  
         newgal%bulge%speed   = sqrt(0.177d0*gravconst*(tot_mass)/newgal%bulge%rgal)
         ! size of burst prop to mass ratio (isothermal sphere approx)
         newgal%burst%rgal    = max(newgal%bulge%rgal*newgal%burst%mgal/(newgal%burst%mgal+newgal%bulge%mgal),min_size) 
         ! velocity of burst is the same as that of bulge --> isothermal sphere approx
         newgal%burst%speed   = newgal%bulge%speed

      else ! happens if there is a minor merger of a pure disc with a pure disc

         ! assume that the disc which is the more extended survives but in the same time do not allow 
         ! the new disc to be smaller than the bulge (for most of the cold gas in the bulge should be in the burst)
         newgal%disc%rgal    = max(newgal%disc%rgal,g%disc%rgal,newgal%bulge%rgal*bulge_params%s_to_m/disc_params%s_to_m)
         if (.not. ss) then
            ! compute new core properties before recomputing disc properties if the galaxy is central
            rhalfmgal       = (newgal%disc%mgal*newgal%disc%rgal*disc_params%s_to_m                                           &
                  &         + (newgal%bulge%mgal*newgal%bulge%rgal+newgal%burst%mgal*newgal%burst%rgal)*bulge_params%s_to_m)  &
                  &         / total_mass(newgal)
            call compute_core_mass(newgal,rhalfmgal,core_mass,'calc')
            newgal%core%mass = min(newgal%core%mass,sum_coremass)
         endif
         call disc_dynamics(newgal,0.0d0,h)

      endif

      ! assign a new orbital radius
      if (ss) then 
         ! if this is a satellite-satellite merger, get the new radius from orbital energy 
         ! but dont let galaxy fall to the center
         newgal%r = min(newgal%r*g%r*(m1+m2)/(m1*g%r+m2*newgal%r),h%rfof)
         rgal     = max(disc_params%s_to_m * newgal%disc%rgal, &
              & bulge_params%s_to_m * newgal%bulge%rgal, bulge_params%s_to_m * newgal%burst%rgal)
         newgal%r = max(newgal%r,rgal)
      else
         ! conversely, for a dynamical friction merger, energy comes from the energy of the two objects orbiting 
         ! around each other at r minimum, and the new radius = 0.          
         newgal%r         = 0.0d0
      end if

      ! compute burst/bulge dynamical timescales (disk has been taken care of in recompute_disc subroutine)
      newgal%bulge%tdyn   = bulge_tdyn(newgal%bulge)       
      newgal%burst%tdyn   = bulge_tdyn(newgal%burst)

      ! determine QS0 mass
      newgal%QSO%mass     = newgal%QSO%mass + g%QSO%mass

      call check_newgal(newgal%disc,'disc')
      call check_newgal(newgal%bulge,'bulg')
      call check_newgal(newgal%burst,'burs')

      return

    end subroutine modif_gal_simple

!--------------------------------------------------------------------------------------------------
  end subroutine evolve_halo

!******************************************************************************************************************
  subroutine star_feed(gal,t1,t2,h,st)

    ! sub to integrate SFR, feedback, metallicity evolution of a galaxy from t1 to t2.
    ! returns .true. if gal is erased because too small ... (this does happen at high res...)
    
    implicit none

    type(halo)                :: h
    type(galaxy)              :: gal
    integer(kind=4)           :: st                ! timestep
    real(kind=8)              :: t1,t2,delta_t
    
    ! time span for evolution 
    delta_t = t2 - t1

    ! compute star formation and feedback during delta_t by substepping each component to ensure proper numerical convergence 
    call sfh_integrator(h,gal,t1,t2,st)

    ! compute amount of stars transfered from burst to bulge
#ifdef DUAL_IMF
    if ((gal%burst%minstar > 0.0d0 .or. gal%burst%minstar2 > 0.0d0) .and. gal%bulge%rgal > gal%burst%rgal) then 
       call burst_transport(gal,delta_t)
    end if
#else
    if (gal%burst%minstar > 0.0d0 .and. gal%bulge%rgal > gal%burst%rgal) call burst_transport(gal,delta_t)
#endif

    ! update total galaxy mass
    call get_mgal(gal,'in star_feed')

    ! remove components which have too little a mass
    if (gal%burst%mgal < 1.d-10) then 
       call erase_comp(gal%burst)
       if (gal%bulge%mgal < 1.d-10) call erase_comp(gal%bulge)  ! erase bulge only if no burst is going to feed it soon 
    end if
    if (gal%disc%mgal  < 1.d-10) call erase_comp(gal%disc)

    ! update total baryonic mass of halo h  
    call get_mhalo(h)

    return

  end subroutine star_feed
  
!*****************************************************************************************************************
  subroutine sfh_integrator(h,gal,t1,t2,st)

    implicit none

    type(halo)                 :: h
    type(galaxy)               :: gal
    type(halo_return)          :: eject
    real(kind=8),intent(in)    :: t1,t2
    integer(kind=4),intent(in) :: st                ! timestep
    real(kind=8)               :: delta_t,subt1,subt2,tlastacc,ej_disc,ej_bulg,ej_burst,mfs_disc,mfs_bulg,mfs_burst,ej_tot,mfs_tot
    integer(kind=4)            :: nsubsteps,i,mytime,itime
#ifdef CLUMPY_SF
    real(kind=8)               :: macc_true_burst, maccz_true_burst
#endif
    real(kind=8)               :: macc_wind, maccz_wind, macc_true, maccz_true,macc_true_disc, maccz_true_disc
    real(kind=8)               :: global_sf,metalsformed
#ifdef DUAL_IMF
    real(kind=8)               :: global_sf2
#endif
    real(kind=8),parameter     :: faccftn2bst = 0.1 !0.5 ! 0.
#ifdef SHELL
    type(shell)                :: shell
#endif

    real(kind=8)               :: zmigr

    ! total mass of gas accreted onto gal from t1 to t2 (include cold streams and ejecta falling back)
    call accretion_mass(gal%acc,t1,t2,macc_true,maccz_true)
    call accretion_mass(gal%wind,t1,t2,macc_wind,maccz_wind)
    if (total_mass(gal) + macc_true + macc_wind <= 0.d0) return ! nothing much to do here ...

    ! define integration timestep to be 1Myr. 
    ! NB: This is ok for SP99 feedback scheme. S01 scheme requires an adaptive step for numerical convergence... 
    !  -> see old version of the code. 
    delta_t   = delta_timetab(1)  ! 1 Myr (exactly!)
    nsubsteps = nint((t2-t1)/delta_t) + 1
    delta_t   = (t2-t1) / dble(nsubsteps) ! a 'chouille' less than 1 Myr ... 


    ! reallocate sfhtabs to hold forthcoming evolution
    call resize_sfhtab(gal%disc,t1,t2,gal%tbirth)
    call resize_sfhtab(gal%bulge,t1,t2,gal%tbirth)
    call resize_sfhtab(gal%burst,t1,t2,gal%tbirth)
#ifdef RECORD_SFR
    call resize_sfrtab(gal%disc,t1,t2,gal%tbirth)
    call resize_sfrtab(gal%bulge,t1,t2,gal%tbirth)
    call resize_sfrtab(gal%burst,t1,t2,gal%tbirth)
#endif

    global_sf   = 0.0d0
#ifdef DUAL_IMF
    global_sf2  = 0.0d0
#endif
    eject%mhot   = 0.0d0
    eject%mhotz  = 0.0d0
    eject%mout   = 0.0d0
    eject%moutz  = 0.0d0
    metalsformed = 0.0d0  ! metals formed by stars and released to ism during timestep
    tlastacc     = t1
    
    gal%nsubsteps = nsubsteps 

    do i = 1, nsubsteps 

       ! jeje : general note about the two lines below ... 
       ! It turns out that numerical precision is much better when subt1 and subt2 are defined as below
       ! (instead of subt2 = subt1 + delta_t which i originally used), "much better" meaning a few digits better (2-3)!!! 
       ! also, in the test below (macc_true > ...), the factor nsubsteps is somewhat minimal, and could be turned higher for better
       ! results ... although it then means that we never allow accretion of more than a thousandth of the present mass, which 
       ! becomes physically borderline. 
       subt1 = t1 + (i-1) * delta_t   
       subt2 = t1 + i * delta_t

       ! accrete onto the disc component the mass of gas that reaches the galaxy from tlastacc to subt2
       call accretion_mass(gal%acc,subt1,subt2,macc_true,maccz_true)
       call accretion_mass(gal%wind,subt1,subt2,macc_wind,maccz_wind)
#ifdef CLUMPY_SF

#ifdef Z_MIGRATION
       if (global_redshift .ge. 6.) then
          zmigr = 1.0
       else
! 3.57 is the normalization at z=6, so that zmigr = 1.0
          zmigr = ((1.+global_redshift)/3.)**1.5  / 3.57
       endif
       if ((zmigr * frac_acc_to_burst) .gt. 1.0) then
          stop
       end if
      ! print*,'Redshift ',global_redshift,' Zmigr = ',zmigr,' Zmigr x fracbst = ',zmigr*frac_acc_to_burst
#else ! no Z_MIGRATION
       zmigr = 1.0
#endif ! end Z_MIGRATION
       
#ifdef ALL2BURST
       macc_true_disc   = macc_true  * (1.d0 - frac_acc_to_burst * zmigr) + macc_wind * (1.d0 - faccftn2bst)
       maccz_true_disc  = maccz_true * (1.d0 - frac_acc_to_burst * zmigr) + maccz_wind * (1.d0 - faccftn2bst)
       macc_true_burst  = macc_true  * frac_acc_to_burst * zmigr + macc_wind * faccftn2bst
       maccz_true_burst = maccz_true * frac_acc_to_burst * zmigr + maccz_wind * faccftn2bst
#else
       macc_true_disc   = macc_wind  + macc_true  * (1.d0 - frac_acc_to_burst * zmigr)
       maccz_true_disc  = maccz_wind + maccz_true * (1.d0 - frac_acc_to_burst * zmigr)
       macc_true_burst  = macc_true  * frac_acc_to_burst * zmigr
       maccz_true_burst = maccz_true * frac_acc_to_burst * zmigr
#endif
#else ! no CLUMPY_SF
       macc_true_disc   = macc_wind + macc_true
       maccz_true_disc  = maccz_wind + maccz_true
#endif

       ! add gas only if numerically significant (with respect to mcold). 

!++++++++++++++++++++++++++++ Comment for analytical gas evolution +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef CLUMPY_SF
       if (macc_true_disc > rel_prec * gal%disc%mcold * nsubsteps &
            & .or. macc_true_burst > rel_prec * gal%burst%mcold * nsubsteps ) then  ! NB: factor nsubsteps could be higher to ensure even better accuracy ... 
#else
       if (macc_true_disc > rel_prec * gal%disc%mcold * nsubsteps) then 
#endif
          call disc_dynamics(gal,macc_true_disc,h)  ! update dynamics before adding the mass (need mass-weighted radii averages here)
          gal%disc%mcold   = gal%disc%mcold  + macc_true_disc
          gal%disc%mcoldz  = gal%disc%mcoldz + maccz_true_disc
          gal%disc%mgal    = gal%disc%mgal   + macc_true_disc
#ifdef CLUMPY_SF
          gal%burst%mcold  = gal%burst%mcold  + macc_true_burst
          gal%burst%mcoldz = gal%burst%mcoldz + maccz_true_burst
          gal%burst%mgal   = gal%burst%mgal   + macc_true_burst
          call buxxx_dynamics(gal,h)
#endif
          tlastacc         = subt2
       end if

!++++++++++++++++++++++++++++ STOP Comment for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! evolve the stars during delta_t and add metal ejecta to cold phase
       call evolve_stars(gal%disc,delta_t,metalsformed)
       call evolve_stars(gal%bulge,delta_t,metalsformed)
       call evolve_stars(gal%burst,delta_t,metalsformed)
       
#ifdef RECORD_SFR
       ! evolve SFR table (simple book-keeping)
       call evolve_sfr_tab(gal%disc,delta_t)
       call evolve_sfr_tab(gal%bulge,delta_t)
       call evolve_sfr_tab(gal%burst,delta_t)
#endif

!++++++++++++++++++++++++++++ Comment for analytical gas evolution +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ej_disc   = 0.0d0
       ej_bulg   = 0.0d0
       ej_burst  = 0.0d0
       mfs_disc  = 0.0d0
       mfs_bulg  = 0.0d0
       mfs_burst = 0.0d0

       mytime = 0.5 * (subt1 + subt2)
       itime = int(mytime)

       ! compute star formation and feedback in each component during delta_t
       if (gal%disc%mcold > 0.0d0) then 
#ifdef DUAL_IMF
          call comp_sf_n_fb(gal,'disc',subt1,delta_t,global_sf,global_sf2,eject,st,h,ej_disc,mfs_disc,i)
#else
          call comp_sf_n_fb(gal,'disc',subt1,delta_t,global_sf,eject,st,h,ej_disc,mfs_disc,i)
#endif
       end if
       if (gal%bulge%mcold > 0.0d0) then 
#ifdef DUAL_IMF
          call comp_sf_n_fb(gal,'bulg',subt1,delta_t,global_sf,global_sf2,eject,st,h,ej_bulg,mfs_bulg,i)
#else
          call comp_sf_n_fb(gal,'bulg',subt1,delta_t,global_sf,eject,st,h,ej_bulg,mfs_bulg,i)
#endif         
       end if
       if (gal%burst%mcold > 0.0d0) then 
#ifdef DUAL_IMF
          call comp_sf_n_fb(gal,'burs',subt1,delta_t,global_sf,global_sf2,eject,st,h,ej_burst,mfs_burst,i)
#else
          call comp_sf_n_fb(gal,'burs',subt1,delta_t,global_sf,eject,st,h,ej_burst,mfs_burst,i)
#endif
       end if
 

!tibo: bookkeep substeps' ejecta
!!$       ej_tot  = 0.0d0
!!$       mfs_tot = 0.0d0
!!$
!!$       ej_tot  =  ej_burst !ej_disc + ej_bulg
!!$       mfs_tot = mfs_disc + mfs_bulg + mfs_burst
       !print*, ej_disc, ej_bulg,  ej_burst
       
!!$       if (st == 19) then 
!!$          print*,ej_tot,mfs_tot      
!!$       endif
!!$       if (st == 19 .and. i == nsubsteps) then
!!$          print*,gal%disc%mgal
!!$          stop
!!$       endif

!++++++++++++++++++++++++++++ STOP Comment for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
#ifndef NO_DISC_INSTAB
       if (gal%disc%mgal > 0.0d0) then 
          call check_for_disc_instability(h,gal) ! and actually do the transfers if there is one ... 
       end if
#endif    
   
!tibo: output ej_cold at each substep (sum of the 3 components)
       
!!$       if (ISNAN(ej_tot) .or. ISNAN(mfs_tot)) then
!!$          print*,'ej_tot or mfs_tot is NaN'
!!$          stop
!!$       endif
!!$
!!$       if (subt1 >= tsno(st)%age_univm100) then 
!!$          gal%wind100 = gal%wind100 + ej_tot 
!!$          if (subt1 >= tsno(st)%age_univm50) then 
!!$             gal%wind50 = gal%wind50 + ej_tot 
!!$             if (subt1 >= tsno(st)%age_univm20) then 
!!$                gal%wind20 = gal%wind20 + ej_tot 
!!$                if (subt1 >= tsno(st)%age_univm10) then 
!!$                   gal%wind10 = gal%wind10 + ej_tot 
!!$                   if (subt1 >= tsno(st)%age_univm2) then
!!$                      gal%wind2 = gal%wind2 + ej_tot
!!$                      gal%wind1  = ej_tot
!!$                   endif
!!$                endif
!!$             endif
!!$          end if
!!$       end if
!!$         
!!$       gal%galej(i) = ej_tot
!!$       gal%burst_ej(i) = ej_burst

    end do

    ! finish up remaining accretion if any left 
    if (tlastacc < t2) then 
       call accretion_mass(gal%acc,tlastacc,t2,macc_true,maccz_true)
       call accretion_mass(gal%wind,tlastacc,t2,macc_wind,maccz_wind)
#ifdef CLUMPY_SF
#ifdef ALL2BURST
       gal%disc%mcold  = gal%disc%mcold  + macc_true  * (1.d0 - frac_acc_to_burst * zmigr) + macc_wind * (1.d0 - faccftn2bst)
       gal%disc%mcoldz = gal%disc%mcoldz + maccz_true * (1.d0 - frac_acc_to_burst * zmigr) + maccz_wind * (1.d0 - faccftn2bst)
       gal%disc%mgal   = gal%disc%mgal   + macc_true  * (1.d0 - frac_acc_to_burst * zmigr)  + macc_wind * (1.d0 - faccftn2bst)
       gal%burst%mcold  = gal%burst%mcold  + macc_true  * frac_acc_to_burst * zmigr + macc_wind * faccftn2bst
       gal%burst%mcoldz = gal%burst%mcoldz + maccz_true * frac_acc_to_burst * zmigr + maccz_wind * faccftn2bst
       gal%burst%mgal   = gal%burst%mgal   + macc_true  * frac_acc_to_burst * zmigr + macc_wind * faccftn2bst
#else
       gal%disc%mcold  = gal%disc%mcold  + macc_true  * (1.d0 - frac_acc_to_burst * zmigr) + macc_wind
       gal%disc%mcoldz = gal%disc%mcoldz + maccz_true * (1.d0 - frac_acc_to_burst * zmigr) + maccz_wind
       gal%disc%mgal   = gal%disc%mgal   + macc_true  * (1.d0 - frac_acc_to_burst * zmigr) + macc_wind
       gal%burst%mcold  = gal%burst%mcold  + macc_true  * frac_acc_to_burst * zmigr
       gal%burst%mcoldz = gal%burst%mcoldz + maccz_true * frac_acc_to_burst * zmigr
       gal%burst%mgal   = gal%burst%mgal   + macc_true  * frac_acc_to_burst * zmigr
       call buxxx_dynamics(gal,h)
#endif
#else
       gal%disc%mcold  = gal%disc%mcold  + macc_true  + macc_wind
       gal%disc%mcoldz = gal%disc%mcoldz + maccz_true + maccz_wind
       gal%disc%mgal   = gal%disc%mgal   + macc_true  + macc_wind
#endif
    end if
    tsno(st)%global_sf  = tsno(st)%global_sf  + global_sf
#ifdef DUAL_IMF
    tsno(st)%global_sf2 = tsno(st)%global_sf2 + global_sf2
#endif

    ! update baryon content of halo after full timestep
    h%datas%mhotz     = h%datas%mhotz     + eject%mhotz
    h%datas%mhotgaz   = h%datas%mhotgaz   + eject%mhot
    h%datas%mcoldgaz  = h%datas%mcoldgaz  - eject%mhot  - eject%mout
    h%datas%mcoldz    = h%datas%mcoldz    - eject%mhotz - eject%moutz + metalsformed
    h%datas%mgazout   = h%datas%mgazout   + eject%mout
    h%datas%metalsout = h%datas%metalsout + eject%moutz
    
    ! clean up SFHTABs allocated uselessly. 
    call clean_up_sfhtab(gal%disc)
    call clean_up_sfhtab(gal%bulge)
    call clean_up_sfhtab(gal%burst)

    return

  end subroutine sfh_integrator

!*****************************************************************************************************************
  subroutine clean_up_sfhtab(comp)

    implicit none
    
    type(gal_comp) :: comp

    if (comp%minstar <= 0.0d0) then 
       if (associated(comp%sfh_tab)) deallocate(comp%sfh_tab)
       comp%minstar = 0.0d0
    end if
#ifdef DUAL_IMF 
    if (comp%minstar2 <= 0.0d0) then 
       if (associated(comp%sfh_tab2)) deallocate(comp%sfh_tab2)
       comp%minstar2 = 0.0d0
    end if
#endif
    
    return
    
  end subroutine clean_up_sfhtab

!*****************************************************************************************************************
#ifdef DUAL_IMF
  function chose_imf(comp)
    
    implicit none
    
    type(gal_comp) :: comp
    logical(kind=4) :: chose_imf
    
#if (DUAL_IMF == 1)
    if ((comp%tdyn < tdyn_threshold).and.(global_redshift.gt.z_threshold)) then ! use IMF2
       chose_imf = .true.
    else ! use IMF1
       chose_imf = .false.
    end if
#elif (DUAL_IMF == 2)
    if ((comp%mcold .gt. m_threshold) .and. (global_redshift.gt.z_threshold)) then ! use IMF2
       chose_imf = .true.
    else ! use IMF1
       chose_imf = .false.
    end if
#elif (DUAL_IMF == 3) 
    chose_imf = .true.
#endif
    
    return
    
  end function chose_imf
#endif
!*****************************************************************************************************************
  subroutine calculate_disturbance(gal,dt)

    implicit none
 
    type (galaxy)   :: gal
    real(kind=8)    :: dt,tdyn,mass
  
    ! calculate disturbance.  also ensure that gals keep correct dynamical timescale.  
    mass = total_mass(gal)
    if (mass > 0.0d0) then 
       tdyn = gal%disc%tdyn*gal%disc%mgal + gal%bulge%tdyn*gal%bulge%mgal + gal%burst%tdyn*gal%burst%mgal 
       tdyn = tdyn/mass
       if (dt/tdyn <= 20.0d0) then
          gal%disturb = gal%disturb * exp(-dt/tdyn)      
       else
          gal%disturb = 0.0d0
       end if
    else 
       ! it may now happen that a galaxy has zero mass (but some accretion is on its way...). 
       ! -> can't really be disturbed, then ... 
       gal%disturb = 0.0d0 
    end if

    return

  end subroutine calculate_disturbance

!*****************************************************************************************************************

  subroutine evolve_stars(comp,delta_t,metalsformed)

    ! routine which evolves a population of stars, i.e. which changes the distribution of masses in the grid 
    ! age vs metallicity.

    implicit none

    type(gal_comp)           :: comp      ! component of gal to evolve
    real(kind=8)             :: delta_t   ! during delta_t
    integer(kind=4)          :: i,j,ind_star
    real(kind=8),allocatable :: stars(:,:)
    real(kind=8)             :: mej,ej_inst,mold_ejz,mass_stars,mejz,gcalc,metalsformed
#ifdef DUAL_IMF
    real(kind=8)             :: mej2,mold_ejz2,mejz2
#endif

    mej   = 0.0d0 ; mejz = 0.0d0 ; mold_ejz = 0.0d0
    if (comp%minstar > 0.0d0) then 
       ind_star = size(comp%sfh_tab,dim=1)
       allocate(stars(ind_star,nfile))
       stars = comp%sfh_tab
       do i = 1, nfile           ! loop on metallicity cells
          do j = 2, ind_star     ! loop on time cells
             ! computation of ejecta
             ej_inst             = stars(j-1,i) * delta_t
             gcalc               = ej_inst * gaztab(j-1,i)
             mejz                = mejz         + ej_inst * ztab(j-1,i)
             mej                 = mej          + gcalc
             mold_ejz            = mold_ejz     + gcalc * tabmetspec(i)
             ! stellar mass loss due to ejecta
             stars(j-1,i)        = stars(j-1,i) - gcalc
             comp%sfh_tab(j-1,i) = stars(j-1,i)
             ! stellar mass loss due to evolution (stars evolve from cell j-1 to cell j)
             stars(j-1,i)        = stars(j-1,i) * (1.0d0 - delta_t / delta_timetab(j-1))
          end do
       end do
       do i = 1, nfile             ! loop on metallicity cells
          do j = 2, ind_star       ! loop on time cells
             ! stellar gain due to evolution (stars that were in cell j-1 are now in cell j)
             stars(j,i) = stars(j,i) + comp%sfh_tab(j-1,i) * delta_t / delta_timetab(j-1)
          end do
       end do
       comp%sfh_tab = stars
       mass_stars   = sum(stars)
       deallocate(stars)
       ! update minstars 
       comp%minstar = mass_stars
    end if

#ifdef DUAL_IMF
    mej2   = 0.0d0 ; mejz2 = 0.0d0 ; mold_ejz2 = 0.0d0
    if (comp%minstar2 > 0.0d0) then 
       ind_star = size(comp%sfh_tab2,dim=1)
       allocate(stars(ind_star,nfile))
       stars = comp%sfh_tab2
       do i = 1, nfile           ! loop on metallicity cells
          do j = 2, ind_star     ! loop on time cells
             ! computation of ejecta
             ej_inst              = stars(j-1,i) * delta_t
             gcalc                = ej_inst * gaztab2(j-1,i)
             mejz2                = mejz2         + ej_inst * ztab2(j-1,i)
             mej2                 = mej2          + gcalc
             mold_ejz2            = mold_ejz2     + gcalc * tabmetspec(i)
             ! stellar mass loss due to ejecta
             stars(j-1,i)        = stars(j-1,i) - gcalc
             comp%sfh_tab2(j-1,i) = stars(j-1,i)
             ! stellar mass loss due to evolution (stars evolve from cell j-1 to cell j)
             stars(j-1,i)        = stars(j-1,i) * (1.0d0 - delta_t / delta_timetab(j-1))
          end do
       end do
       do i = 1, nfile                  ! loop on metallicity cells
          do j = 2, ind_star       ! loop on time cells
             ! stellar gain due to evolution (stars that were in cell j-1 are now in cell j)
             stars(j,i) = stars(j,i) + comp%sfh_tab2(j-1,i) * delta_t / delta_timetab(j-1)
          end do
       end do
       comp%sfh_tab2 = stars
       mass_stars    = sum(stars)
       deallocate(stars)
       ! update minstars 
       comp%minstar2 = mass_stars
    end if
#endif

    ! add new metals to cold phase
#ifdef DUAL_IMF 
    comp%mcold   = comp%mcold   + (mej + mej2)
    comp%mcoldz  = comp%mcoldz  + (mold_ejz + mejz + mold_ejz2 + mejz2)
    metalsformed = metalsformed + (mold_ejz + mejz + mold_ejz2 + mejz2)
#else 
    comp%mcold   = comp%mcold   + mej
    comp%mcoldz  = comp%mcoldz  + (mold_ejz + mejz)
    metalsformed = metalsformed + (mold_ejz + mejz)
#endif

    call get_mcomp(comp)

    return

  end subroutine evolve_stars

!*****************************************************************************************************************
!!$  subroutine ram_pres_strip(h,nbins,delta_t)
!!$! B. Lanzoni rampressure stripping code (slightly modified by J.D. :), & re-modified by BL 22/1/04)
!!$
!!$    implicit none
!!$
!!$    type(halo)                  :: h,h_temp
!!$    integer(kind=4)             :: i,j,ii,nbins,nbgal,nbgal_wiped
!!$    real(kind=8)                :: delta_t
!!$    real(kind=8)                :: r,r_gas(nbins),alpha_gas,c_gas,rhoicm,r_disc,den
!!$    real(kind=8)                :: sig0_stgas,sig0_gas,exp_rstr,mass_gal
!!$    real(kind=8)                :: r_str,mstrip,strip_met
!!$    real(kind=8)                :: vel_disp,vel_perp
!!$    logical(kind=4)             :: gal_wipe
!!$    integer(kind=4),allocatable :: wipe_gal_number(:)
!!$    type(galaxy), pointer       :: gal
!!$
!!$    ! No ICM ==> stripping impossible ==> do not even start:
!!$    if (h%datas%mhotgaz == 0.) then
!!$      return 
!!$    endif
!!$
!!$    ! build a radius array the same size as the array containing the tabulated hot gas density profile 
!!$    do i=1,nbins
!!$       r_gas(i) = real(i-1)*min_size
!!$    enddo
!!$
!!$    ! Compute terms of ICM density that are independent of galaxy position:
!!$    ! in this way all gals of the current halo will see the same hot gas 
!!$    ! distribution. If I compute these terms inside the gals loop, in case of
!!$    ! stripping of one gal, the others gals in the halo would see more hot 
!!$    ! gas ==> would suffer more stripping. But since gals are analysed in a 
!!$    ! order that has no physical meaning, this would not be reasonable. 
!!$    ! Halo 1D velocity dispersion:
!!$    vel_disp    = sqrt(0.5*gravconst*h%mfof/h%rfof)  
!!$!!    if (h%datas%rvir /= 0)  vel_disp = sqrt(0.5*gravconst*h%datas%mvir/h%datas%rvir)
!!$
!!$    ! Start stripping gals if needed:
!!$    gal_wipe        = .false.
!!$    nbgal_wiped     = 0
!!$    allocate(wipe_gal_number(h%datas%nbgal))
!!$    wipe_gal_number = 0
!!$    do_gals: do i=1,h%datas%nbgal
!!$
!!$       gal => h%datas%liste_galaxies(i)
!!$
!!$       ! If the gal has no cold gas, skip it:
!!$       if (gal%disc%mcold <= 0.0) cycle do_gals
!!$       ! If the gal is in the halo centre, skip it:
!!$       call check_orbit_radius(gal)
!!$       r  = gal%r    ! gal position
!!$       if (r == 0.) cycle do_gals
!!$
!!$       ! Otherwise, go on:
!!$       r_str     = 0.0
!!$       mstrip    = 0.0
!!$       strip_met = 0.0
!!$
!!$       ! Compute rho_ICM at gal position r:
!!$       if (r < r_gas(nbins)) then 
!!$           ! galaxy is located within r_vir --> no problem hot gas profile is defined 
!!$           call locate(r_gas,nbins,r,ii)
!!$           rhoicm     = rho_gas(ii) 
!!$       else
!!$           ! galaxy is between r_fof and r_vir --> have to extrapolate hot gas profile
!!$           ! we choose to do so with a power law here: rho_gas = c_gas r_gas^alpha_gas
!!$            alpha_gas  = (alog10(rho_gas(nbins))-alog10(rho_gas(nbins-1))) / &
!!$                & (alog10(r_gas(nbins))-alog10(r_gas(nbins-1)))
!!$            c_gas      = 10**(alog10(rho_gas(nbins-1)) - alpha_gas*alog10(r_gas(nbins-1)))
!!$            rhoicm     = c_gas*r**alpha_gas
!!$       endif
!!$       if (rhoicm <= 0.) then 
!!$           write(errunit,*)'rhoicm <=0 !!',ii,rhoicm,h%datas%mhotgaz 
!!$           write(errunit,*)'r,rgas',r,r_gas(ii),rho_gas(ii)
!!$           write(errunit,*)'rgasmax',r_gas(nbins),rho_gas(nbins)
!!$           write(errunit,*)'rvir',h%datas%rvir,c_gas,alpha_gas
!!$           stop
!!$       endif   
!!$
!!$       ! Compute gal velocity perpendicular to the disc: - 1 dimensional velocity (absolute).
!!$       vel_perp   = abs(gasdev(iseed) * vel_disp)
!!$       if (vel_perp <= 0.0) write(errunit,*) 'vel_perp <= 0 in ram_pres_strip',vel_perp
!!$
!!$       ! Compute central surface brightness of the exponential disc:
!!$       r_disc     = gal%disc%rgal
!!$       den        = 2.0*pi*r_disc*r_disc
!!$       sig0_stgas = gal%disc%mgal/den  ! --> cold_gas+stars
!!$       sig0_gas   = gal%disc%mcold/den ! --> cold gas
!!$
!!$       ! Compute stripping radius (normalized to r_disc):
!!$       den        = 2.0*pi*gravconst*sig0_stgas*sig0_gas
!!$       exp_rstr   = sqrt((rhoicm/den)*vel_perp*vel_perp)
!!$       r_str      = -log(exp_rstr)
!!$
!!$       ! Compute effects of ram pressure on the galaxy:
!!$       ! --> r_strip_if, CASE 1: so strong ram pressure that ALL gas is stripped out
!!$       !     NB: r_str <= 0.15 is already sufficient for ALL the gas to
!!$       !         be stripped, since M(<r_str) = Exp[-0.15]*(1+0.15) ~= 0.99
!!$       if (r_str <= 0.15) then 
!!$           mstrip          = gal%disc%mcold
!!$           strip_met       = gal%disc%mcoldz
!!$           ! Take off stripped gas and metals from the disc:
!!$           gal%disc%mcold  = 0.0
!!$           gal%disc%mcoldz = 0.0
!!$           gal%disc%rstrip = 0.0
!!$       !--> r_strip_if, CASE 2: partial stripping
!!$       else if (r_str < gal%disc%rstrip) then 
!!$           mstrip          = gal%disc%mcold * exp_rstr * (1.0+r_str)
!!$           mstrip          = min(gal%disc%mcold,mstrip)
!!$           strip_met       = (mstrip/gal%disc%mcold) * gal%disc%mcoldz
!!$           strip_met       = min(strip_met,gal%disc%mcoldz)
!!$           gal%disc%rstrip = r_str ! --> normalized to r_disc 
!!$           ! Take off stripped gas and metals from the disc:
!!$           gal%disc%mcold  = gal%disc%mcold  - mstrip
!!$           gal%disc%mcoldz = gal%disc%mcoldz - strip_met
!!$       endif
!!$
!!$       ! If almost all the mass in the disc has been stripped out (pure gas disc, without stars)
!!$       !   ==> the whole disc disappears:
!!$       if (gal%disc%mgal - mstrip <= 1.0e-10) then
!!$           write(errunit,*) '> Warning: disc vanished bc ram press stripping', &
!!$                 & h%my_timestep, h%my_number, gal%my_number
!!$           mstrip = max(gal%disc%mgal, mstrip)
!!$           call erase_comp(gal%disc) ! ==> gal%disc%mgal = 0
!!$           ! If the galaxy was a pure disc, put to zero all its components and
!!$           ! wipe it out cleanly from the halo galaxy list at the end of the routine...
!!$           mass_gal = gal%bulge%mgal + gal%burst%mgal
!!$           if (mass_gal <= 1.0e-10) then
!!$               write(errunit,*) '> Warning: galaxy vanished bc ram press stripping', &
!!$                     & h%my_timestep, h%my_number, gal%my_number
!!$               mstrip = mstrip + mass_gal
!!$               strip_met = strip_met + gal%bulge%mcoldz + gal%burst%mcoldz
!!$               call erase_comp(gal%disc)
!!$               call erase_comp(gal%bulge)
!!$               call erase_comp(gal%burst)
!!$               gal_wipe           = .true.
!!$               nbgal_wiped        = nbgal_wiped + 1
!!$               wipe_gal_number(i) = gal%my_number 
!!$           endif 
!!$       else
!!$           gal%disc%mgal = gal%disc%mgal - mstrip
!!$       endif
!!$
!!$       ! Compute effects of ram pressure on the halo:
!!$       ! Add stripped gas and metals to hot gas and metals in the halo
!!$       mstrip           = min(h%datas%mcoldgaz,mstrip)
!!$       h%datas%mcoldgaz = h%datas%mcoldgaz - mstrip 
!!$       h%datas%mhotgaz  = h%datas%mhotgaz  + mstrip
!!$       strip_met        = min(h%datas%mcoldz,strip_met)
!!$       h%datas%mcoldz   = h%datas%mcoldz   - strip_met
!!$       h%datas%mhotz    = h%datas%mhotz    + strip_met
!!$
!!$    end do do_gals
!!$
!!$    ! remove galaxies from halo galaxy list if they have been wiped out by RAM pres stripping 
!!$    if (gal_wipe) then 
!!$       call copy_halo_props(h_temp,h)
!!$       nbgal = h_temp%datas%nbgal - nbgal_wiped
!!$       call alloc_hlist(h,nbgal)
!!$       ! if all gals in the halo have been completely wiped out, exactly set its gas properties:
!!$       if (nbgal == 0) then
!!$           h%datas%mcoldgaz = 0.0
!!$           h%datas%mgaz = h%datas%mhotgaz
!!$           h%datas%mcoldz = 0.0
!!$           h%datas%cooling_rate=0.
!!$           h%datas%zcooling_rate=0.
!!$       else 
!!$         j = 1
!!$         do i = 1,h_temp%datas%nbgal
!!$            if (h_temp%datas%liste_galaxies(i)%my_number == wipe_gal_number(i)) then 
!!$               if (h_temp%datas%liste_galaxies(i)%cooling_frac == 1.0) then
!!$                  h%datas%mcoldgaz = h%datas%mcoldgaz - h%datas%cooling_rate*delta_t
!!$                  h%datas%cooling_rate=0.
!!$                  h%datas%mcoldz = h%datas%mcoldz - h%datas%zcooling_rate*delta_t
!!$                  h%datas%zcooling_rate=0.
!!$               endif
!!$            else
!!$               call copy_gal_props(h%datas%liste_galaxies(j),h_temp%datas%liste_galaxies(i))
!!$               j = j+1
!!$            endif
!!$         end do
!!$         ! now renumber the galaxies between 1 and h%datas%nbgal to avoid confusion
!!$         do i = 1,h%datas%nbgal
!!$            h%datas%liste_galaxies(i)%my_number = i
!!$         end do   
!!$      endif
!!$   endif
!!$
!!$   deallocate(wipe_gal_number)
!!$
!!$    return
!!$      
!!$  end subroutine ram_pres_strip
!!$
!*****************************************************************************************************************
  
#ifdef DUAL_IMF
  subroutine comp_sf_n_fb(gal,comp_type,subt1,delta_t,global_sf,global_sf2,ejectsave,st,h,ej_comp,mfs_comp,i)
#else
  subroutine comp_sf_n_fb(gal,comp_type,subt1,delta_t,global_sf,ejectsave,st,h,ej_comp,mfs_comp,i)
#endif

    implicit none

    type(galaxy)      :: gal
    type(halo)        :: h
    type(ejecta)      :: eject
    type(halo_return) :: ejectsave
    character(4)      :: comp_type
    real(kind=8)      :: delta_t,mfs,metfs,subt1,delay_fountain2
    real(kind=8)      :: mscale,global_sf,ej_comp,mfs_comp
    integer(kind=4)   :: st,i,itime
    real(kind=8)      :: accrate,accratez,t1wind,t2wind,htdyn
#ifdef DUAL_IMF 
    real(kind=8)      :: global_sf2
    logical(kind=4)   :: use_second_imf
#endif

!tibo: modify delay_fountain
    
    !if (gal%mgal < 0.1) then 
    !   delay_fountain2 = 0.1
    !else
    !   delay_fountain2 = 1.
    !end if

    if (comp_type == 'disc') then 
#ifdef DUAL_IMF
       ! chose an IMF ... 
       use_second_imf = chose_imf(gal%disc)
#endif          

       ej_comp  = 0.0d0
       mfs_comp = 0.0d0

!++++++++++++++++++++++++++++ Comment for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!       call compute_mass_for_stars(gal%disc,comp_type,delta_t,mfs,metfs)
!#ifdef DUAL_IMF
!       call compute_ejecta(gal%disc,comp_type,eject,mfs,use_second_imf)
!#else
!       call compute_ejecta(gal%disc,comp_type,eject,mfs)
!#endif

!++++++++++++++++++++++++++++ STOP Comment for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++ CALL the EVOLVE_GAS routine for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++



       !call evolve_gas(comp_type,mfs,metfs,eject)
       call evolve_gas(gal%disc,comp_type,delta_t,mfs,metfs,eject)
       
       ej_comp  = eject%mcold
       mfs_comp = mfs 

!------------------------------------- STOP ANALYTICAL GAS EVOLUTION -------------------------------------------------

       ! if mass of formed stars + ejected gas > mass of ism, re-normalize
       if (mfs + eject%mcold + eject%mhot > gal%disc%mcold) then
          write(errunit,*) 're-normalization of feedback in disc... fishy ... '
          mscale      = gal%disc%mcold / (mfs + eject%mcold + eject%mhot) 
          mfs         = mfs * mscale
          metfs       = metfs * mscale
          eject%mcold  = eject%mcold  * mscale
          eject%mcoldz = eject%mcoldz * mscale
          eject%mhot   = eject%mhot   * mscale
          eject%mhotz  = eject%mhotz  * mscale
          if ((mfs + eject%mcold + eject%mhot - gal%disc%mcold)/gal%disc%mcold > rel_prec) then 
             write(errunit,*) 'oh no for gas in sf_step... disc'
             stop
          end if
          if ((metfs + eject%mcoldz + eject%mhotz - gal%disc%mcoldz)/gal%disc%mcoldz > rel_prec) then 
             write(errunit,*) 'oh no for metals in sf_step... disc'
             stop
          end if
       end if

       ! update mcold(z) according to mass of ejected gas and mass of stars formed
       gal%disc%mcold    = max(0.0d0,gal%disc%mcold  - eject%mcold  - eject%mhot  - mfs)
       gal%disc%mcoldz   = max(0.0d0,gal%disc%mcoldz - eject%mcoldz - eject%mhotz - metfs)
       ! The cold phase of ejecta either mixes with hot halo phase if it exists, or 
       ! they will fall back onto the galaxy. The hot phase will escape the halo if it has no
       ! hot phase already in place, or mix with it otherwise.
! jeje 
       if (h%datas%mhotgaz > h%mfof) then !0.0d0) then 
          print*,'this should never happen... '
          ! The whole wind (hot+cold) mixes with the halo's hot phase
          ejectsave%mhot   = ejectsave%mhot  + eject%mcold  + eject%mhot
          ejectsave%mhotz  = ejectsave%mhotz + eject%mcoldz + eject%mhotz
          call inc_gbk_mej(gal%gbk,0.0d0,eject%mcold  + eject%mhot,0.0d0)
       else
          ! the ejected cold gas leaves the galaxy but returns later
          htdyn    = h%datas%rvir / h%datas%cvel * 978.1d0 ! in Gyr
          !tibo:
          t1wind   = subt1 + delay_fountain * htdyn + 2.0d0*accretion_deltat
          !t1wind   = subt1 + 0.2 + 2.0d0*accretion_deltat

          t2wind   = t1wind + duration_fountain * htdyn 
          accrate  = eject%mcold  / (t2wind - t1wind)
          accratez = eject%mcoldz / (t2wind - t1wind)
          call add_accretion_event(gal%wind,accrate,accratez,t1wind,t2wind)
          ! the ejected hot gas leaves the halo forever
          ejectsave%mout  = ejectsave%mout  + eject%mhot
          ejectsave%moutz = ejectsave%moutz + eject%mhotz
          call inc_gbk_mej(gal%gbk,eject%mcold,0.0d0,eject%mhot)
       end if

#ifdef DUAL_IMF
       if (mfs > 0.0) then
          call add_new_stars_to_sfhtab(gal%disc,mfs,metfs,use_second_imf)
#ifdef RECORD_SFR           
          call add_new_stars_to_sfrtab(gal%disc,mfs,metfs,use_second_imf)
#endif          
          ! increment SFR counters 
          if (use_second_imf) then 
             if (subt1 >= tsno(st)%age_univm100) then 
                gal%disc%sfr2100 = gal%disc%sfr2100 + mfs
                if (subt1 >= tsno(st)%age_univm10) then 
                   gal%disc%sfr210 = gal%disc%sfr210 + mfs
                   gal%disc%sfr21  = mfs
                end if
             end if
             global_sf2 = global_sf2 + mfs
          else
             if (subt1 >= tsno(st)%age_univm100) then 
                gal%disc%sfr100 = gal%disc%sfr100 + mfs
                if (subt1 >= tsno(st)%age_univm10) then 
                   gal%disc%sfr10 = gal%disc%sfr10 + mfs
                   gal%disc%sfr1  = mfs
                end if
             end if
             global_sf = global_sf + mfs
          end if
       end if
#else
       if (mfs > 0.0) then 
          call add_new_stars_to_sfhtab(gal%disc,mfs,metfs)
#ifdef RECORD_SFR
          call add_new_stars_to_sfrtab(gal%disc,mfs,metfs)
#endif
          ! increment SFR counters 
          if (subt1 >= tsno(st)%age_univm100) then 
             gal%disc%sfr100 = gal%disc%sfr100 + mfs 
             if (subt1 >= tsno(st)%age_univm10) then 
                gal%disc%sfr10 = gal%disc%sfr10 + mfs 
                gal%disc%sfr1  = mfs
             end if
          end if
          global_sf = global_sf + mfs 
       endif
#endif
       call get_mcomp(gal%disc) ! update mgal with new minstar and mcold          
       ! update bookkeeping
       call inc_gbk_mstarform(gal%gbk,mfs)


    elseif (comp_type == 'bulg') then 
#ifdef DUAL_IMF
       ! chose an IMF ... 
       use_second_imf = chose_imf(gal%bulge)
#endif          
       !++++++++++++++++++++++++++++ STOP Comment for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++ CALL the EVOLVE_GAS routine for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++

       !call evolve_gas(comp_type,mfs,metfs,eject)
       call evolve_gas(gal%bulge,comp_type,delta_t,mfs,metfs,eject)
       ej_comp  = eject%mcold
       mfs_comp = mfs 

!------------------------------------- STOP ANALYTICAL GAS EVOLUTION -------------------------------------------------

       !call compute_mass_for_stars(gal%bulge,comp_type,delta_t,mfs,metfs)
#ifdef DUAL_IMF
       !call compute_ejecta(gal%bulge,comp_type,eject,mfs,use_second_imf)
#else
       !call compute_ejecta(gal%bulge,comp_type,eject,mfs)
#endif
       if (mfs + eject%mcold + eject%mhot > gal%bulge%mcold) then 
          write(errunit,*) 're-normalization of feedback in bulge... fishy ... '
          mscale       = gal%bulge%mcold / (mfs + eject%mcold + eject%mhot)
          mfs          = mfs * mscale
          metfs        = metfs * mscale
          eject%mcold  = eject%mcold * mscale
          eject%mcoldz = eject%mcoldz * mscale
          eject%mhot   = eject%mhot * mscale
          eject%mhotz  = eject%mhotz * mscale
          if ((mfs + eject%mcold + eject%mhot - gal%bulge%mcold)/gal%bulge%mcold > rel_prec) then 
             write(errunit,*) 'oh no for gas in sf_step... bulge'
             stop
          end if
          if ((metfs + eject%mcoldz + eject%mhotz - gal%bulge%mcoldz)/gal%bulge%mcoldz > rel_prec) then 
             write(errunit,*) 'oh no for metals in sf_step... bulge'
             stop
          end if
       end if

       ! update mcold(z) according to mass of ejected gas and mass of stars formed
       gal%bulge%mcold    = max(0.0d0,gal%bulge%mcold  - eject%mcold  - eject%mhot  - mfs)
       gal%bulge%mcoldz   = max(0.0d0,gal%bulge%mcoldz - eject%mcoldz - eject%mhotz - metfs)
       ! The cold phase of ejecta either mixes with hot halo phase if it exists, or 
       ! they will fall back onto the galaxy. The hot phase will escape the halo if it has no
       ! hot phase already in place, or mix with it otherwise.
! jeje 
       if (h%datas%mhotgaz > h%mfof) then !0.0d0) then 
          print*,'this should never happen... '
          ! The whole wind (hot+cold) mixes with the halo's hot phase
          ejectsave%mhot   = ejectsave%mhot  + eject%mcold  + eject%mhot
          ejectsave%mhotz  = ejectsave%mhotz + eject%mcoldz + eject%mhotz
          call inc_gbk_mej(gal%gbk,0.0d0,eject%mcold+eject%mhot,0.0d0)
       else
          ! the ejected cold gas leaves the galaxy but returns later (to the disc)
          htdyn    = h%datas%rvir / h%datas%cvel * 978.1d0 ! in Gyr
          !tibo
          t1wind   = subt1 + delay_fountain * htdyn + 2.0d0*accretion_deltat
          !t1wind   = subt1 + 0.2 + 2.0d0*accretion_deltat

          t2wind   = t1wind + duration_fountain * htdyn 
          accrate  = eject%mcold / (t2wind - t1wind)
          accratez = eject%mcoldz / (t2wind - t1wind)
          call add_accretion_event(gal%wind,accrate,accratez,t1wind,t2wind)
          ! the ejected hot gas leaves the halo forever
          ejectsave%mout  = ejectsave%mout  + eject%mhot
          ejectsave%moutz = ejectsave%moutz + eject%mhotz
          call inc_gbk_mej(gal%gbk,eject%mcold,0.0d0,eject%mhot)
       end if

#ifdef DUAL_IMF
       if (mfs > 0.0) then 
          call add_new_stars_to_sfhtab(gal%bulge,mfs,metfs,use_second_imf)
#ifdef RECORD_SFR
          call add_new_stars_to_sfrtab(gal%bulge,mfs,metfs,use_second_imf)
#endif
          ! increment SFR counters 
          if (use_second_imf) then 
             if (subt1 >= tsno(st)%age_univm100) then 
                gal%bulge%sfr2100 = gal%bulge%sfr2100 + mfs 
                if (subt1 >= tsno(st)%age_univm10) then 
                   gal%bulge%sfr210 = gal%bulge%sfr210 + mfs
                   gal%bulge%sfr21  = mfs 
                end if
             end if
             global_sf2 = global_sf2 + mfs 
          else
             if (subt1 >= tsno(st)%age_univm100) then 
                gal%bulge%sfr100 = gal%bulge%sfr100 + mfs 
                if (subt1 >= tsno(st)%age_univm10) then 
                   gal%bulge%sfr10 = gal%bulge%sfr10 + mfs 
                   gal%bulge%sfr1  = mfs 
                end if
             end if
             global_sf = global_sf + mfs 
          end if
       end if
#else
       if (mfs > 0.0) then 
          call add_new_stars_to_sfhtab(gal%bulge,mfs,metfs)
#ifdef RECORD_SFR
          call add_new_stars_to_sfrtab(gal%bulge,mfs,metfs)
#endif
          ! increment SFR counters 
          if (subt1 >= tsno(st)%age_univm100) then 
             gal%bulge%sfr100 = gal%bulge%sfr100 + mfs 
             if (subt1 >= tsno(st)%age_univm10) then 
                gal%bulge%sfr10 = gal%bulge%sfr10 + mfs 
                gal%bulge%sfr1  = mfs 
             end if
          end if
          global_sf = global_sf + mfs 
       end if
#endif
       call get_mcomp(gal%bulge) ! update mgal with new minstar and mcold          

       ! update bookkeeping
       call inc_gbk_mstarform(gal%gbk,mfs)
       ! Tibo adds: update bulge bookkeeping
       call inc_gbk_mstarform_bulge(gal%gbk,mfs)
       
    elseif (comp_type == 'burs') then 
#ifdef DUAL_IMF
       ! chose an IMF ... 
       use_second_imf = chose_imf(gal%burst)
#endif          

!++++++++++++++++++++++++++++ STOP Comment for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++ CALL the EVOLVE_GAS routine for analytical gas evolution ++++++++++++++++++++++++++++++++++++++++++++++++++

       !call evolve_gas(comp_type,mfs,metfs,eject)
       call evolve_gas(gal%burst,comp_type,delta_t,mfs,metfs,eject)
       ej_comp  = eject%mcold
       mfs_comp = mfs
!------------------------------------- STOP ANALYTICAL GAS EVOLUTION -------------------------------------------------

       !call compute_mass_for_stars(gal%burst,comp_type,delta_t,mfs,metfs)
#ifdef DUAL_IMF
       !call compute_ejecta(gal%burst,comp_type,eject,mfs,use_second_imf)
#else
       !call compute_ejecta(gal%burst,comp_type,eject,mfs)
#endif
       if (mfs + eject%mcold + eject%mhot > gal%burst%mcold) then 
          write(errunit,*) 're-normalization of feedback in burst... fishy ... '
          mscale       = gal%burst%mcold / (mfs + eject%mcold + eject%mhot)
          mfs          = mfs          * mscale
          metfs        = metfs        * mscale
          eject%mcold  = eject%mcold  * mscale
          eject%mcoldz = eject%mcoldz * mscale
          eject%mhot   = eject%mhot   * mscale
          eject%mhotz  = eject%mhotz  * mscale
          if ((mfs + eject%mcold + eject%mhot - gal%burst%mcold)/gal%burst%mcold > rel_prec) then 
             write(errunit,*) 'oh no for gas in sf_step... burst'
             stop
          end if
          if ((metfs + eject%mcoldz + eject%mhotz - gal%burst%mcoldz)/gal%burst%mcoldz > rel_prec) then 
             write(errunit,*) 'oh no for metals in sf_step... burst '
             stop
          end if
       end if
       ! update mcold(z) according to mass of ejected gas and mass of stars formed
       gal%burst%mcold    = max(0.0d0,gal%burst%mcold  - eject%mcold  - eject%mhot  - mfs)
       gal%burst%mcoldz   = max(0.0d0,gal%burst%mcoldz - eject%mcoldz - eject%mhotz - metfs)
       ! The cold phase of ejecta either mixes with hot halo phase if it exists, or 
       ! they will fall back onto the galaxy. The hot phase will escape the halo if it has no
       ! hot phase already in place, or mix with it otherwise.
! jeje 
       if (h%datas%mhotgaz > h%mfof) then !0.0d0) then 
          print*,'this should never happen... '
          ! The whole wind (hot+cold) mixes with the halo's hot phase
          ejectsave%mhot   = ejectsave%mhot  + eject%mcold  + eject%mhot
          ejectsave%mhotz  = ejectsave%mhotz + eject%mcoldz + eject%mhotz
          call inc_gbk_mej(gal%gbk,0.0d0,eject%mcold+eject%mhot,0.0d0)
       else
          ! the ejected cold gas leaves the galaxy but returns later (to the disc)
          htdyn    = h%datas%rvir / h%datas%cvel * 978.1d0 ! in Gyr
          !tibo
          t1wind   = subt1 + delay_fountain * htdyn + 2.0d0*accretion_deltat
          !t1wind   = subt1 + 0.2 + 2.0d0*accretion_deltat

          t2wind   = t1wind + duration_fountain * htdyn 
          accrate  = eject%mcold / (t2wind - t1wind)
          accratez = eject%mcoldz / (t2wind - t1wind)
          call add_accretion_event(gal%wind,accrate,accratez,t1wind,t2wind)
          ! the ejected hot gas leaves the halo forever
          ejectsave%mout  = ejectsave%mout  + eject%mhot
          ejectsave%moutz = ejectsave%moutz + eject%mhotz
          call inc_gbk_mej(gal%gbk,eject%mcold,0.0d0,eject%mhot)
       end if

#ifdef DUAL_IMF
       if (mfs > 0.0) then 
          call add_new_stars_to_sfhtab(gal%burst,mfs,metfs,use_second_imf)
#ifdef RECORD_SFR 
          call add_new_stars_to_sfrtab(gal%burst,mfs,metfs,use_second_imf)
#endif
          ! increment SFR counters 
          if (use_second_imf) then 
             if (subt1 >= tsno(st)%age_univm100) then 
                gal%burst%sfr2100 = gal%burst%sfr2100 + mfs 
                if (subt1 >= tsno(st)%age_univm10) then 
                   gal%burst%sfr210 = gal%burst%sfr210 + mfs
                   gal%burst%sfr21  = mfs 
                end if
             end if
             global_sf2 = global_sf2 + mfs 
          else
             if (subt1 >= tsno(st)%age_univm100) then 
                gal%burst%sfr100 = gal%burst%sfr100 + mfs 
                if (subt1 >= tsno(st)%age_univm10) then 
                   gal%burst%sfr10 = gal%burst%sfr10 + mfs 
                   gal%burst%sfr1  = mfs
                end if
             end if
             global_sf = global_sf + mfs 
          end if
       end if
#else
       if (mfs > 0.0) then 
          call add_new_stars_to_sfhtab(gal%burst,mfs,metfs)
#ifdef RECORD_SFR
          call add_new_stars_to_sfrtab(gal%burst,mfs,metfs)
#endif
          ! increment SFR counters 
          if (subt1 >= tsno(st)%age_univm100) then 
             gal%burst%sfr100 = gal%burst%sfr100 + mfs 
             if (subt1 >= tsno(st)%age_univm10) then 
                gal%burst%sfr10 = gal%burst%sfr10 + mfs
                gal%burst%sfr1  = mfs 
             end if
          end if
          global_sf = global_sf + mfs 
       end if
#endif
       call get_mcomp(gal%burst) ! update mgal with new minstar and mcold          

       ! update bookkeeping
       call inc_gbk_mstarform(gal%gbk,mfs)


    end if
  
    return

  end subroutine comp_sf_n_fb

!*****************************************************************************************************************
#ifdef RECORD_SFR
  subroutine evolve_sfr_tab(comp,delta_t)

    ! age the sfr_tab content of delta_t

    implicit none

    type(gal_comp)           :: comp      ! component of gal to evolve
    real(kind=8)             :: delta_t   ! during delta_t
    integer(kind=4)          :: i,j,ind_star
    real(kind=8),allocatable :: stars(:,:)
    real(kind=8)             :: frac,totsfr

    if (comp%totsfr > 0.0d0) then 
       ind_star = size(comp%sfr_tab,dim=1)
       allocate(stars(ind_star,nfile))
       stars = comp%sfr_tab
       do i = 1, nfile           ! loop on metallicity cells
          do j = 2, ind_star     ! loop on time cells
             stars(j-1,i)        = stars(j-1,i) * (1.0d0 - delta_t / delta_timetab(j-1))
          end do
       end do
       do i = 1, nfile             ! loop on metallicity cells
          do j = 2, ind_star       ! loop on time cells
             ! stellar gain due to evolution (stars that were in cell j-1 are now in cell j)
             stars(j,i) = stars(j,i) + comp%sfr_tab(j-1,i) * delta_t / delta_timetab(j-1)
          end do
       end do
       comp%sfr_tab = stars
       ! total sfr should not have changed : check that 
       totsfr = sum(comp%sfr_tab)
       if (abs(totsfr - comp%totsfr)/totsfr > rel_prec) then 
          write(errunit,*) '> Warning : sfr changed more than relative precision: '
          write(errunit,*) '> ',abs(totsfr - comp%totsfr)/totsfr, rel_prec,sum(comp%sfr_tab),comp%totsfr
          stop
       end if
       comp%totsfr = totsfr
       deallocate(stars)
    end if

#ifdef DUAL_IMF
    if (comp%totsfr2 > 0.0d0) then 
       ind_star = size(comp%sfr_tab2,dim=1)
       allocate(stars(ind_star,nfile))
       stars = comp%sfr_tab2
       do i = 1, nfile           ! loop on metallicity cells
          do j = 2, ind_star     ! loop on time cells
             stars(j-1,i)        = stars(j-1,i) * (1.0d0 - delta_t / delta_timetab(j-1))
          end do
       end do
       do i = 1, nfile             ! loop on metallicity cells
          do j = 2, ind_star       ! loop on time cells
             ! stellar gain due to evolution (stars that were in cell j-1 are now in cell j)
             stars(j,i) = stars(j,i) + comp%sfr_tab2(j-1,i) * delta_t / delta_timetab(j-1)
          end do
       end do
       comp%sfr_tab2 = stars
       ! total sfr should not have changed : check that 
       totsfr = sum(comp%sfr_tab2)
       if (abs(totsfr - comp%totsfr2)/totsfr > rel_prec) then 
          write(errunit,*) '> Warning : sfr changed more than relative precision: '
          write(errunit,*) '> ',abs(totsfr - comp%totsfr2)/totsfr > rel_prec
          stop
       end if
       comp%totsfr2 = totsfr
       deallocate(stars)
    end if
#endif

    return

  end subroutine evolve_sfr_tab
#endif

  subroutine evolve_gas(comp,comp_type,delta_t,mfs,metfs,ej)

    implicit none
    
    type(gal_comp)          :: comp
    real(kind=8)            :: nhav,coldens,cdthr
    character(4),intent(in) :: comp_type
    real(kind=8)            :: mfs,metfs,mcold,delta_t
    type(ejecta)            :: ej
    real(kind=8),parameter  :: vhot = 650.0d0  ! km/s -> a few 10^6 K
    real(kind=8)            :: vesc,esc_para,feed_fac
    real(kind=8)            :: mej_cold,mej_hot,mej,scale!,dmstar
    real(kind=8)            :: A,B

    ej%mcold  = 0.0d0
    ej%mcoldz = 0.0d0
    ej%mhot   = 0.0d0
    ej%mhotz  = 0.0d0
    
    if (comp_type == 'disc') then 
       esc_para  = disc_params%esc_para
    else 
       esc_para  = bulge_params%esc_para
    endif
    vesc       = root2 * comp%speed * esc_para * 0.3 ! vesc
   
    !  delta_t ~ 9.95e-4  in Gyr, rgal in Mpc and mfs in 10^11 Msun
    
! analytic evolution
    A = (2.0d0 * epsilon * eta_sn * 1.0e8 / 1.989d0 / vesc**2 + 1.0d0) * 0.0328 * alphapar / comp%rgal**0.8 
    B = 2.0d0 * epsilon * eta_sn * 1.0e8 / 1.989d0 / vesc**2

!bug? mcold should be comp%mcold, no? and shouldn't it be after mfs and ej?
!    mcold = ( 0.4 * A * delta_t + 1. / comp%mcold**0.4 )**-2.5 
    
    if (comp_type == 'disc') then 
       ! mass weighted average for column density
       nhav = disc_params%obsc_rad 
    endif
    
    if (comp_type == 'bulg') then       
       nhav = bulge_params%obsc_rad 
    end if

    if (comp_type == 'burs') then
       nhav = bulge_params%obsc_rad
    end if

    if (comp%rgal > 0.0d0) then
       coldens   = log10(comp%mcold/(pi*(nhav*comp%rgal)**2))+18.95d0 ! 10^18.95 --> 10^11/1.4 M_sun/Mpc^2 in at/cm^2  
    else
       coldens   = 0.0d0
    end if

    cdthr            = 20.0d0 ! 20.0d0 = fiducial value  
    if (coldens > cdthr) then ! column density for star formation is high enough 
       mfs   = 1.0d0 / (B + 1.0d0) * (comp%mcold - 1.0d0 / ( 0.4 * A * delta_t + 1.0d0 / comp%mcold**0.4)**2.5)
       metfs = mfs * comp%mcoldz / comp%mcold
       !mfs   = alphapar * 0.0328 * delta_t * comp%mcold**1.4 / comp%rgal**0.8
       !metfs = mfs * comp%mcoldz / comp%mcold
    else
       mfs   = 0.0d0
       metfs = 0.0d0
    end if

! Numeric                                                                                                                                            
    !if (mfs > comp%mcold) then 
    !   write(errunit,*) '> rescaling SFR ... '
    !   mfs   = comp%mcold
    !   metfs = comp%mcoldz
    !end if
    !end Numeric                                                                                                                                          

!analytic evolution    
!    if (coldens > cdthr) then
!       mej = B / (B + 1.0d0) * (comp%mcold - 1.0d0 / ( 0.4 * A * delta_t + 1.0d0 / comp%mcold**0.4)**2.5)
!       comp%mcold  = ( 0.4 * A * delta_t + 1. / comp%mcold**0.4 )**-2.5
!       comp%mcoldz = ( 0.4 * A * delta_t + 1. / comp%mcoldz**0.4 )**-2.5
!    endif

    mej = B * mfs

! Numeric
    !dmstar = mfs
    !feed_fac   = 2.0d0*epsilon*(eta_sn*dmstar*1.0e8/1.989d0)
    !mej = feed_fac / vesc**2
!end Numeric                                                                                                                                          
 
    mej_cold = mej * vhot / (vhot + vesc)
    mej_hot  = mej_cold * vesc / vhot ! momentum equality : mcold vcold = mhot vhot ... 
    
! Numeric                                                                                                                                                      
    !if (mej > comp%mcold) then
    !   write(errunit,*) '> rescaling feedback... '
    !   scale = comp%mcold / mej
    !   mej_cold = mej_cold * scale
    !   mej_hot  = mej_hot  * scale 
    !end if
!end Numeric

    ej%mcold   = mej_cold
    ej%mcoldz  = (mej_cold / comp%mcold) * comp%mcoldz
    ej%mhot    = mej_hot
    ej%mhotz   = (mej_hot / comp%mcold) * comp%mcoldz
    
    
    return
         
  end subroutine evolve_gas

!*****************************************************************************************************************
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

end module HAL_N_GAL_BARYONS
