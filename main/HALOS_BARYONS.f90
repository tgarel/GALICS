module HALOS_BARYONS

  use IO_BARYONS
  use HAL_N_GAL_BARYONS

  public

  ! type used for convenience to pass info from compute_accretion_events to compute_accretion
  type local_accretion
     integer(kind=4)      :: NAccretion_events
     real(kind=8),pointer :: accretion_rate(:),accretion_ratez(:),accretion_start(:),accretion_stop(:)
  end type local_accretion

  contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!*****************************************************************************************************************
  subroutine build_halotree

  ! This subroutine builds the simple halo-tree out of that read in tree.dat
  ! It calls the subroutine 'tree' which does all the job.

    implicit none

    type (halo),pointer :: h
    integer(kind=4)     :: i,st

    do st = 1,nsteps_do                            ! loop on all timesteps
       do i = 1,tsno(st)%nb_of_halos               ! loop on all halos in a given timestep
          h => tsno(st)%liste_halos(i)             ! h is an alias for current halo
          call tree(h)                             ! find progs of h, and define family links. 
#ifndef READ_HALO_ACCRETION
          h%macc = get_halo_macc(h,st)
#endif
       end do
    end do

    return

  end subroutine build_halotree

!*****************************************************************************************************************
  subroutine tree(h) 

  ! This subroutine is called by 'build_halo_tree' in put_baryons. It builds a simpler halo tree from the complete
  ! tree obtained with build_tree.F. Note that %my_fathers refers to the complete tree's fathers,
  ! and %tree%dads refers to the simple-tree's fathers.
  ! The routine proceeds in 3 steps: 
  !   - determine main fathers of h (i.e. halo progenitors which contribute most to the mass of current halo 'h'
  !   (this is done by 'det_main_fathers') and basically removes accretion of background particles)
  !   - check whether h is main son of its main fathers (i.e. it gets most of the dark matter mass they have
  !   (this is done with 'test_cooling'))
  !   - define the simpler tree with the following rules: 
  !      + fathers of h are main fathers and h is their main son 
  !      + these fathers have only one son: h (i.e. if they fragment link to other sons are voided).
  !      + the non-main sons are "castrated" (i.e. they're not allowed to form new galaxies ANYMORE after this)

    implicit none
    
    type (halo)                 :: h              ! this is the current halo which we link to the halo list
    type (halo_ptr),allocatable :: list_h(:)      ! temporary copy of main fathers
    integer(kind=4)             :: st,ierr        ! timestep of h and allocation error variable 
    logical(kind=4)             :: merging        ! do the fathers merge ? 
    logical(kind=4),allocatable :: t(:)           ! will be set in 'test_cooling' (see comment where called)
    integer(kind=4)             :: nbf,nbs        ! number of main fathers, main sons
    real(kind=8)                :: acc_mass       ! quantity of DM mass coming from background accretion  
    integer(kind=4),allocatable :: list_f(:)      ! indices of main fathers
    real(kind=8),allocatable    :: mass_f(:)      ! masses of main fathers
    integer(kind=4)             :: dadno,i        ! number of dads, and dumb loop variable
#ifdef RENUMBERING
    integer(kind=4), allocatable :: renum(:)      ! indices of main fathers
#endif
!RENUMBERING

    st  = h%my_timestep  ! define current timestep

    if (st > 1) then
  
       nbf = h%my_fathers%nb_fathers
       allocate(list_f(nbf),mass_f(nbf))  
#ifdef RENUMBERING  
       allocate(renum(nbf))
#endif
!RENUMBERING


       call det_fathers(h,nbf,list_f,mass_f,merging) 
       ! nbf is now the number of main progenitors (meaning of true halos not background particle accretion)
       ! list_f is the list of these fathers (indices in list of halos tsno(st-1)%liste_halos)
       ! mass_f is the mass transfered from progenitor(s) to halo h
       ! merging is true if at least two fathers merge into h

#ifdef RENUMBERING  
       do i=1,nbf
          renum(i) = list_f(i)
          call renumbering(renum(i),st-1)
       enddo
#endif
!RENUMBERING

       allocate(list_h(nbf),t(nbf),stat=ierr)  
       if (ierr /= 0) then 
          write(errunit,*) '> not enough memory to allocate list_h and t in tree'
       endif
       nbs = 0 
       do i=1,nbf          ! define these temp variables as the fathers - simpler notation
#ifdef RENUMBERING
          list_h(i)%p   => tsno(st-1)%liste_halos(renum(i))
#else
          list_h(i)%p   => tsno(st-1)%liste_halos(list_f(i))
#endif
          ! here t(i) is set to true if h is the main son of list_h(i) (i.e. gets most of its DM parts)
          t(i)          = test_cooling(list_h(i)%p,h)
          ! nbs is the total number of fathers for which h is main son
          if (t(i)) nbs = nbs + 1  
       end do

       ! case 1 --> h is not a main son of any of its fathers (or it has no father)
       if (nbs == 0) then 
          acc_mass = h%mfof
          do i=1,nbf
             acc_mass = acc_mass - mass_f(i)    
          enddo
          ! subcase 1: most of h DM mass comes from accretion: its a genuine new halo
          if (acc_mass/h%mfof >= 0.5d0) then
             h%tree%frag = 0
          else ! subcase 2: most of h DM mass comes from halos: it comes from fragmentation
             h%tree%frag = 1
          endif
       else   
       ! case 2 --> h is a main son of at least one of its fathers
          h%tree%frag = 1
          do i = 1,nbf     
             ! h has fragmented flag if ALL the fathers for which it is main son also have it
             if (t(i)) h%tree%frag = min(list_h(i)%p%tree%frag,h%tree%frag)
          end do
       endif
       
       dadno = 0
       do i = 1,nbf             ! loop on main fathers 
          if (t(i)) then        ! if h is main son of list_h(i), declare it as the only son of list_h(i)
             list_h(i)%p%tree%nsons   = 1  
             allocate(list_h(i)%p%tree%sons(1),stat=ierr)
             if (ierr /= 0) then 
                write(errunit,*) '> not enough memory to allocate list_h%tree%sons in tree'
             endif
             list_h(i)%p%tree%sons(1) = h%my_number
             dadno                    = dadno + 1  ! count the number of real (ie simple tree-based) fathers this halo has.
          end if
       end do

       ! note that previous definition => h can have MORE THAN ONE less father in the simple tree than in the full tree
       h%tree%ndads = dadno
       allocate(h%tree%dads(dadno),stat=ierr)
       if (ierr /= 0) then 
          write(errunit,*) '> not enough memory to allocate list_h%tree%dads in tree'
       endif
       dadno        = 0
       do i = 1,nbf  
          if (t(i)) then        ! if h is main son of list_h(i), declare list_h(i) as one of the fathers of h
             dadno               = dadno + 1 
             h%tree%dads(dadno)  = list_f(i)
          end if
       end do

       deallocate(list_f,mass_f,list_h,t)

#ifdef RENUMBERING
       deallocate(renum)
#endif
!RENUMBERING

    endif  ! if (st > 1)

    return

  end subroutine tree

!*****************************************************************************************************************
  subroutine det_fathers(h,nbf,list_f,mass_f,merging)

  ! This routine determines the progenitors of a given halo 'h'.
  ! Basically, we remove the background from the list of progenitors ...

    implicit none

    type (halo)                 :: h            ! halo for which you want to determine main progs 
    integer(kind=4)             :: nbf,st       ! number of main progs (0,1,2 ...) and timestep of h  
    integer(kind=4)             :: list_f(nbf)  ! list of indices of h's fathers 
    real(kind=8)                :: mass_f(nbf)  ! masses transfered to h from main progs          
    logical(kind=4)             :: merging      ! is there a merging between st-1 and st ?        
    integer(kind=4)             :: j
    real(kind=8),allocatable    :: masses(:)    ! used for sorting
    integer(kind=4),allocatable :: ind(:)       ! used for sorting

    st     = h%my_timestep            ! define current timestep (that of h)
    list_f = 0
    mass_f = 0.d0
    ! now deal with 3 cases : 1 progs, 2 progs, more progs.

    if (h%my_fathers%nb_fathers == 1) then
       ! if there is one father, there are actually two possibilities: 
       !   1/ halo has just formed from accretion of background particles
       !   2/ it has truly one unique dad 
       if (h%my_fathers%list_fathers(1) == 0) then   
          ! number == 0 means that father is environment... i.e. halo h just formed (1st possibility)
          ! so set number of dads to zero
          nbf    = 0
       else
          ! if %number is not 0, then the father is a halo which was well defined at previous timestep
          ! so put it into list_f and mass_f, and set nbf to 1.
          nbf       = 1
          list_f(1) = h%my_fathers%list_fathers(1)
          mass_f(1) = h%my_fathers%mass_fathers(1)*h%mfof*0.01d0 ! go from % to actual mass
       endif
       merging = .false.               ! one single halo just can't merge with itself... 

    else if (h%my_fathers%nb_fathers == 2) then
       ! if there are exactly two fathers which are halos, take both of them into the simple tree. 
       ! it is a merging only if none of the fathers is the environment (otherwise it is accretion)
       if ((h%my_fathers%list_fathers(1) /= 0) .and. (h%my_fathers%list_fathers(2) /= 0)) then
          ! Merging          
          nbf       = 2
          list_f(1) = h%my_fathers%list_fathers(1)
          list_f(2) = h%my_fathers%list_fathers(2)
          mass_f(1) = h%my_fathers%mass_fathers(1)*h%mfof*0.01d0
          mass_f(2) = h%my_fathers%mass_fathers(2)*h%mfof*0.01d0
          merging   = .true.
       else if ((h%my_fathers%list_fathers(1) == 0) .or. (h%my_fathers%list_fathers(2) == 0)) then
          ! Accretion
          merging = .false.
          ! So only one real prog:
          nbf     = 1
          if (h%my_fathers%list_fathers(1) == 0) then
             list_f(1) = h%my_fathers%list_fathers(2)
             mass_f(1) = h%my_fathers%mass_fathers(2)*h%mfof*0.01d0
          else
             list_f(1) = h%my_fathers%list_fathers(1)
             mass_f(1) = h%my_fathers%mass_fathers(1)*h%mfof*0.01d0
          endif
       endif

    else if (h%my_fathers%nb_fathers > 2)  then 
       ! case multiple (>2) merger ...
       ! once again one does not count background accretion.
       nbf = h%my_fathers%nb_fathers
       ! now, sort fathers with the one that contributes most to h's mass first.
       ! this is easy because ..%mass_fathers is percentage of mass given to h.
       allocate(masses(h%my_fathers%nb_fathers),ind(h%my_fathers%nb_fathers))
       do j=1,h%my_fathers%nb_fathers 
          masses(j) = - h%my_fathers%mass_fathers(j) ! NB: masses(j) -ve bc of sub indexx !!
          if (h%my_fathers%list_fathers(j) == 0) then
             masses(j) = 0.0d0    ! get rid of background
             nbf       = nbf - 1
          endif   
       end do
       ! h%my_fathers%list_fathers(ind) is sorted list of fathers in descending order of contributed mass
       call indexx(h%my_fathers%nb_fathers,masses,ind)   
       ! set fathers in list_f
       do j=1,nbf
          list_f(j) = h%my_fathers%list_fathers(ind(j))
          mass_f(j) = h%my_fathers%mass_fathers(ind(j))*h%mfof*0.01d0
       enddo
       ! progs merge into h 
       merging   = .true.
       deallocate(masses,ind)

    endif  

    return

  end subroutine det_fathers

!*****************************************************************************************************************
  function test_cooling(hf,hs)

! checks if halo 'hs' is the main son of halo 'hf', i.e. checks if, among all hf's sons, hs gets most of its 
! DM particles from hf AND most of the DM particles of hf. If it is the case, then hs will also get all the 
! baryons (gas and galaxies) contained in hf. If it ain't, hs don't get nothin'...
! This function therefore proceeds in 4 steps: 
!   1/ determine all brothers of hs (i.e. all sons of hf)
!   2/ among them, select those for which hf is a main progenitor (i.e. hf gives them most of their DM mass)
!   3/ among these, find the one which gets most of hf's dark matter mass
!   4/ if this one is hs, return TRUE, else return FALSE

    implicit none

    logical(kind=4)             :: test_cooling
    type (halo)                 :: hf,hs
    type (halo),pointer         :: hst
    integer(kind=4),allocatable :: liste_brothers(:)                         ! list of all sons of hf
    integer(kind=4),allocatable :: liste_main_brothers(:)                    ! list of main sons of hf, i.e., 
    integer(kind=4)             :: nb_sons,i,inds,ptr,j,nb_main_brothers     ! sons for which hf is actually 
    integer(kind=4)             :: nbf,imaxs,indmaxs,indf,ierr            ! one of the two main progenitors
    integer(kind=4),allocatable :: list_f(:)                                 ! list of fathers
    real(kind=8),allocatable    :: mass_f(:)                                 ! list of fathers' masses 
    real(kind=8)                :: mmaxs,mtemp
    logical(kind=4)             :: merging,main_brother

    allocate(liste_brothers(1:hf%my_sons%nb_sons),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for liste_brothers allocation (in test_cooling)'
       stop
    endif
    allocate(liste_main_brothers(1:hf%my_sons%nb_sons),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for liste_main_brothers allocation (in test_cooling)'
       stop
    endif

    ! 1st step: determine all brothers of hs (ie all sons of hf)
    nb_sons = 0
    do i = 1,hf%my_sons%nb_sons                            ! loop over hf's sons 
       if (hf%my_sons%list_sons(i) == -1) then             ! -1 indicates end of list
          exit
       endif
       inds                    = hf%my_sons%list_sons(i)   ! store son's index
       nb_sons                 = nb_sons+1                 ! increment nb of sons
       liste_brothers(nb_sons) = inds                      ! add son to list
    end do

    ! now, list_brothers contains numbers of all the sons of hf. So 2nd step is to select among these only those
    ! for which hf is a main progenitor.
    ptr = 0  !initialize number of main brothers
    do i = 1,nb_sons    

       inds = liste_brothers(i)
       if (inds == hs%my_number) then
       ! hs is automatically put on the list (since hf is known to be a main prog of hs)
          ptr                      = ptr+1
#ifdef RENUMBERING
          call renumbering(inds,hf%my_timestep+1)
#endif
          liste_main_brothers(ptr) = inds
          cycle 
       endif
#ifdef RENUMBERING
          call renumbering(inds,hf%my_timestep+1)
#endif
       hst => tsno(hf%my_timestep+1)%liste_halos(inds)    ! temporary copy of current son (except hs)
       ! now determine main progs of halo hst:
       nbf = hst%my_fathers%nb_fathers
       allocate(list_f(nbf),mass_f(nbf))
       call det_fathers(hst,nbf,list_f,mass_f,merging)
       !list_f is the list of indices of main progs of hst. If hf among them, set 'main_brother' to true:
       main_brother = .false.
       do j = 1,nbf
          if (list_f(j) == hf%my_number) then
             main_brother = .true.
             exit
          endif
       end do
       deallocate(list_f,mass_f)
      !if we found a main brother, add it to list and increment ptr (ie number of main bros)
       if (main_brother) then 
          ptr                      = ptr+1
          liste_main_brothers(ptr) = inds
          cycle        
       endif

    end do
    nb_main_brothers = ptr
    
    ! now that the list of main sons of hf is established, 3rd step decides where the gas and galaxies will go. 
    ! Recipe used here is: son which gets the most DM particles of hf gets all the baryons.
    mmaxs   = 0.0d0
    imaxs   = 0
    indmaxs = 0
    do i=1,nb_main_brothers

       inds  = liste_main_brothers(i)
       hst   => tsno(hf%my_timestep+1)%liste_halos(inds)
       call det_index_father(hf,hst,indf)
       ! indf is such that hst%my_fathers%list_fathers(indf) == hf%my_number
       ! define mtemp as the mass coming from hf 
       mtemp = hst%my_fathers%mass_fathers(indf)*hst%mfof*0.01d0
       ! as hst%my_fathers%mass_fathers is a mass percentage need to *0.01 to get mass 
       if (mtemp > mmaxs) then  ! store the index of son which gets most of the mass:
          mmaxs   = mtemp
          imaxs   = i
#ifdef RENUMBERING 
          indmaxs = hst%my_number
#else
          indmaxs = inds
#endif
!RENUMBERING
       endif

    end do

    if (indmaxs == hs%my_number) then   ! 4th step: if hs does get the baryons, set 'test_cooling' to true
       test_cooling = .true.
    else                                ! else, hs gets nothing and 'test_cooling' is set to false.
       test_cooling = .false.
    endif

    deallocate(liste_brothers)
    deallocate(liste_main_brothers)

    return

  end function test_cooling

!*****************************************************************************************************************
  subroutine det_index_father(hf,hs,indxf)

  ! finds index of halo hf (father of hs) in the list of progenitors of hs (hs%my_fathers%liste_fathers)
  ! also checks that halo hs is indeed a son of hf.

    implicit none

    type (halo)     :: hf,hs                                           ! halos father and son
    integer(kind=4) :: indxf,i   
    logical(kind=4) :: found_it 

    found_it = .false.
    do i = 1,hs%my_fathers%nb_fathers                                 ! loop on progs of hs
       if (hs%my_fathers%list_fathers(i) == hf%my_number) then
          indxf    = i        ! found it so return the index
          found_it = .true.   
          return
       endif
    end do

    if (.not. found_it) then
       write(errunit,*) '> Problem in det_index_father'
       stop
    end if

    return

  end subroutine det_index_father

!******************************************************************************************************************
  subroutine merge_halos(nbmerg,hs,hf)

  ! this routine does the halo merging: it basically dumps the gas and the galaxies 
  ! of the ith halo progenitor hs(i) into final halo hf for all the nbmerg progenitors 

    implicit none

    integer(kind=4)             :: nbmerg,st,nbf,i,j
    integer(kind=4)             :: nb(nbmerg),nb_cum(0:nbmerg)
    real(kind=8)                :: rn,dt
    real(kind=8)                :: d(nbmerg)
    type (halo)                 :: hs(nbmerg)
    type (halo)                 :: hf
    logical(kind=4)             :: center_occupied
    real(kind=8)                :: rgal
    
    ! hot and cold gas, hot and cold metals accumulate following the fraction of DM passed from hs 
    ! to its son hf.     
    hf%datas%mcoldgaz   = sum(hs(:)%datas%mcoldgaz)
    hf%datas%mcoldz     = sum(hs(:)%datas%mcoldz)
    hf%datas%mhotgaz    = sum(hs(:)%datas%mhotgaz)
    hf%datas%mhotz      = sum(hs(:)%datas%mhotz)
    hf%datas%mgazout    = sum(hs(:)%datas%mgazout)
    hf%datas%metalsout  = sum(hs(:)%datas%metalsout)

    ! now deal with the galaxies
    st             = hf%my_timestep
    nb             = hs(:)%datas%nbgal
    nbf            = sum(nb)
    call alloc_hlist(hf,nbf)

    nb_cum(0)      = 0
    do i=1,nbmerg
        nb_cum(i) = nb_cum(i-1) + nb(i)
    enddo

    do i=1,nbf
        do j=1,nbmerg
           if (nb(j) > 0 .and. i > nb_cum(j-1) .and. i <= nb_cum(j)) then 
              call copy_gal_props(hf%datas%liste_galaxies(i),hs(j)%datas%liste_galaxies(i-nb_cum(j-1)))
           end if
        enddo
    end do
   
    ! now, what radii do galaxies start at?  1st, calc the displacement of the halo CM's 
    center_occupied = .false.
    dt = tsno(st)%age_univ - tsno(hs(1)%my_timestep)%age_univ
    do i=1,nbmerg
        d(i) = cm_disp(hs(i),hf,dt)
    enddo
    do i=1,nbf
       do j=1,nbmerg
          if (nb(j) > 0 .and. i > nb_cum(j-1) .and. i <= nb_cum(j)) then
             rn = starting_radius(hs(j)%datas%liste_galaxies(i-nb_cum(j-1))%r,d(j))
             if (rn == 0.0 .and. center_occupied) then
                ! if there is already a galaxy at the center of the halo, leave it alone and place other candidates
                ! for central position at twice their radius from the center.
                rgal = max(disc_params%s_to_m * hf%datas%liste_galaxies(i)%disc%rgal, &
                     & bulge_params%s_to_m * hf%datas%liste_galaxies(i)%bulge%rgal, &
                     & bulge_params%s_to_m * hf%datas%liste_galaxies(i)%burst%rgal)
                rgal = 2.0d0 * rgal
                hf%datas%liste_galaxies(i)%r   = min(rgal,hf%rfof)
             else
                if (rn == 0.0d0) center_occupied = .true.
                hf%datas%liste_galaxies(i)%r   = min(rn,hf%rfof)
             end if
             call check_orbit_radius(hf%datas%liste_galaxies(i))
          endif   
       end do
    end do

    call get_mhalo(hf)

!!$    ! jeje 
!!$    ! give a merging time to each galaxy (but the central one)
!!$    do i = 1,nbf
!!$       
!!$    end do
!!$    ! end jeje 

    return

    contains

!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    function cm_disp(hi,hf,delta_t)

      ! calculate the sudden displacemnt in halo CM after a merger.  
      ! ie (initial position + velocity * delta_t) - final position
      ! do the calc in comoving co-ords.  to do this, divide the 
      ! physical co-ords by aexp(i).  NB velocities always physical, 
      ! so have to convert these to comving as well.  

      implicit none 

      integer(kind=4) :: k,jj
      real(kind=8)    :: cm_disp,delta_t,f_over_i
      type(halo)      :: hi,hf
      
      f_over_i = tsno(hf%my_timestep)%aexp/tsno(hi%my_timestep)%aexp
      cm_disp  = sqrt((hi%p%x*f_over_i + hi%v%x*delta_t*f_over_i/977.9d0  - hf%p%x)**2 + &
      &               (hi%p%y*f_over_i + hi%v%y*delta_t*f_over_i/977.9d0  - hf%p%y)**2 + &
      &               (hi%p%z*f_over_i + hi%v%z*delta_t*f_over_i/977.9d0  - hf%p%z)**2 )
      do k=1,hf%my_fathers%nb_fathers
         if (hf%my_fathers%list_fathers(k) == hi%my_number) jj = k
      enddo
      cm_disp  = cm_disp * (1- hf%my_fathers%mass_fathers(jj)*0.01d0)

    end function cm_disp

!*******************************************************************************************************************
    function starting_radius(r,d)

      ! given an initial radius, r, and halo CM displacement, d,
      ! select a random position in the original halo and thus 
      ! calc what the new radius is.  

      implicit none

      real(kind=8) :: r,d,x,starting_radius

      if (r == 0.0d0 .and. d == 0.0d0) then 
         starting_radius = 0.0d0
      else     
         ! apply cosine rule, x = cos(theta) has flat distribution between 1 and -1 
         ! for spherical symmetry.  
         call ran1(iseed,x)
         x               = (x - 0.5d0)*2.0d0
         starting_radius = (r**2 + d**2 - 2.0d0*r*d*x)
         starting_radius = sqrt(starting_radius)
!         write(912,'(3(f6.4,1x),2x,3(e10.4,1x),2x,e10.4)') mfof, mfath,percmfath, d,r,starting_radius,rfof
      endif

      return

    end function starting_radius
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  end subroutine merge_halos

!******************************************************************************************************************
  real function time_of_merging(h)
! define time of halo-merging as half time between h's progenitors' time step and h's time step.

    implicit none

    integer(kind=4)         :: st  
    type (halo) ,intent(in) :: h

    st               =  h%my_timestep   
    time_of_merging  =  (tsno(st)%age_univ + tsno(st-1)%age_univ)/2.0d0
    ! could come up with a more complicated recipe, e.g. random number between timesteps.  
    ! in the limit of good time resolution in the simulation, this does not matter.  

    return

  end function time_of_merging

!******************************************************************************************************************

  subroutine compute_accretion(h)

    ! h is not a fragged halo. This routine computes the evolution of h's contents from st-1 to st.

    implicit none

    type(halo)              :: h,htf
    type (halo),allocatable :: list_h(:),ht(:)
    type(local_accretion)   :: acc
    integer(kind=4)         :: st,nbf,igal,renum,i,j,hno_prog,gno_prog
    real(kind=8)            :: macc_hot,i_t(2),tm,total_mass_core
    real(kind=8)            :: error

    ! Some initialisations first 
    ! --------------------------
    st  = h%my_timestep
    nbf = h%tree%ndads
    ! create local copies of h's progenitors
    if (nbf > 0) then 
       allocate(list_h(nbf),ht(nbf))
       do i=1,nbf
          renum = h%tree%dads(i)
#ifdef RENUMBERING
          call renumbering(renum,st-1)
#endif
          call init_halo(list_h(i))
          call copy_halo_props(list_h(i),tsno(st-1)%liste_halos(renum))
          call init_halo(ht(i))
          call copy_halo_props(ht(i),tsno(st-1)%liste_halos(renum))
#ifndef BIG_RUN    
          ! if one halo progenitor is contaminated, so is his son halo 
          ! NB(jeje): this should be done in the tree routine, really ...  
          if (h%ncont == 0) h%ncont =  h%ncont + ht(i)%ncont 
#endif
       end do
    end if
    ! copy of h used after merger or in case nbf < 2
    call init_halo(htf) 
    call copy_halo_props(htf,h)


    ! Now the work : compute accretion and use it
    ! -------------------------------------------
    ! compute accretion budget for h and its possible progs through current timestep
    call compute_accretion_events(h,acc,macc_hot)

    ! evolve h (and/or its progs) from st-1 to st
    select case (nbf)
    case (0) ! new halo -> use htf and copy it back to h at the end

       if (acc%NAccretion_events > 0) then  ! gas indeed gets into h...  
          h%datas%mhotgaz  = macc_hot
          h%datas%mcoldgaz = acc%accretion_rate(1)  * (acc%accretion_stop(1) - acc%accretion_start(1))
          h%datas%mcoldz   = acc%accretion_ratez(1) * (acc%accretion_stop(1) - acc%accretion_start(1))
          h%datas%mgaz     = h%datas%mhotgaz + h%datas%mcoldgaz
          ! create a galaxy
          call alloc_hlist(h,1)
          tsno(st)%n_form                        = tsno(st)%n_form +1 
          h%datas%liste_galaxies(1)%hno          = h%my_number
          h%datas%liste_galaxies(1)%my_number    = 1
          h%datas%liste_galaxies(1)%tgal         = tsno(st)%age_univ
          h%datas%liste_galaxies(1)%tbirth       = tsno(st)%age_univ
          ! define DM core as host halo 
          h%datas%liste_galaxies(1)%core%fdm     = 1.0d0
          h%datas%liste_galaxies(1)%core%rho_0   = h%halo_profile%rho_0
          h%datas%liste_galaxies(1)%core%r_0     = h%halo_profile%r_c
          h%datas%liste_galaxies(1)%core%mass    = h%mfof
          ! add forthcoming accretion event
          ! NB: we smooth accretion over accretion_deltat in add_accretion_event. To avoid mass losses, add
          ! a couple of deltat's delay for new galaxies. (same trick later)
          call add_accretion_event(h%datas%liste_galaxies(1)%acc,acc%accretion_rate(1),acc%accretion_ratez(1),&
               & acc%accretion_start(1)+2.0d0*accretion_deltat,acc%accretion_stop(1)+2.0d0*accretion_deltat)
          ! no evolution to be done because accretion is just reaching the galaxy at ts.
       else ! only hot gas accretion 
          h%datas%mhotgaz  = macc_hot
          h%datas%mgaz     = macc_hot
          h%datas%mcoldgaz = 0.0d0
          h%datas%mcoldz   = 0.0d0
       end if
       
    case (1) ! only one prog
       ! 1/ update accretion events
       ! 1.a/ cold accretion
       if (acc%NAccretion_events > 0) then  ! some cold gas is accreted
          if (list_h(1)%datas%nbgal > 0) then  
             ! if there is a galaxy, update it with new accretion
             call get_central_galaxy(list_h(1),igal)
             call add_accretion_event(list_h(1)%datas%liste_galaxies(igal)%acc,acc%accretion_rate(1),acc%accretion_ratez(1),&
                  & acc%accretion_start(1)+2.0d0*accretion_deltat,acc%accretion_stop(1)+2.0d0*accretion_deltat)
             list_h(1)%datas%mcoldgaz = list_h(1)%datas%mcoldgaz + acc%accretion_rate(1) &
                  & * (acc%accretion_stop(1) - acc%accretion_start(1))
             list_h(1)%datas%mcoldz   = list_h(1)%datas%mcoldz + acc%accretion_ratez(1) &
                  & * (acc%accretion_stop(1) - acc%accretion_start(1))
          else  
             ! if there is no galaxy, create one in list_h(1)
             call alloc_hlist(list_h(1),1)
             tsno(st)%n_form                                = tsno(st)%n_form +1 
             list_h(1)%datas%liste_galaxies(1)%hno          = list_h(1)%my_number
             list_h(1)%datas%liste_galaxies(1)%my_number    = 1
             list_h(1)%datas%liste_galaxies(1)%tgal         = tsno(st-1)%age_univ
             list_h(1)%datas%liste_galaxies(1)%tbirth       = tsno(st-1)%age_univ
             list_h(1)%datas%liste_galaxies(1)%core%fdm     = 1.0d0
             list_h(1)%datas%liste_galaxies(1)%core%rho_0   = h%halo_profile%rho_0
             list_h(1)%datas%liste_galaxies(1)%core%r_0     = h%halo_profile%r_c
             list_h(1)%datas%liste_galaxies(1)%core%mass    = h%mfof
             call add_accretion_event(list_h(1)%datas%liste_galaxies(1)%acc,acc%accretion_rate(1),acc%accretion_ratez(1),&
                  & acc%accretion_start(1)+2.0d0*accretion_deltat,acc%accretion_stop(1)+2.0d0*accretion_deltat)
             list_h(1)%datas%mcoldgaz = list_h(1)%datas%mcoldgaz + acc%accretion_rate(1) &
                  & * (acc%accretion_stop(1) - acc%accretion_start(1))
             list_h(1)%datas%mcoldz   = list_h(1)%datas%mcoldz + acc%accretion_ratez(1) &
                  & * (acc%accretion_stop(1) - acc%accretion_start(1))
          end if
       end if
       ! 1.b/ hot accretion and global bookkeeping 
       list_h(1)%datas%mhotgaz  = list_h(1)%datas%mhotgaz  + macc_hot
       list_h(1)%datas%mgaz     = list_h(1)%datas%mhotgaz  + list_h(1)%datas%mcoldgaz
       
       ! 2/ evolve the halo's galaxies 
       if (list_h(1)%datas%nbgal > 0) then 
          i_t(1) = tsno(st-1)%age_univ
          i_t(2) = tsno(st)%age_univ
          call evolve_halo(list_h(1),h,i_t,st)
       end if
       call modif_halogas_props(list_h(1),h) ! copy gas props of list_h(1) into h

    case default  ! halo merger(s)

       ! I. Evolve all progs from ts-1 to merging time
       ! ---------------------------------------------
       tm     = time_of_merging(h)   ! time of merging
       ! 1/ update accretion events for the progenitors
       if (acc%NAccretion_events > 1) then 
          do i = 1,nbf
             ! cold phase (hot phase added at merging time only).
             if (acc%accretion_rate(i+1) > 0.0d0) then 
                if (list_h(i)%datas%nbgal <= 0) then 
                   write(errunit,'(a)') '> Oh my god, Jeje screwed up ... '
                end if
                call get_central_galaxy(list_h(i),igal)
                call add_accretion_event(list_h(i)%datas%liste_galaxies(igal)%acc,acc%accretion_rate(i+1),&
                     & acc%accretion_ratez(i+1), acc%accretion_start(i+1)+2.0d0*accretion_deltat, &
                     & acc%accretion_stop(i+1)+2.0d0*accretion_deltat)
                list_h(i)%datas%mcoldgaz = list_h(i)%datas%mcoldgaz + acc%accretion_rate(i+1) &
                     & * (acc%accretion_stop(i+1) - acc%accretion_start(i+1))
                list_h(i)%datas%mcoldz   = list_h(i)%datas%mcoldz + acc%accretion_ratez(i+1) &
                     & * (acc%accretion_stop(i+1) - acc%accretion_start(i+1))
             end if
          end do
       end if
       ! 2/ evolve the halo's galaxies
       i_t(1) = tsno(st-1)%age_univ
       i_t(2) = tm
       do i = 1,nbf
          if (list_h(i)%datas%nbgal > 0) then 
             call evolve_halo(list_h(i),ht(i),i_t,st)
          end if
          call modif_halogas_props(list_h(i),ht(i))
       end do
       
       ! II. Merge halos at tm and evolve remnant to ts
       ! ----------------------------------------------
       call merge_halos(nbf,ht,htf)
       call check_orbit_rfof(htf)
       ! 1/ accrete stuff 
       ! add in the hot accretion now... for some arbitrary reason
       htf%datas%mhotgaz = htf%datas%mhotgaz + macc_hot
       if (acc%NAccretion_events > 0) then 
          ! add in the cold flow
          if (acc%accretion_rate(1) > 0.0d0) then
             if (htf%datas%nbgal > 0) then  ! it's likely htf has galaxies
                call get_central_galaxy(htf,igal)
                call add_accretion_event(htf%datas%liste_galaxies(igal)%acc,acc%accretion_rate(1),acc%accretion_ratez(1),&
                     & acc%accretion_start(1)+2.0d0*accretion_deltat,acc%accretion_stop(1)+2.0d0*accretion_deltat)
                htf%datas%mcoldgaz = htf%datas%mcoldgaz + acc%accretion_rate(1)  * (acc%accretion_stop(1) - acc%accretion_start(1))
                htf%datas%mcoldz   = htf%datas%mcoldz   + acc%accretion_ratez(1) * (acc%accretion_stop(1) - acc%accretion_start(1))
             else                           ! but it ain't sure nor necessary... 
                call alloc_hlist(htf,1)
                tsno(st)%n_form                          = tsno(st)%n_form +1 
                htf%datas%liste_galaxies(1)%hno          = htf%my_number
                htf%datas%liste_galaxies(1)%my_number    = 1
                htf%datas%liste_galaxies(1)%tgal         = tm
                htf%datas%liste_galaxies(1)%tbirth       = tm
                htf%datas%liste_galaxies(1)%core%fdm     = 1.0d0
                htf%datas%liste_galaxies(1)%core%rho_0   = h%halo_profile%rho_0
                htf%datas%liste_galaxies(1)%core%r_0     = h%halo_profile%r_c
                htf%datas%liste_galaxies(1)%core%mass    = h%mfof
                call add_accretion_event(htf%datas%liste_galaxies(1)%acc,acc%accretion_rate(1),acc%accretion_ratez(1),&
                     & acc%accretion_start(1)+2.0d0*accretion_deltat,acc%accretion_stop(1)+2.0d0*accretion_deltat)
                htf%datas%mcoldgaz = htf%datas%mcoldgaz + acc%accretion_rate(1)  * (acc%accretion_stop(1) - acc%accretion_start(1))
                htf%datas%mcoldz   = htf%datas%mcoldz   + acc%accretion_ratez(1) * (acc%accretion_stop(1) - acc%accretion_start(1))
             end if
          end if
          ! bookkeep
          htf%datas%mgaz = htf%datas%mcoldgaz + htf%datas%mhotgaz
       end if
       ! 2/ evolve htf into h
       i_t(1) = tm 
       i_t(2) = tsno(st)%age_univ
       if (htf%datas%nbgal > 0) then 
          call evolve_halo(htf,h,i_t,st)
       end if
       call modif_halogas_props(htf,h)

    end select

    ! simple check of baryonic budget
    ! NB: i made some checks with star formation and feedback turned
    ! off. In that case, we should strictly have mcoldgaz = mass in
    ! galaxies (including accretion flows.). However, the mass in
    ! galaxies is incremented during substeps, which leads to small
    ! numerical errors. I left a  factor "nsubsteps" (~100-200) in
    ! sfh_integrator, which ensures better accuracy (but could be made
    ! 1000 for even better results),  and bellow, i use a factor 20
    ! (which is typically a tenth of the number of substeps). This
    ! seems to work, empirically, and avoids  stopping the code for
    ! numerical leaks... 
    if (h%datas%mcoldgaz > 0.0d0) then 
       error = (h%datas%mcoldgaz - cold_gas_in_halo(h))/h%datas%mcoldgaz
    else 
       error = cold_gas_in_halo(h)       
    end if
    if (abs(error) > rel_prec * 1.0e2) then ! tibo: 20 => 1.0e4
!       write(errunit,*) '> error in compute_accretion: mcoldgaz not equal to mass in galaxies ... ',error,rel_prec,h%datas%mcoldgaz
       write(errunit,*) '> error in compute_accretion: mcoldgaz not equal to mass in galaxies ... ',error, cold_gas_in_halo(h),h%datas%mcoldgaz

      ! print*,nbf,h%datas%nbgal
       !print*,h%datas%nbgal,h%datas%mcoldgaz,cold_gas_in_halo(h),total_mass(h%datas%liste_galaxies(1))
       ! stop
       if (error > maxerror) then 
          maxerror = error
       endif
       h%datas%mcoldgaz = cold_gas_in_halo(h)
    else
       h%datas%mcoldgaz = cold_gas_in_halo(h)
    end if

    ! some cleaning 
    ! -------------
    call clear_halo(htf)
    if (acc%naccretion_events > 0) then 
       deallocate(acc%accretion_rate,acc%accretion_ratez,acc%accretion_start,acc%accretion_stop)
       acc%naccretion_events = 0
    end if
    if (nbf > 0) then 
       do i = 1,nbf
          call clear_halo(list_h(i)) 
          call clear_halo(ht(i))
       end do
       deallocate(list_h,ht)
    end if
    
    ! simple check 
    if (h%datas%nbgal /= 0 .and. h%tree%frag == 1) then 
       write(errunit,*) '> error in compute_accretion: fragmented halo contains galaxies',h%my_number
       stop
    endif   


    ! final step : link galaxies to their new haloes
    ! ----------------------------------------------
    if (h%datas%nbgal /= 0) then
       total_mass_core = 0.0d0   
       do i=1,h%datas%nbgal        
          total_mass_core                      = total_mass_core + h%datas%liste_galaxies(i)%core%fdm * &
               & h%datas%liste_galaxies(i)%core%mass
          h%datas%liste_galaxies(i)%hno        = h%my_number
          h%datas%liste_galaxies(i)%my_number  = i
          do j = 1,h%datas%liste_galaxies(i)%my_progs%nb_prog
             hno_prog = h%datas%liste_galaxies(i)%my_progs%hno(j)
#ifdef RENUMBERING
             call renumbering(hno_prog,st-1)
#endif
!RENUMBERING
             gno_prog = h%datas%liste_galaxies(i)%my_progs%gno(j)
             tsno(st-1)%liste_halos(hno_prog)%datas%liste_galaxies(gno_prog)%my_girls%hno = h%my_number
             tsno(st-1)%liste_halos(hno_prog)%datas%liste_galaxies(gno_prog)%my_girls%gno = i
          end do   
       end do
       if (total_mass_core .gt. h%mfof) then 
          write(logerrunit,*) '> Warning sum of core masses > halo mass for halo ',h%my_timestep,h%my_number
          write(logerrunit,*) '> ',total_mass_core,h%mfof,h%datas%nbgal
          write(logerrunit,*) '> fragmentation flag and mass of DM core of central gal'
          write(logerrunit,*) '> ',h%tree%frag,h%datas%liste_galaxies(1)%core%fdm*h%datas%liste_galaxies(1)%core%mass
          write(logerrunit,*) 
       endif   
    end if

    ! normalize SFR counters (so far they are a mass in 10^11Msun integrated over 1,10, and 100Myr -> 
    ! this corresponds to fators 1d5,1d4,1d3 to convert to Msun/Gyr
    do i = 1,h%datas%nbgal
       h%datas%liste_galaxies(i)%disc%sfr100  = h%datas%liste_galaxies(i)%disc%sfr100  * 1.d3 ! Msun/yr
       h%datas%liste_galaxies(i)%disc%sfr10   = h%datas%liste_galaxies(i)%disc%sfr10   * 1.d4 ! Msun/yr
       h%datas%liste_galaxies(i)%disc%sfr1    = h%datas%liste_galaxies(i)%disc%sfr1    * 1.d5 ! Msun/yr
       h%datas%liste_galaxies(i)%bulge%sfr100 = h%datas%liste_galaxies(i)%bulge%sfr100 * 1.d3 ! Msun/yr
       h%datas%liste_galaxies(i)%bulge%sfr10  = h%datas%liste_galaxies(i)%bulge%sfr10  * 1.d4 ! Msun/yr
       h%datas%liste_galaxies(i)%bulge%sfr1   = h%datas%liste_galaxies(i)%bulge%sfr1   * 1.d5 ! Msun/yr
       h%datas%liste_galaxies(i)%burst%sfr100 = h%datas%liste_galaxies(i)%burst%sfr100 * 1.d3 ! Msun/yr
       h%datas%liste_galaxies(i)%burst%sfr10  = h%datas%liste_galaxies(i)%burst%sfr10  * 1.d4 ! Msun/yr
       h%datas%liste_galaxies(i)%burst%sfr1   = h%datas%liste_galaxies(i)%burst%sfr1   * 1.d5 ! Msun/yr
#ifdef DUAL_IMF 
       h%datas%liste_galaxies(i)%disc%sfr2100  = h%datas%liste_galaxies(i)%disc%sfr2100  * 1.d3 ! Msun/yr
       h%datas%liste_galaxies(i)%disc%sfr210   = h%datas%liste_galaxies(i)%disc%sfr210   * 1.d4 ! Msun/yr
       h%datas%liste_galaxies(i)%disc%sfr21    = h%datas%liste_galaxies(i)%disc%sfr21    * 1.d5 ! Msun/yr
       h%datas%liste_galaxies(i)%bulge%sfr2100 = h%datas%liste_galaxies(i)%bulge%sfr2100 * 1.d3 ! Msun/yr
       h%datas%liste_galaxies(i)%bulge%sfr210  = h%datas%liste_galaxies(i)%bulge%sfr210  * 1.d4 ! Msun/yr
       h%datas%liste_galaxies(i)%bulge%sfr21   = h%datas%liste_galaxies(i)%bulge%sfr21   * 1.d5 ! Msun/yr
       h%datas%liste_galaxies(i)%burst%sfr2100 = h%datas%liste_galaxies(i)%burst%sfr2100 * 1.d3 ! Msun/yr
       h%datas%liste_galaxies(i)%burst%sfr210  = h%datas%liste_galaxies(i)%burst%sfr210  * 1.d4 ! Msun/yr
       h%datas%liste_galaxies(i)%burst%sfr21   = h%datas%liste_galaxies(i)%burst%sfr21   * 1.d5 ! Msun/yr
#endif
    end do

    return
    
  end subroutine compute_accretion

!*****************************************************************************************************************

  subroutine compute_accretion_events(h,acc,macc_hot)
    
    ! compute accretion onto h and its progs into acc (acc%acc(1) = onto h after merger, the others are onto progs)
          
    ! give half the accretion to h in first accretion event -> tstart = tm, tstop=ts (+delays)
    ! give a fraction to each prog according to its mass (nbf events) -> tstart=ts-1, tstop=tm (+delays)
    ! -> NB : don't give useless accretion to a prog which doesn't have a galaxy in ... 
    ! -> also consider only non fragged progenitors, etc... (test on containing a galaxy should be ok except for tvir thing)

    implicit none
    
    type(halo)               :: h
    type(local_accretion)    :: acc
    real(kind=8)             :: mbarprog,mbardesc,macctot,macc_hot,macc_cold,tstart,tstop,delay,dt
    real(kind=8)             :: mprogfof,tm,adiab_fact,cfrac,global_rate
    real(kind=8),allocatable :: accfrac(:)
    integer(kind=4)          :: i,st,nbf,renum
!!$    ! jeje 
!!$    real(kind=8)             :: mainmass
!!$    integer(kind=4)          :: imain
!!$    ! end jeje

    ! initialize accretion events
    acc%NAccretion_events = 0
    nullify(acc%accretion_rate,acc%accretion_ratez,acc%accretion_start,acc%accretion_stop)
    macc_hot  = 0.0d0
    macc_cold = 0.0d0
    
#ifndef NO_REIONIZATION    
    ! check that gas actually wants to fall in potential well 
    call compute_adiab_fact(h,adiab_fact)
#else
    adiab_fact = 0.0d0
#endif

    st  = h%my_timestep

    if (h%datas%tvir > tsno(st)%t_igm*adiab_fact .and. &
         & h%tree%frag == 0 .and. h%et < 0.0d0 .and. h%spin <= 0.5d0) then 
   
#ifdef READ_HALO_ACCRETION 
       ! we know the mass that h accreted from the background from ts-1 to ts.
       ! To this mass we add that of "small" accreted haloes (not fragged). These small haloes are
       ! defined as those which do not contain a galaxy (for some reason, which is most
       ! likely that they are below t_igm...). We have to add these explicitly because they 
       ! have no hot phase and no cold phase defined in the model... they do come with baryons, though... 
       macctot  = h%macc * omega_b / omega_0  ! mass accreted from the diffuse background
       nbf      = h%tree%ndads
       mbarprog = 0.0d0
       do i = 1,nbf
          renum = h%tree%dads(i)
#ifdef RENUMBERING
          call renumbering(renum,st-1)
#endif 
          ! count accretion of small haloes as diffuse accretion (spread equaly in time...)
          if (tsno(st-1)%liste_halos(renum)%datas%nbgal == 0) mbarprog = mbarprog + &
               & tsno(st-1)%liste_halos(renum)%mfof * omega_b / omega_0
       end do
       macctot = macctot + mbarprog
#else
       ! In case accretion onto halos is not computed directly from the simulation, make 
       ! an estimate : macc = mass(h) - mass(progs)
       ! NB: take into account all gas that has ever been accreted for mass(h), ie what is there plus mgazout.
       ! compute mbarprog : mass of baryons in the progenitors 
       nbf      = h%tree%ndads
       mbarprog = 0.0d0
       do i=1,nbf 
          renum = h%tree%dads(i)
#ifdef RENUMBERING
          call renumbering(renum,st-1)
#endif
          mbarprog = mbarprog + tsno(st-1)%liste_halos(renum)%datas%mgaz + tsno(st-1)%liste_halos(renum)%datas%mgazout
       end do
       ! compute mbardesc : mass of baryons in h (defined as mfof * fb)
       mbardesc = h%mfof * omega_b / omega_0
       ! mass available for accretion : 
       macctot  = mbardesc - mbarprog
#endif

       if (macctot > 0) then 
          ! Part of this mass is shock-heated, part flows in cold
          ! NB: could do better and interpolate using mass of each prog for half timestep... 
          cfrac     = cold_fraction(h%mfof)
          macc_cold = cfrac * macctot
          macc_hot  = macctot - macc_cold 
       else 
          return ! no gas to accrete : get out of here
       end if
       
       ! distribute this accretion into accretion events for the different progenitors of h
       if (macc_cold <= 0.0d0) return ! return if there is only hot gas accretion 
       select case (nbf)

       case(0)      ! new halo
          
          acc%NAccretion_events = 1
          allocate(acc%accretion_rate(1),acc%accretion_ratez(1),acc%accretion_start(1),acc%accretion_stop(1))
          ! some dynamical time as an indicator of accretion duration and halo formation 
          delay = h%datas%rvir / h%datas%cvel * 978.1d0 ! in Gyr
          !tibo_fast_acc
          delay = delay / 5.d0
          acc%accretion_rate(1)  = macc_cold / delay
          acc%accretion_ratez(1) = macc_cold / delay * ColdAccretionZ
          acc%accretion_start(1) = tsno(st)%age_univ         ! this assumes that the halo formed a dynamical time ago, which cancels 
          acc%accretion_stop(1)  = tsno(st)%age_univ + delay ! out with the dyn time it takes for gas to free fall ... 

       case(1)      ! single prog
          
          acc%NAccretion_events = 1
          allocate(acc%accretion_rate(1),acc%accretion_ratez(1),acc%accretion_start(1),acc%accretion_stop(1))
          ! some dynamical time as an indicator of accretion duration and halo formation 
          delay = h%datas%rvir / h%datas%cvel * 978.1d0 ! in Gyr
          !tibo_fast_acc
          delay = delay / 5.d0
          dt    = tsno(st)%age_univ - tsno(st-1)%age_univ
          acc%accretion_rate(1)  = macc_cold / dt
          acc%accretion_ratez(1) = macc_cold / dt * ColdAccretionZ
          acc%accretion_start(1) = tsno(st-1)%age_univ + delay ! gas takes delay to fall 
          acc%accretion_stop(1)  = tsno(st)%age_univ   + delay

       case default ! mergers

          ! if h has several progs we want to distribute accretion among them, as they probably grow too before merging.
          ! However, we do not want to create a galaxy which will live only half a timestep (ie no accretion for progs which
          ! dont already have a galaxy). Each prog (with a gal) gets part of the accretion budget depending on its mass. 
          allocate(accfrac(nbf))
!!$          ! jeje 
!!$          mainmass = 0.0d0
!!$          imain    = -1
!!$          ! end jeje 
          do i = 1, nbf 
             renum = h%tree%dads(i)
#ifdef RENUMBERING
             call renumbering(renum,st-1)
#endif
             if (tsno(st-1)%liste_halos(renum)%datas%nbgal > 0) then 
                accfrac(i) = tsno(st-1)%liste_halos(renum)%mfof
!!$                ! jeje 
!!$                if (accfrac(i) > mainmass) then 
!!$                   mainmass = accfrac(i)
!!$                   imain    = i
!!$                end if
!!$                ! end jeje
             else
                accfrac(i) = 0.d0
             end if
          end do

!!$          ! jeje  -> give all to one prog.
!!$          accfrac(1:nbf) = 0.0d0
!!$          if (imain >0) accfrac(imain) = 1.0d0
!!$          ! end jeje

          mprogfof = sum(accfrac(:))
          if (mprogfof > 0.0d0) then 
             accfrac = accfrac / mprogfof   ! fraction of mass in each accreting prog (which has a galaxy)
             acc%NAccretion_events = nbf+1
             allocate(acc%accretion_rate(nbf+1),acc%accretion_ratez(nbf+1),acc%accretion_start(nbf+1),&
                  acc%accretion_stop(nbf+1))
             ! accretion rate (the whole macc_cold) happens during a timestep dt : 
             dt          = (tsno(st)%age_univ - tsno(st-1)%age_univ)
             global_rate = macc_cold / dt             

             ! deal with halo h first
             delay    = h%datas%rvir / h%datas%cvel * 978.1d0 ! in Gyr
             !tibo_fast_acc
             delay = delay / 5.d0
             tm       = time_of_merging(h)
             acc%accretion_rate(1)  = global_rate
             acc%accretion_ratez(1) = global_rate * ColdAccretionZ
             acc%accretion_start(1) = tm + delay
             acc%accretion_stop(1)  = tsno(st)%age_univ + delay

             ! deal with progs
             tstart = tsno(st-1)%age_univ
             tstop  = tm 
             do i = 1, nbf
                acc%accretion_rate(i+1)  = accfrac(i) * global_rate
                acc%accretion_ratez(i+1) = ColdAccretionZ * accfrac(i) * global_rate
                renum = h%tree%dads(i)
#ifdef RENUMBERING
                call renumbering(renum,st-1)
#endif
                delay = tsno(st-1)%liste_halos(renum)%datas%rvir / tsno(st-1)%liste_halos(renum)%datas%cvel * 978.1d0 !Gyr
                !tibo_fast_acc
                delay = delay / 5.d0
                acc%accretion_start(i+1) = tstart + delay
                acc%accretion_stop(i+1)  = tstop  + delay
             end do

          else

             ! no progenitor has a galaxy ... All accretion goes to h, starting at time of merging.
             acc%NAccretion_events = 1
             allocate(acc%accretion_rate(1),acc%accretion_ratez(1),acc%accretion_start(1),acc%accretion_stop(1))
             delay = h%datas%rvir / h%datas%cvel * 978.1d0 ! in Gyr
             !tibo_fast_acc
             delay = delay / 5.d0
             tm    = time_of_merging(h)
             dt    = (tsno(st)%age_univ - tm)
             acc%accretion_rate(1)  = macc_cold / dt
             acc%accretion_ratez(1) = ColdAccretionZ * macc_cold / dt
             acc%accretion_start(1) = tm + delay
             acc%accretion_stop(1)  = tsno(st)%age_univ + delay

          end if

       end select
          
    end if
       
    return
    
  end subroutine compute_accretion_events

!*****************************************************************************************************************

  function get_halo_macc(h,st)
    
    implicit none
    
    type(halo)      :: h
    real(kind=8)    :: get_halo_macc,mprog
    integer(kind=4) :: nbf,i,renum,st

    if (st > 1) then 
       nbf   = h%tree%ndads
       mprog = 0.0d0
       do i=1,nbf 
          renum = h%tree%dads(i)
#ifdef RENUMBERING
          call renumbering(renum,st-1)
#endif
          mprog = mprog + tsno(st-1)%liste_halos(renum)%mfof
       end do
    else
       get_halo_macc = h%mfof
    end if
    get_halo_macc = h%mfof - mprog

    return

  end function get_halo_macc

!*****************************************************************************************************************
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

end module HALOS_BARYONS

