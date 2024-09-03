module BOOKKEEP_BARYONS

  use CHECK_BARYONS
  use INIT_BARYONS
  use UTILS_BARYONS
  
  public

  ! new type, used for ease as argument in compute_ejecta
  type ejecta
     real(kind=8) :: mcold      ! gas ejected from ISM in the form of cold clumps
     real(kind=8) :: mcoldz     ! metals ejected from ISM with these cold clups
     real(kind=8) :: mhot       ! gas ejected from the halo in the form of hot phase
     real(kind=8) :: mhotz      ! metals ejected from the halo with this hot phase
  end type ejecta

contains

!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!******************************************************************************************************************
  subroutine stats_st_update(h,st)

    implicit none

    type(halo)      :: h
    integer(kind=4) :: st,j

    ! update timestep statistics
    tsno(st)%mass_bary   = tsno(st)%mass_bary   + h%datas%mgaz
    tsno(st)%mass_halo   = tsno(st)%mass_halo   + h%mfof
    tsno(st)%mass_cold   = tsno(st)%mass_cold   + h%datas%mcoldgaz
    do j = 1,h%datas%nbgal
       tsno(st)%mass_gal_gas               = tsno(st)%mass_gal_gas + total_gal_gas_mass(h%datas%liste_galaxies(j))
       tsno(st)%mass_star                  = tsno(st)%mass_star    + total_stellar_mass(h%datas%liste_galaxies(j))
#ifdef DUAL_IMF
       tsno(st)%mass_star2                 = tsno(st)%mass_star2   + total_stellar_mass2(h%datas%liste_galaxies(j))
#endif
       tsno(st)%mass_metals                = tsno(st)%mass_metals  + total_metal_mass(h%datas%liste_galaxies(j))
       tsno(st)%n_unst                     = tsno(st)%n_unst       + h%datas%liste_galaxies(j)%inst_bulg
       h%datas%liste_galaxies(j)%inst_bulg = 0
    end do
    tsno(st)%n_total     = tsno(st)%n_total     + h%datas%nbgal
    
    ! the only way we can loose galaxies is if their host halo does not have a descendent at the next timestep, i.e.
    ! it is dissolved in the background (goes under the DM mass resolution limit for a halo)
    if (st < nsteps_do .and. h%tree%nsons == 0 .and. h%datas%nbgal > 0) then
       tsno(st+1)%n_lost = tsno(st+1)%n_lost + h%datas%nbgal
       if (h%mfof > 4.d0) tsno(st+1)%n_biglost = tsno(st+1)%n_biglost + h%datas%nbgal 
    endif

    ! count the number of halo merging events in the box 
    tsno(st)%nb_halo_mergs = tsno(st)%nb_halo_mergs + max(0,h%tree%ndads - 1)
    ! mass of hot stuff
    tsno(st)%mass_hot_gas  = tsno(st)%mass_hot_gas  + h%datas%mhotgaz
    tsno(st)%mass_hot_gas  = tsno(st)%mass_hot_mets + h%datas%mhotz

#ifndef BIG_RUN   
    if (h%datas%nbgal > 0) then 
       if (h%ncont > 0) then 
          h%datas%liste_galaxies(1:h%datas%nbgal)%contam = 1
       else
          h%datas%liste_galaxies(1:h%datas%nbgal)%contam = 0
       end if
    endif
#endif

    return

  end subroutine stats_st_update

!******************************************************************************************************************
  subroutine timestep_stats(st)

    implicit none

    integer(kind=4) :: st

    ! finish off with timestep statistics
    tsno(st)%average_tau    = tsno(st)%average_tau / max(tsno(st)%tau_count,1)
    ! global star formation rate in the box
    if (st > 1) then 
       tsno(st)%global_SF   = tsno(st)%global_SF   / (tsno(st)%age_univ - tsno(st-1)%age_univ)
#ifdef DUAL_IMF
       tsno(st)%global_SF2  = tsno(st)%global_SF2  / (tsno(st)%age_univ - tsno(st-1)%age_univ)
#endif
    end if
    ! baryonic fractions
    if (tsno(st)%mass_halo /= 0.0d0) then 
       tsno(st)%mass_bary   = omega_0 * tsno(st)%mass_bary   / tsno(st)%mass_halo
       tsno(st)%mass_cold   = omega_0 * tsno(st)%mass_cold   / tsno(st)%mass_halo
       tsno(st)%mass_star   = omega_0 * tsno(st)%mass_star   / tsno(st)%mass_halo
#ifdef DUAL_IMF
       tsno(st)%mass_star2  = omega_0 * tsno(st)%mass_star2  / tsno(st)%mass_halo
#endif
       tsno(st)%mass_metals = omega_0 * tsno(st)%mass_metals / tsno(st)%mass_halo
    else
       tsno(st)%mass_bary   = omega_b
       tsno(st)%mass_cold   = omega_b*nu
       tsno(st)%mass_star   = 0.0d0
#ifdef DUAL_IMF 
       tsno(st)%mass_star2  = 0.0d0
#endif
       tsno(st)%mass_metals = 0.0d0
    end if

    return

  end subroutine timestep_stats    

!******************************************************************************************************************
  subroutine compute_mass_for_stars(gal,comp_type,delta_t,mfs,metfs)

    implicit none  

    type(gal_comp)          :: gal
    real(kind=8)            :: delta_t,coldensweight,mfs,metfs
    character(4),intent(in) :: comp_type  

    if (comp_type == 'disc') then 
       ! mass weighted average for column density
       coldensweight = disc_params%obsc_rad 
       call comp_sfr(coldensweight,delta_t,gal,'disc',mfs,metfs)
    endif 

    if (comp_type == 'bulg') then       
       coldensweight = bulge_params%obsc_rad 
       call comp_sfr(coldensweight,delta_t,gal,'bulg',mfs,metfs)
    end if

    if (comp_type == 'burs') then
       coldensweight = bulge_params%obsc_rad
       call comp_sfr(coldensweight,delta_t,gal,'burs',mfs,metfs)
    end if
  
    return
  
    contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    subroutine comp_sfr(nhav,delta_t,comp,comp_type,mfs,metfs)

      implicit none

      type(gal_comp)          :: comp
      real(kind=8)            :: nhav,delta_t,coldens,cdthr
      character(4),intent(in) :: comp_type
      real(kind=8)            :: mfs,metfs
      !real(kind=8)            :: metallicity
      !integer(kind=4)         :: met_ind

!!$      ! 21.0: column density threshold for star formation as reported in Kennicutt 1989 in log10(at.cm^2)
!!$      ! 20.0: corresponds to 1/ gal%mcold ~ 10^5 M_sun in a galaxy component of characteristic size gal%rgal ~  0.1 kpc
!!$      !                      2/ gal%mcold ~ 10^7 M_sun in a galaxy component of characteristic size gal%rgal ~  1.0 kpc
!!$      !                      3/ gal%mcold ~ 10^9 M_sun in a galaxy component of characteristic size gal%rgal ~ 10.0 kpc
!!$      !       therefore a negligible fraction of either luminosity or mass of the component and is just here to avoid 
!!$      !       numerical errors / underflows without having an impact on the results ... Probably have to make it bigger
!!$      !       if you want it to have a physical meaning (density threshold for star formation exists in real life)        
!!$      cdthr            = 20.0d0  
!!$
!!$      ! deal with the burst differently than the bulge and disk: no column density threshold rule
!!$      if (comp_type == 'burs') then   
!!$         mfs      = alphapar * comp%mcold  * (delta_t/comp%tdyn) 
!!$         metfs    = alphapar * comp%mcoldz * (delta_t/comp%tdyn)
!!$      else
!!$         coldens       = log10(comp%mcold/(pi*(nhav*comp%rgal)**2))+18.95d0 ! 10^18.95 --> 10^11/1.4 M_sun/Mpc^2 in at/cm^2  
!!$         if (coldens > cdthr) then ! column density for star formation is high enough 
!!$            mfs   = alphapar * comp%mcold  * (delta_t/comp%tdyn)
!!$            metfs = alphapar * comp%mcoldz * (delta_t/comp%tdyn)
!!$         else
!!$            mfs   = 0.0d0
!!$            metfs = 0.0d0
!!$         endif
!!$      end if

      ! new implementation of the Kennicutt law
      ! Simga_sfr (in Msun/yr/kpc^2) = 2.5e-4 * (Sigma_gas [in Msun/pc^2] )^1.4
      ! => SFR (in Msun/yr) = 2.5e-4 / (1000.)^2.8 * Mgas[Msun]^1.4 * r[kpc]^-0.8
      ! where r is the radius to which we average the SFR and gas surface densities ... 
      ! We take r = r_90 which is conviniently 4 * r_scale for both disc and bulge profiles... and also 
      ! should match more or less the scale on which K98 take the average ... 
      ! 
      ! in GALICS units, Sigma_SFR = 0.1 * Sigma_gas^1.4
      ! -> Mstar = 0.1 * delta_t * Mgas^1.4 / R^0.8
      ! -> taking R = 4 * scale_length -> Mstar = 0.0328 * dt * Mgas^1.4 / Rd^0.8
      coldens   = log10(comp%mcold/(pi*(nhav*comp%rgal)**2))+18.95d0 ! 10^18.95 --> 10^11/1.4 M_sun/Mpc^2 in at/cm^2  
      cdthr            = 20.0d0  
      if (coldens > cdthr) then ! column density for star formation is high enough 
         mfs   = alphapar * 0.0328 * delta_t * comp%mcold**1.4 / comp%rgal**0.8
         metfs = mfs * comp%mcoldz / comp%mcold
      else
         mfs   = 0.0d0
         metfs = 0.0d0
      end if
      
      ! add the following lines to be fair with feedback ... 
      if (mfs > comp%mcold) then 
         write(errunit,*) '> rescaling SFR ... '
         mfs   = comp%mcold
         metfs = comp%mcoldz
      end if

      return

  end subroutine comp_sfr

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
end subroutine compute_mass_for_stars

!*****************************************************************************************************************
!!$  subroutine compute_QSO_accretion(gal,delta_t)
!!$
!!$    implicit none  
!!$
!!$    type(galaxy)    :: gal
!!$    real(kind=8)    :: delta_t,t1,t,l_edd,l_qso,l_max,m_eff,m_acc,mmet_acc,tau
!!$    integer(kind=4) :: i
!!$
!!$    ! if there is no bulge, we have a disc with central BH, accrete passively onto it.  
!!$    ! if there is a bulge, we have either a) a starburst or b) very little gas. 
!!$    ! in the latter case, shouldnt matter too much what we do.  In the former, 
!!$    ! pick a time of merging + exponential decay time.  
!!$    if (gamma == 0.0) return
!!$    if (gal%bulge%mgal == 0.0) then 
!!$
!!$       m_acc = gamma * gal%disc%mfs
!!$       m_acc = min(m_acc, gal%disc%mcold)
!!$       if (gal%disc%mcold /= 0.) then
!!$          mmet_acc = (m_acc/gal%disc%mcold) * gal%disc%mcoldz
!!$       else
!!$          mmet_acc = 0.
!!$       endif
!!$       mmet_acc = min(mmet_acc,gal%disc%mcoldz)
!!$
!!$       gal%disc%mcold     = gal%disc%mcold     - m_acc 
!!$       gal%disc%mcoldz    = gal%disc%mcoldz    - mmet_acc
!!$       gal%QSO%mass       = gal%QSO%mass       + m_acc
!!$
!!$       l_edd = mu * gal%QSO%mass* 2.22e11      ! 10^11 M_sun (kms^-1)^2  Gyr^{-1}
!!$       l_qso = mu * m_acc * 8.988e10 / delta_t ! 10^11 M_sun (km s^{-1})^2 Gyr^{-1}       
!!$       m_eff = m_acc*l_edd/max(l_qso,l_edd)    ! if l_qso > l_edd, effective mass is lower than m_acc
!!$       l_qso = min(l_qso,l_edd)
!!$
!!$    else
!!$
!!$       m_acc = gamma * gal%bulge%mfs
!!$       m_acc = min(m_acc,gal%bulge%mcold)
!!$       if (gal%bulge%mcold /= 0.) then
!!$          mmet_acc = (m_acc/gal%bulge%mcold) * gal%bulge%mcoldz
!!$       else
!!$          mmet_acc = 0.
!!$       endif
!!$       mmet_acc = min(mmet_acc,gal%bulge%mcoldz)
!!$
!!$       gal%bulge%mcold     = gal%bulge%mcold     - m_acc 
!!$       gal%bulge%mcoldz    = gal%bulge%mcoldz    - mmet_acc
!!$       gal%QSO%mass        = gal%QSO%mass       + m_acc
!!$
!!$       call ran1(iseed,t)
!!$       t     = t*delta_t ! random number of Gyears ago.  
!!$       l_edd = mu * gal%QSO%mass* 2.22e11
!!$       tau   = t_acc / (gal%QSO%mass / 0.001)**0.2
!!$       l_max = mu * (m_acc * 8.988e10 / tau)
!!$       l_qso = l_max * exp( -t/tau)
!!$       l_qso = min(l_qso,l_edd)
!!$       ! if l_edd < l_max, find the time t1 at which l = l_edd
!!$       if (l_edd < l_max) then 
!!$          t1 = tau * log(l_max/l_edd)
!!$       else 
!!$          t1 = 0.0 
!!$       end if
!!$       m_eff = m_acc*( (t1 * l_edd) + (tau*(exp(-t1/tau) *l_max)))/(tau*l_max) 
!!$
!!$    end if
!!$
!!$    gal%QSO%lum   = l_qso / (6.20e-5)
!!$    gal%QSO%m_eff = m_eff
!!$
!!$    return
!!$  
!!$  end subroutine compute_QSO_accretion

!*****************************************************************************************************************
  subroutine compute_halo_escape_velocity(r,v_esc,h)

    implicit none

    type (halo)  :: h
    real(kind=8) :: r,s,v_esc,rho_0,r_0,r_02,mvir,rvir,phi

    mvir  = h%datas%mvir
    rvir  = h%datas%rvir
    rho_0 = h%halo_profile%rho_0
    r_0   = h%halo_profile%r_c
    r_02  = r_0*r_0
    phi   = 0.0d0

    ! check if the halo is at least partially virialized
    if (rvir > 0.0d0) then 
       if (profile == 'TSIS') then

          ! NB: in the case of a DM TSIS (Truncated Singular Isothermal Sphere), rho_0 is rho_vir and r_0 is r_vir 
          !     as the density in the center is formally infinite  [ rho(r) = rho_vir * r/r_vir ]
          !     and phi as well so for r=0 we use r=1e-4 Mpc=0.1 kpc to be consistent with the routine
          !     compute_pot_gasprof_n_entropy where we apply the same "trick" to phi
          if (r > 0.0d0 .and. r <= rvir) then
             phi = 4*pi*gravconst*rho_0*r_0**2*log(r/r_0) - gravconst*mvir/r_0
          else if (r < 0.0d0) then 
             write(errunit,*) '> error in compute_halo_escape_velocity radius < 0. for halo ',h%my_number 
             stop
          else if (r == 0.0d0) then
             phi = 4*pi*gravconst*rho_0*r_0**2*log(min_size/r_0) - gravconst*mvir/r_0
          else if (r > rvir) then
             ! r is between rvir and rfof ... assume mass outside rvir is negligible
             phi = - gravconst*mvir/h%rfof
          endif

       else if (profile == 'NSTIS') then
          
          if (r > 0.0d0 .and. r <= rvir) then
             ! first compute the constant for the halo potential by equating its expression @ R_vir to that 
             ! of a point mass (the NSTIS is truncated there !!!)
             s          = rvir/r_0
             phi        = - 2.d0*pi*gravconst*rho_0*r_02*(ais*log(1.d0+s*s/a2is) - bis*log(1.d0+s*s/b2is))
             ! now get the absolute halo potential @ r_g            
             s          = r/r_0
             phi        = phi - 4.d0*pi*gravconst*rho_0*r_02/s*((ais-bis)*s - sqrt(a2is)*ais*atan(s/sqrt(a2is))   &
                  & + sqrt(b2is)*bis*atan(s/sqrt(b2is))) + 2.d0*pi*gravconst*rho_0*r_02*(ais*log(1.d0+s*s/a2is) &
                  & - bis*log(1.d0+s*s/b2is))
          else if (r < 0.0d0) then 
             write(errunit,*) '> error in compute_halo_escape_velocity radius < 0. for halo ',h%my_number 
             stop
          else if (r == 0.0d0) then  
             ! phi is actually finite in 0. as atan(r) goes to 0. like r-r**3/3. so following formula is just a 
             ! first order expansion of the r > 0. case above 
             s          = rvir/r_0
             phi        = - 2.d0*pi*gravconst*rho_0*r_02*(ais*log(1.d0+s*s/a2is) - bis*log(1.d0+s*s/b2is))
          else if (r > rvir) then
             ! r is between rvir and rfof ... assume mass outside rvir is negligible
             phi = - gravconst*mvir/h%rfof   
          endif

       else if (profile == 'NFW') then
          
          write(errunit,*) '> error in compute_halo_escape_velocity: profile ',profile,' not implemented yet ' 
          stop
          
       endif

    else  ! halo is not virialized so use potential of homogeneous sphere truncated at rfof to 
          ! compute escape velocity

       phi = -1.d0/2.d0*gravconst*h%mfof/h%rfof*(3.0d0-(r/h%rfof)**2)

    endif

    v_esc  = sqrt(2.0d0*abs(phi)) ! this is the definition of the escape velocity

    return

  end subroutine compute_halo_escape_velocity

!******************************************************************************************************************
  subroutine burst_transport(gal,dt)

    ! move stars older than age_index_transp from burst to bulge (if any such stars)

    implicit none

    type(galaxy)    :: gal
    real(kind=8)    :: dt
    real(kind=8)    :: transp
    integer(kind=4) :: size_bur, size_bul
    real(kind=8),allocatable :: star_temp(:,:)

    if (dt > 0.0d0) then
       if (gal%burst%minstar > 0.0d0) then 
          size_bur = size(gal%burst%sfh_tab,dim=1)
          if (associated(gal%bulge%sfh_tab)) then
             size_bul = size(gal%bulge%sfh_tab,dim=1)
             if (size_bur > size_bul) then
                allocate(star_temp(size_bur,nfile))
                star_temp = 0.0d0
                star_temp(1:size_bul,:) = gal%bulge%sfh_tab
                deallocate(gal%bulge%sfh_tab)
                allocate(gal%bulge%sfh_tab(size_bur,nfile))
                gal%bulge%sfh_tab = star_temp
                deallocate(star_temp)
             end if
          else
             allocate(gal%bulge%sfh_tab(size_bur,nfile))
             gal%bulge%sfh_tab = 0.0d0
          end if
          if (size_bur >= age_index_transp) then
             gal%bulge%sfh_tab(age_index_transp:size_bur,:) = gal%bulge%sfh_tab(age_index_transp:size_bur,:) + &
                  & gal%burst%sfh_tab(age_index_transp:size_bur,:)
             transp                                         = sum(gal%burst%sfh_tab(age_index_transp:size_bur,:))
             allocate(star_temp(age_index_transp-1,nfile))
             star_temp = gal%burst%sfh_tab(1:age_index_transp-1,:)
             deallocate(gal%burst%sfh_tab)
             allocate(gal%burst%sfh_tab(age_index_transp-1,nfile))
             gal%burst%sfh_tab = star_temp
             deallocate(star_temp)
             gal%burst%minstar               = sum(gal%burst%sfh_tab)
             gal%burst%transp                = gal%burst%transp + transp -  gal%burst%minstar ! this is wrong ?
             gal%bulge%minstar               = sum(gal%bulge%sfh_tab)
          end if
       endif
#ifdef DUAL_IMF 
       if (gal%burst%minstar2 > 0.0d0) then 
          size_bur = size(gal%burst%sfh_tab2,dim=1)
          if (associated(gal%bulge%sfh_tab2)) then
             size_bul = size(gal%bulge%sfh_tab2,dim=1)
             if (size_bur > size_bul) then
                allocate(star_temp(size_bur,nfile))
                star_temp = 0.0d0
                star_temp(1:size_bul,:) = gal%bulge%sfh_tab2
                deallocate(gal%bulge%sfh_tab2)
                allocate(gal%bulge%sfh_tab2(size_bur,nfile))
                gal%bulge%sfh_tab2 = star_temp
                deallocate(star_temp)
             end if
          else
             allocate(gal%bulge%sfh_tab2(size_bur,nfile))
             gal%bulge%sfh_tab2 = 0.0d0
          end if
          if (size_bur >= age_index_transp) then
             gal%bulge%sfh_tab2(age_index_transp:size_bur,:) = gal%bulge%sfh_tab2(age_index_transp:size_bur,:) + &
                  & gal%burst%sfh_tab2(age_index_transp:size_bur,:)
             transp                                         = sum(gal%burst%sfh_tab2(age_index_transp:size_bur,:))
             allocate(star_temp(age_index_transp-1,nfile))
             star_temp = gal%burst%sfh_tab2(1:age_index_transp-1,:)
             deallocate(gal%burst%sfh_tab2)
             allocate(gal%burst%sfh_tab2(age_index_transp-1,nfile))
             gal%burst%sfh_tab2 = star_temp
             deallocate(star_temp)
             gal%burst%minstar2               = sum(gal%burst%sfh_tab2)
             gal%burst%transp2                = gal%burst%transp2 + transp -  gal%burst%minstar2 ! this is wrong ?
             gal%bulge%minstar2               = sum(gal%bulge%sfh_tab2)
          end if
       endif
#endif
       call get_mcomp(gal%bulge)
       call get_mcomp(gal%burst)
      
#ifdef RECORD_SFR
       ! repeat same transport for SFR_TAB if they exist.
       ! ------------------------------------------------
       if (gal%burst%totsfr > 0.0d0) then 
          size_bur = size(gal%burst%sfr_tab,dim=1)
          if (associated(gal%bulge%sfr_tab)) then
             size_bul = size(gal%bulge%sfr_tab,dim=1)
             if (size_bur > size_bul) then
                allocate(star_temp(size_bur,nfile))
                star_temp = 0.0d0
                star_temp(1:size_bul,:) = gal%bulge%sfr_tab
                deallocate(gal%bulge%sfr_tab)
                allocate(gal%bulge%sfr_tab(size_bur,nfile))
                gal%bulge%sfr_tab = star_temp
                deallocate(star_temp)
             end if
          else
             allocate(gal%bulge%sfr_tab(size_bur,nfile))
             gal%bulge%sfr_tab = 0.0d0
          end if
          if (size_bur >= age_index_transp) then
             gal%bulge%sfr_tab(age_index_transp:size_bur,:) = gal%bulge%sfr_tab(age_index_transp:size_bur,:) + &
               & gal%burst%sfr_tab(age_index_transp:size_bur,:)
             allocate(star_temp(age_index_transp-1,nfile))
             star_temp = gal%burst%sfr_tab(1:age_index_transp-1,:)
             deallocate(gal%burst%sfr_tab)
             allocate(gal%burst%sfr_tab(age_index_transp-1,nfile))
             gal%burst%sfr_tab = star_temp
             deallocate(star_temp)
             gal%burst%totsfr               = sum(gal%burst%sfr_tab)
             gal%bulge%totsfr               = sum(gal%bulge%sfr_tab)
          end if
       endif
#ifdef DUAL_IMF 
       if (gal%burst%totsfr2 > 0.0d0) then 
          size_bur = size(gal%burst%sfr_tab2,dim=1)
          if (associated(gal%bulge%sfr_tab2)) then
             size_bul = size(gal%bulge%sfr_tab2,dim=1)
             if (size_bur > size_bul) then
                allocate(star_temp(size_bur,nfile))
                star_temp = 0.0d0
                star_temp(1:size_bul,:) = gal%bulge%sfr_tab2
                deallocate(gal%bulge%sfr_tab2)
                allocate(gal%bulge%sfr_tab2(size_bur,nfile))
                gal%bulge%sfr_tab2 = star_temp
                deallocate(star_temp)
             end if
          else
             allocate(gal%bulge%sfr_tab2(size_bur,nfile))
             gal%bulge%sfr_tab2 = 0.0d0
          end if
          if (size_bur >= age_index_transp) then
             gal%bulge%sfr_tab2(age_index_transp:size_bur,:) = gal%bulge%sfr_tab2(age_index_transp:size_bur,:) + &
                  & gal%burst%sfr_tab2(age_index_transp:size_bur,:)
             allocate(star_temp(age_index_transp-1,nfile))
             star_temp = gal%burst%sfr_tab2(1:age_index_transp-1,:)
             deallocate(gal%burst%sfr_tab2)
             allocate(gal%burst%sfr_tab2(age_index_transp-1,nfile))
             gal%burst%sfr_tab2 = star_temp
             deallocate(star_temp)
             gal%burst%totsfr2               = sum(gal%burst%sfr_tab2)
             gal%bulge%totsfr2               = sum(gal%bulge%sfr_tab2)
          end if
       endif
#endif
#endif
    end if
    
    return

  end subroutine burst_transport

!*****************************************************************************************************************
  subroutine compute_core_mass(gal,r_half_mass,core_mass,compute_flag)

    implicit none

    type(galaxy)            :: gal  
    real(kind=8)            :: core_mass,r_half_mass
    character(4),intent(in) :: compute_flag
    real(kind=8)            :: rho_0_cor,r_0_cor,max_mass

    rho_0_cor   = gal%core%rho_0
    r_0_cor     = gal%core%r_0

    if (profile == 'TSIS') then 

       ! rho_0 is actually rho_vir and r_0 is r_vir as the density profile is singular in r=0
       ! we have to use the same "trick" as for the potential i.e. we truncate the density profile @ 0.1 kpc
       ! idem for the satellite: rho_0_sat is rho_vir_sat and r_0_sat is r_vir_sat in the TSIS case
       if (r_half_mass > min_size) then 
          core_mass = 4*pi*rho_0_cor*r_0_cor**2*r_half_mass
       else
          core_mass = 4*pi*rho_0_cor*r_0_cor**2*min_size
       endif
       ! maximal mass in the core (truncation at rvir)
       max_mass = 4*pi*rho_0_cor*r_0_cor**3.d0

    else if (profile == 'NSTIS') then

       core_mass = 4*pi*rho_0_cor*r_0_cor**3*((ais-bis)*r_half_mass    &
               & - sqrt(a2is)*ais*atan(r_half_mass/sqrt(a2is))         &
               & + sqrt(b2is)*bis*atan(r_half_mass/sqrt(b2is)))

    else if (profile == 'NFW') then

       write(errunit,*) '> error in compute_core_mass: profile ',profile,' not implemented yet'
       stop

    endif

    ! if routine called to recompute core mass --> get new core mass 
    if (compute_flag == 'calc') then
       ! total DM core mass is set to twice the amount of DM mass within the half mass radius 
       gal%core%mass = min(2.d0*core_mass,max_mass)
       core_mass     = min(core_mass,max_mass)
    else
       ! DM mass within half mass radius ONLY cannot be greater than what is left of core mass itself
       core_mass     = min(gal%core%mass,core_mass,max_mass)
    endif

    return 


  end subroutine compute_core_mass

!******************************************************************************************************************
  subroutine compute_circular_velocity(r,vc,h)

  ! compute the circular velocity vc of halo h at radius r   

    implicit none 

    type (halo)  :: h
    real(kind=8) :: r,vc,rho_0,r_0,r_02,mvir,rvir

    mvir  = h%datas%mvir
    rvir  = h%datas%rvir
    rho_0 = h%halo_profile%rho_0
    r_0   = h%halo_profile%r_c
    r_02  = r_0*r_0

    ! check if halo is at least partially virialized  
    if (rvir > 0.0d0) then 
       if (profile == 'TSIS') then

          ! NB: in the case of a DM TSIS (Truncated Singular Isothermal Sphere), rho_0 is rho_vir and r_0 is r_vir 
          !     as the density in the center is formally infinite  [ rho(r) = rho_vir * r/r_vir ]
          !     and vc is independent of radius r
          if (r <= rvir) then 
             vc = sqrt(gravconst*mvir/r_0)
          else
             ! radius is between rfof and rvir so interpolate to get vc
             vc = sqrt(gravconst*(mvir+(h%mfof-mvir)*(r-rvir)/h%rfof)/r)
          endif

       else if (profile == 'NSTIS') then

          if (r > 0.0d0 .and. r <= rvir) then 
             vc = sqrt(4*pi*gravconst*rho_0*r_02*(ais-bis-r_0/r*(sqrt(a2is)* &
                  & ais*atan(r/r_0/sqrt(a2is))-sqrt(b2is)*bis*atan(r/r_0/sqrt(b2is))))) 
          else if (r < 0.0d0) then 
             write(errunit,*) '> error in compute_circular_velocity radius < 0. for halo ',h%my_number 
             stop
          else if (r == 0.0d0) then
             ! vc is actually finite in 0. as atan(r) goes to 0. like r-r**3/3. so following formula is just a 
             ! first order expansion of the r > 0. case above
             vc = 0.0d0
          else if (r > rvir) then 
             ! radius is between rfof and rvir so interpolate to get vc
             vc = sqrt(gravconst*(mvir+(h%mfof-mvir)*(r-rvir)/h%rfof)/r)
          endif

       else if (profile == 'NFW') then
       
          write(errunit,*) '> error in compute_circular_velocity: profile ',profile,' not implemented yet ' 
          stop

       endif

    else ! halo not virialized: assume homogeneous sphere profile to compute vc 

       vc = sqrt(gravconst*h%mfof/h%rfof)*r/h%rfof

    endif   

    return

  end subroutine compute_circular_velocity

!******************************************************************************************************************
#ifdef DUAL_IMF 
  subroutine recompute_bulge(gal,h,sta_transf,sfh_tab_trans,gas_transf,met_transf,size_disc,sta_transf2,sfh_tab2_trans,size_disc2)
#else
  subroutine recompute_bulge(gal,h,sta_transf,sfh_tab_trans,gas_transf,met_transf,size_disc)
#endif

    implicit none

    type(galaxy)    :: gal
    type(halo)      :: h
    integer(kind=4) :: size_disc
    real(kind=8)    :: sfh_tab_trans(size_disc,nfile)
    real(kind=8)    :: sta_transf,gas_transf,met_transf
    real(kind=8)    :: rhalf_disc,rhalf_bulge,mtot_bulge_init,e_orbit,mtot_bulge,core_mass
    real(kind=8),allocatable :: star_temp(:,:)
#ifdef DUAL_IMF 
    integer(kind=4) :: size_disc2
    real(kind=8)    :: sfh_tab2_trans(size_disc2,nfile)
    real(kind=8)    :: sta_transf2
#endif
    
    ! transfer star formation history (from disc to burst)
    if (gal%burst%minstar > 0.0d0) then
       if (size(gal%burst%sfh_tab,dim=1) < size_disc) then
          allocate(star_temp(size_disc,nfile)) 
          star_temp = 0.0d0
          star_temp(1:size(gal%burst%sfh_tab,dim=1),:) = gal%burst%sfh_tab
          deallocate(gal%burst%sfh_tab)
          allocate(gal%burst%sfh_tab(size_disc,nfile))
          gal%burst%sfh_tab = star_temp
          deallocate(star_temp)
          gal%burst%sfh_tab = gal%burst%sfh_tab + sfh_tab_trans
       else
          gal%burst%sfh_tab(1:size_disc,:) = gal%burst%sfh_tab(1:size_disc,:) + sfh_tab_trans
       end if
    else
       allocate(gal%burst%sfh_tab(size_disc,nfile))
       gal%burst%sfh_tab = sfh_tab_trans
    end if
#ifdef DUAL_IMF 
    if (gal%burst%minstar2 > 0.0d0) then
       if (size(gal%burst%sfh_tab2,dim=1) < size_disc2) then
          allocate(star_temp(size_disc2,nfile)) 
          star_temp = 0.0d0
          star_temp(1:size(gal%burst%sfh_tab2,dim=1),:) = gal%burst%sfh_tab2
          deallocate(gal%burst%sfh_tab2)
          allocate(gal%burst%sfh_tab2(size_disc2,nfile))
          gal%burst%sfh_tab2 = star_temp
          deallocate(star_temp)
          gal%burst%sfh_tab2 = gal%burst%sfh_tab2 + sfh_tab2_trans
       else
          gal%burst%sfh_tab2(1:size_disc2,:) = gal%burst%sfh_tab2(1:size_disc2,:) + sfh_tab2_trans
       end if
    else
       allocate(gal%burst%sfh_tab2(size_disc2,nfile))
       gal%burst%sfh_tab2 = sfh_tab2_trans
    end if
#endif

    ! compute initial properties of disc and bulge
    rhalf_disc        = gal%disc%rgal     * disc_params%s_to_m
    ! for radius, take max of burst/bulge, because bulge may be empty and will then have rgal==0
    rhalf_bulge       = max(gal%bulge%rgal,gal%burst%rgal)   * bulge_params%s_to_m
    mtot_bulge_init   = gal%bulge%mgal    + gal%burst%mgal

    ! transfer gas and metals from disc to burst
    gal%burst%mcold   = gal%burst%mcold   + gas_transf 
    gal%burst%mcoldz  = gal%burst%mcoldz  + met_transf
    gal%burst%minstar = gal%burst%minstar + sta_transf
#ifdef DUAL_IMF
    gal%burst%minstar2 = gal%burst%minstar2 + sta_transf2
#endif
    
    ! !!! 
    ! transfer stars from disc to bulge : do it yourself (as above)
    ! !!!
    
    call get_mgal(gal,'in recompute_bulge')

    if (mtot_bulge_init > 0.0d0) then ! there already is a bulge/burst

       ! estimate orbital energy of disk material transfered to bulge/burst 
       ! E_orbit = - G M_transf M_inside / (r_half_disc + r_half_bulge) --> no G because it simplifies out later
#ifdef DUAL_IMF 
       e_orbit          = -2.d0* (sta_transf+sta_transf2+gas_transf)*mtot_bulge_init / (rhalf_disc+rhalf_bulge)
#else
       e_orbit          = -2.d0* (sta_transf+gas_transf)*mtot_bulge_init / (rhalf_disc+rhalf_bulge)
#endif
       ! assume all the stuff is initially and finally virialized + conservation of energy 
       ! pot ener    --> 0.4 factor comes from Hernquist profile for bulges (ApJ 1990)
       mtot_bulge       = gal%bulge%mgal + gal%burst%mgal
#ifdef DUAL_IMF
       gal%bulge%rgal   = mtot_bulge*mtot_bulge / (-e_orbit/0.4d0 + mtot_bulge_init*mtot_bulge_init/rhalf_bulge &
                      & + (sta_transf+sta_transf2+gas_transf)**2.d0/rhalf_disc) / bulge_params%s_to_m
#else
       gal%bulge%rgal   = mtot_bulge*mtot_bulge / (-e_orbit/0.4d0 + mtot_bulge_init*mtot_bulge_init/rhalf_bulge &
                      & + (sta_transf+gas_transf)**2.d0/rhalf_disc) / bulge_params%s_to_m
#endif
       ! do not authorize the bulge radius to be greater than the radius of the host halo
       gal%bulge%rgal   = min(gal%bulge%rgal,h%rfof/bulge_params%s_to_m)
       ! limit the minimum bulge size to min_size = 100 pc
       gal%bulge%rgal   = max(gal%bulge%rgal,min_size)
       rhalf_bulge      = gal%bulge%rgal * bulge_params%s_to_m
       call compute_core_mass(gal,rhalf_bulge,core_mass,'nope')
       core_mass        = gal%core%fdm*min(2.d0*core_mass,gal%core%mass)
       ! 3D vel disp @ half mass radius --> 0.177 factor comes also from Hernquist profile for bulges (ApJ 1990)
       gal%bulge%speed  = sqrt(0.177d0*gravconst*(mtot_bulge+core_mass)/gal%bulge%rgal)
       ! size of burst --> proportional to mass (isothermal sphere) 
       gal%burst%rgal   = max(gal%bulge%rgal*gal%burst%mgal/mtot_bulge,min_size) 
       ! velocity of burst is the same as that of bulge --> isothermal sphere approximation
       gal%burst%speed  = gal%bulge%speed
       
    else ! there is no bulge/burst --> have to create one

       mtot_bulge       = gal%burst%mgal + gal%bulge%mgal
       ! assume that the final stuff is virialized with an average 1D velocity dispersion = V_disc/sqrt(2) (isothermal sphere)   
       gal%bulge%rgal   = gravconst*mtot_bulge / (1.5d0*gal%disc%speed*gal%disc%speed) / bulge_params%s_to_m  
       ! do not authorize the bulge radius to be greater than the radius of the host halo
       gal%bulge%rgal   = min(gal%bulge%rgal,h%rfof/bulge_params%s_to_m)
       ! limit the minimum bulge size to min_size = 100 pc
       gal%bulge%rgal   = max(gal%bulge%rgal,min_size)  
       rhalf_bulge      = gal%bulge%rgal * bulge_params%s_to_m
       call compute_core_mass(gal,rhalf_bulge,core_mass,'nope')
       core_mass        = gal%core%fdm*min(2.d0*core_mass,gal%core%mass)
       ! 3D vel disp @ half mass radius --> 0.177 factor comes also from Hernquist profile for bulges (ApJ 1990)
       gal%bulge%speed  = sqrt(0.177d0*gravconst*(mtot_bulge+core_mass)/gal%bulge%rgal)
       ! size of burst according to eq 4.5 of GalICS I
       gal%burst%rgal   = gal%bulge%rgal
       ! velocity of burst is the same as that of bulge
       gal%burst%speed  = gal%bulge%speed
       ! we have just done burst = bulge except mass in bulge is zero not to count mass twice

    endif

    ! compute new dynamical timescales
    gal%bulge%tdyn   = bulge_tdyn(gal%bulge)       
    gal%burst%tdyn   = bulge_tdyn(gal%burst)
    gal%inst_bulg    = 1

    return

  end subroutine recompute_bulge

!*****************************************************************************************************************
  subroutine mass_recipe_burst(A,B,X,Agas,Astar)

    implicit none
    
    real(kind=8)           :: Agas(3,3),Astar(3,3)
    real(kind=8)           :: add,adb,abd,abb,A,B,X

    
    if (B >= 0.0d0) then 
       X = A/B
       if (X >= 1.0d0+rel_prec) then 
          X = 1.0d0/(1.0d0 + ((chi-1.0d0)/(X-1.0d0))**chi)
       else
          X = 0.0d0
       end if
    else
       write(errunit,*) 'Pb in mass_recipe_burst: mgal = 0.0'
       stop
    end if
    
    if (X > 1.d0) stop 'Pb in mass_recipe_burst: X>1'  
    add = X ; abd = 0.0d0 ; adb = 1.0d0 - X ; abb = 1.0d0  

    call check_eps(add,adb)

    ! matrix is (row, column).  columns must add up to 1.  
    ! 0) now:  how much stuff in the disc stays in the disc?
    Agas(1,1) = add ; Astar(1,1) = add
    ! 1) what was in bulge now goes into disc?  assume this is same for gas and for stars, 
    ! and same fraction for bulge as burst    
    Agas(1,2) = 0.0d0 ; Astar(1,2) = 0.0d0    
    Agas(1,3) = 0.0d0 ; Astar(1,3) = 0.0d0
    ! now have to start dealing with things separately!  
    ! 2) all the stars coming from the disk to the burst.
    Astar(2,1) = 0.0d0 
    Astar(3,1) = adb 
    ! 3) all the disk gas also goes to the burst.  
    Agas(2,1) = 0.0d0 
    Agas(3,1) = adb 
    ! 4) now, between the bulge and burst, the stars wont mix: 
    Astar(2,2) = 1.0d0 ; Astar(3,3) = 1.0d0 
    Astar(2,3) = 0.0d0 ; Astar(3,2) = 0.0d0 
    ! 5) but, gas goes to burst preferentially: 
    Agas(2,2) = add     
    Agas(3,2) = 1.0d0 - Agas(2,2) - Agas(1,2)
    ! 6) burst gas stays in burst:
    Agas(3,3) = 1.0d0 ; Agas(2,3) = 0.0d0  

    return

  end subroutine mass_recipe_burst

!*****************************************************************************************************************
  subroutine add_star_props(comp1,comp2,comp3,Astar)

    implicit none

    type(gal_comp)           :: comp1,comp2,comp3
    real(kind=8)             :: Astar(3,3)
    real(kind=8)             :: avec(3),cvec(3)
    real(kind=8),allocatable :: t1(:,:),t2(:,:),t3(:,:) 
    integer(kind=4)          :: max_dim,dim1,dim2,dim3

    avec            = (/comp1%minstar,comp2%minstar,comp3%minstar/) 
    cvec            = matmul(Astar,avec) 
    comp1%minstar   = cvec(1) 
    comp2%minstar   = cvec(2) 
    comp3%minstar   = cvec(3) 

    dim1 = 0 ; dim2 = 0 ; dim3 = 0
    if (associated(comp1%sfh_tab)) dim1 = size(comp1%sfh_tab,dim=1)
    if (associated(comp2%sfh_tab)) dim2 = size(comp2%sfh_tab,dim=1)
    if (associated(comp3%sfh_tab)) dim3 = size(comp3%sfh_tab,dim=1)
    max_dim = max(dim1,dim2,dim3)

    if (max_dim > 0) then 

       allocate(t1(max_dim,nfile),t2(max_dim,nfile),t3(max_dim,nfile))
       t1 = 0.0d0  ;  t2 = 0.0d0  ;  t3 = 0.0d0

       if (dim1 < max_dim) then
          if (dim1 > 0) then
             t1(1:dim1,:) = comp1%sfh_tab
             deallocate(comp1%sfh_tab)
             allocate(comp1%sfh_tab(max_dim,nfile))
          else
             allocate(comp1%sfh_tab(max_dim,nfile))
          end if
       else
          t1 = comp1%sfh_tab
       end if
          
       if (dim2 < max_dim) then
          if (dim2 > 0) then
             t2(1:dim2,:) = comp2%sfh_tab
             deallocate(comp2%sfh_tab)
             allocate(comp2%sfh_tab(max_dim,nfile))
          else
             allocate(comp2%sfh_tab(max_dim,nfile))
          end if
       else
          t2 = comp2%sfh_tab
       end if
          
       if (dim3 < max_dim) then
          if (dim3 > 0) then
             t3(1:dim3,:) = comp3%sfh_tab
             deallocate(comp3%sfh_tab)
             allocate(comp3%sfh_tab(max_dim,nfile))
          else
             allocate(comp3%sfh_tab(max_dim,nfile))
          end if
       else
          t3 = comp3%sfh_tab
       end if
          
       comp1%sfh_tab = t1 * Astar(1,1) + t2 * Astar(1,2) + t3 * Astar(1,3)
       comp2%sfh_tab = t1 * Astar(2,1) + t2 * Astar(2,2) + t3 * Astar(2,3)
       comp3%sfh_tab = t1 * Astar(3,1) + t2 * Astar(3,2) + t3 * Astar(3,3)

       deallocate(t1,t2,t3)
       
    endif

#ifdef DUAL_IMF
    avec            = (/comp1%minstar2,comp2%minstar2,comp3%minstar2/) 
    cvec            = matmul(Astar,avec) 
    comp1%minstar2   = cvec(1) 
    comp2%minstar2   = cvec(2) 
    comp3%minstar2   = cvec(3) 
    dim1 = 0 ; dim2 = 0 ; dim3 = 0
    if (associated(comp1%sfh_tab2)) dim1 = size(comp1%sfh_tab2,dim=1)
    if (associated(comp2%sfh_tab2)) dim2 = size(comp2%sfh_tab2,dim=1)
    if (associated(comp3%sfh_tab2)) dim3 = size(comp3%sfh_tab2,dim=1)
    max_dim = max(dim1,dim2,dim3)
    if (max_dim > 0) then 
       allocate(t1(max_dim,nfile),t2(max_dim,nfile),t3(max_dim,nfile))
       t1 = 0.0d0  ;  t2 = 0.0d0  ;  t3 = 0.0d0
       if (dim1 < max_dim) then
          if (dim1 > 0) then
             t1(1:dim1,:) = comp1%sfh_tab2
             deallocate(comp1%sfh_tab2)
             allocate(comp1%sfh_tab2(max_dim,nfile))
          else
             allocate(comp1%sfh_tab2(max_dim,nfile))
          end if
       else
          t1 = comp1%sfh_tab2
       end if
       if (dim2 < max_dim) then
          if (dim2 > 0) then
             t2(1:dim2,:) = comp2%sfh_tab2
             deallocate(comp2%sfh_tab2)
             allocate(comp2%sfh_tab2(max_dim,nfile))
          else
             allocate(comp2%sfh_tab2(max_dim,nfile))
          end if
       else
          t2 = comp2%sfh_tab2
       end if
       if (dim3 < max_dim) then
          if (dim3 > 0) then
             t3(1:dim3,:) = comp3%sfh_tab2
             deallocate(comp3%sfh_tab2)
             allocate(comp3%sfh_tab2(max_dim,nfile))
          else
             allocate(comp3%sfh_tab2(max_dim,nfile))
          end if
       else
          t3 = comp3%sfh_tab2
       end if
       comp1%sfh_tab2 = t1 * Astar(1,1) + t2 * Astar(1,2) + t3 * Astar(1,3)
       comp2%sfh_tab2 = t1 * Astar(2,1) + t2 * Astar(2,2) + t3 * Astar(2,3)
       comp3%sfh_tab2 = t1 * Astar(3,1) + t2 * Astar(3,2) + t3 * Astar(3,3)
       deallocate(t1,t2,t3)
    endif
#endif
    return

  end subroutine add_star_props

#ifdef RECORD_SFR
!*****************************************************************************************************************
  subroutine add_sfrtab_props(comp1,comp2,comp3,Astar)

    ! "merge" sfr_tabs the same way as sfh_tabs (done in add_star_props)

    implicit none

    type(gal_comp)           :: comp1,comp2,comp3
    real(kind=8)             :: Astar(3,3)
    real(kind=8)             :: avec(3),cvec(3)
    real(kind=8),allocatable :: t1(:,:),t2(:,:),t3(:,:) 
    integer(kind=4)          :: max_dim,dim1,dim2,dim3

    avec            = (/comp1%totsfr,comp2%totsfr,comp3%totsfr/) 
    cvec            = matmul(Astar,avec) 
    comp1%totsfr   = cvec(1) 
    comp2%totsfr   = cvec(2) 
    comp3%totsfr   = cvec(3) 

    dim1 = 0 ; dim2 = 0 ; dim3 = 0
    if (associated(comp1%sfr_tab)) dim1 = size(comp1%sfr_tab,dim=1)
    if (associated(comp2%sfr_tab)) dim2 = size(comp2%sfr_tab,dim=1)
    if (associated(comp3%sfr_tab)) dim3 = size(comp3%sfr_tab,dim=1)
    max_dim = max(dim1,dim2,dim3)

    if (max_dim > 0) then 

       allocate(t1(max_dim,nfile),t2(max_dim,nfile),t3(max_dim,nfile))
       t1 = 0.0d0  ;  t2 = 0.0d0  ;  t3 = 0.0d0

       if (dim1 < max_dim) then
          if (dim1 > 0) then
             t1(1:dim1,:) = comp1%sfr_tab
             deallocate(comp1%sfr_tab)
             allocate(comp1%sfr_tab(max_dim,nfile))
          else
             allocate(comp1%sfr_tab(max_dim,nfile))
          end if
       else
          t1 = comp1%sfr_tab
       end if
          
       if (dim2 < max_dim) then
          if (dim2 > 0) then
             t2(1:dim2,:) = comp2%sfr_tab
             deallocate(comp2%sfr_tab)
             allocate(comp2%sfr_tab(max_dim,nfile))
          else
             allocate(comp2%sfr_tab(max_dim,nfile))
          end if
       else
          t2 = comp2%sfr_tab
       end if
          
       if (dim3 < max_dim) then
          if (dim3 > 0) then
             t3(1:dim3,:) = comp3%sfr_tab
             deallocate(comp3%sfr_tab)
             allocate(comp3%sfr_tab(max_dim,nfile))
          else
             allocate(comp3%sfr_tab(max_dim,nfile))
          end if
       else
          t3 = comp3%sfr_tab
       end if
          
       comp1%sfr_tab = t1 * Astar(1,1) + t2 * Astar(1,2) + t3 * Astar(1,3)
       comp2%sfr_tab = t1 * Astar(2,1) + t2 * Astar(2,2) + t3 * Astar(2,3)
       comp3%sfr_tab = t1 * Astar(3,1) + t2 * Astar(3,2) + t3 * Astar(3,3)

       deallocate(t1,t2,t3)
       
    endif

#ifdef DUAL_IMF
    avec            = (/comp1%totsfr2,comp2%totsfr2,comp3%totsfr2/) 
    cvec            = matmul(Astar,avec) 
    comp1%totsfr2   = cvec(1) 
    comp2%totsfr2   = cvec(2) 
    comp3%totsfr2   = cvec(3) 
    dim1 = 0 ; dim2 = 0 ; dim3 = 0
    if (associated(comp1%sfr_tab2)) dim1 = size(comp1%sfr_tab2,dim=1)
    if (associated(comp2%sfr_tab2)) dim2 = size(comp2%sfr_tab2,dim=1)
    if (associated(comp3%sfr_tab2)) dim3 = size(comp3%sfr_tab2,dim=1)
    max_dim = max(dim1,dim2,dim3)
    if (max_dim > 0) then 
       allocate(t1(max_dim,nfile),t2(max_dim,nfile),t3(max_dim,nfile))
       t1 = 0.0d0  ;  t2 = 0.0d0  ;  t3 = 0.0d0
       if (dim1 < max_dim) then
          if (dim1 > 0) then
             t1(1:dim1,:) = comp1%sfr_tab2
             deallocate(comp1%sfr_tab2)
             allocate(comp1%sfr_tab2(max_dim,nfile))
          else
             allocate(comp1%sfr_tab2(max_dim,nfile))
          end if
       else
          t1 = comp1%sfr_tab2
       end if
       if (dim2 < max_dim) then
          if (dim2 > 0) then
             t2(1:dim2,:) = comp2%sfr_tab2
             deallocate(comp2%sfr_tab2)
             allocate(comp2%sfr_tab2(max_dim,nfile))
          else
             allocate(comp2%sfr_tab2(max_dim,nfile))
          end if
       else
          t2 = comp2%sfr_tab2
       end if
       if (dim3 < max_dim) then
          if (dim3 > 0) then
             t3(1:dim3,:) = comp3%sfr_tab2
             deallocate(comp3%sfr_tab2)
             allocate(comp3%sfr_tab2(max_dim,nfile))
          else
             allocate(comp3%sfr_tab2(max_dim,nfile))
          end if
       else
          t3 = comp3%sfr_tab2
       end if
       comp1%sfr_tab2 = t1 * Astar(1,1) + t2 * Astar(1,2) + t3 * Astar(1,3)
       comp2%sfr_tab2 = t1 * Astar(2,1) + t2 * Astar(2,2) + t3 * Astar(2,3)
       comp3%sfr_tab2 = t1 * Astar(3,1) + t2 * Astar(3,2) + t3 * Astar(3,3)
       deallocate(t1,t2,t3)
    endif
#endif
    return

  end subroutine add_sfrtab_props
#endif
!*****************************************************************************************************************
  subroutine add_gas_props(a1,a2,a3,Agas)

    implicit none

    real(kind=8) :: Agas(3,3)
    real(kind=8) :: avec(3),cvec(3)
    real(kind=8) :: a1,a2,a3

    avec = (/a1,a2,a3/)
    cvec = matmul(Agas,avec) 
    a1   = cvec(1)
    a2   = cvec(2) 
    a3   = cvec(3) 

    return

  end subroutine add_gas_props

!*****************************************************************************************************************
  subroutine increment_mcold(g1,g2)

    implicit none

    type(galaxy),intent(in)    :: g2
    type(galaxy),intent(inout) :: g1

    g1%disc%mcold       = g1%disc%mcold  + g2%disc%mcold    
    g1%bulge%mcold      = g1%bulge%mcold + g2%bulge%mcold   
    g1%burst%mcold      = g1%burst%mcold + g2%burst%mcold    

    return

  end subroutine increment_mcold

!******************************************************************************************************************
  subroutine increment_mcoldz(g1,g2)

    implicit none

    type(galaxy),intent(in)    :: g2
    type(galaxy),intent(inout) :: g1

    g1%disc%mcoldz       = g1%disc%mcoldz  + g2%disc%mcoldz
    g1%bulge%mcoldz      = g1%bulge%mcoldz + g2%bulge%mcoldz
    g1%burst%mcoldz      = g1%burst%mcoldz + g2%burst%mcoldz

    return

  end subroutine increment_mcoldz

!******************************************************************************************************************
  subroutine increment_minstar(g1,g2)

    ! perform g1 = g1 + g2 ... 

    implicit none

    type(galaxy),intent(in)    :: g2
    type(galaxy),intent(inout) :: g1
    integer(kind=4)            :: size1,size2
    real(kind=8),allocatable   :: star_temp(:,:)

    g1%disc%sfr1   = g1%disc%sfr1   + g2%disc%sfr1
    g1%disc%sfr10  = g1%disc%sfr10  + g2%disc%sfr10
    g1%disc%sfr100 = g1%disc%sfr100 + g2%disc%sfr100
#ifdef DUAL_IMF
    g1%disc%sfr21   = g1%disc%sfr21   + g2%disc%sfr21
    g1%disc%sfr210  = g1%disc%sfr210  + g2%disc%sfr210
    g1%disc%sfr2100 = g1%disc%sfr2100 + g2%disc%sfr2100
#endif
    
#ifdef RECORD_SFR
    call add_sfr_tabs(g1%disc,g2%disc)
    call add_sfr_tabs(g1%bulge,g2%bulge)
    call add_sfr_tabs(g1%burst,g2%burst)
#endif
    
    call add_sfh_tabs(g1%disc,g2%disc)
    call add_sfh_tabs(g1%bulge,g2%bulge)
    call add_sfh_tabs(g1%burst,g2%burst)
    
    return
    
  end subroutine increment_minstar
  
!******************************************************************************************************************
#ifdef DUAL_IMF
  subroutine star_gas_transf(mtot_transf,mstars1,mstars2,mgas)
#else
  subroutine star_gas_transf(mtot_transf,mstars,mgas)
#endif

  ! We want to remove a mass mtot_transf from the disc and put it into a buxxx component. 
  ! This routine computes the *fraction* of stars and gas to remove from the disc, assuming that 
  ! we move stars and gas according to the ratio in the disc.
  ! In case of DUAL_IMF, we also want to move things preserving the population ratio.

  ! USAGE : 
  ! mstars(1,2) and mgas come in as the stellar and gas content of the disc and come out as the mass to remove.

    implicit none

    real(kind=8),intent(in)    :: mtot_transf
    real(kind=8)               :: mstars,mgas
#ifdef DUAL_IMF
    real(kind=8)               :: mstars1,mstars2
    real(kind=8)               :: frac1
#endif
    real(kind=8)               :: star_frac    ! if possible, transfer this % (in mass) of stars (and the rest of gas)
    real(kind=8)               :: m

#ifdef DUAL_IMF
    mstars = mstars1 + mstars2
#endif
    
    if (mtot_transf > mstars + mgas) then
       stop 'mtot_transf > mstars + mgas in BOOKKEEP_BARYONS:star_gas_transf '
    end if

    ! define mass percentage of stars to transfer 
    star_frac = mstars / (mgas + mstars)
        
    ! mass of stars to remove : 
    m    = mtot_transf * star_frac
#ifdef DUAL_IMF 
    frac1   = mstars1 / (mstars1 + mstars2) 
    mstars1 = min(m * frac1,mstars1)
    mstars2 = min(m * (1.d0 - frac1),mstars2)
#else
    mstars  = min(m,mstars)
#endif

    ! mass of gas to remove : 
    m    = mtot_transf  - m
    mgas = min(m,mgas)
    
    return

  end subroutine star_gas_transf
  
!******************************************************************************************************************
#ifdef DUAL_IMF 
  subroutine compute_ejecta(comp,comp_type,ej,dmstar,use_second_imf)
#else
  subroutine compute_ejecta(comp,comp_type,ej,dmstar)
#endif

    ! We use a new feedback model, based on (1) Dekel & Silk (86) energy argument, (2) momentum equipartition, and (3) observations. 
    ! It goes like this : 
    ! A wind has two phases : cold and hot. The hot phase gets energy from supernovae and gives momentum (till equipartition) to the cold.
    ! -> Mcold Vcold = Mhot Vhot
    ! The total energy Etot is the sum of energy in cold and hot phases. We define epsilon_cold such that Ecold = epsilon_cold * Etot 
    ! (idem for hot). Thus, epsilon_hot = Vhot / (Vhot + Vcold) and epsilon_cold = Vcold / (Vcold + Vhot)
    ! Now, Etot is Esn * eta_sn * SFR * epsilon_sn. 
    ! It follows that mcold = 2 epsilon_sn Esn eta_sn SFR / Vcold / (Vhot + Vcold)
    ! all the same mhot = 2 epsilon_sn Esn eta_sn SFR / Vhot / (Vhot + Vcold)
    ! The total mass ejected from the galaxy is then Mej = mhot + mcold = 2 epsilon_sn Esn eta_sn SFR / VhotVcold
    ! We now need to set Vhot and Vcold to some values. From observations, we decide : Vhot ~ 650km/s (independent of galaxy), 
    ! and Vcold ~ Vesc (= root2 * gal%speed * esc_para)
    ! The cold ejecta will fall back (unless they evaporate in massive haloes) and the hot ejecta will leave the halo. 

    implicit none

    type(gal_comp)          :: comp
    type(ejecta)            :: ej
    character(4),intent(in) :: comp_type  
    real(kind=8),parameter  :: vhot = 650.0d0  ! km/s -> a few 10^6 K
    real(kind=8)            :: vesc,esc_para,feed_fac,dmstar
    real(kind=8)            :: mcold,mhot,mej,scale
#ifdef DUAL_IMF
    logical(kind=4)         :: use_second_imf
#endif

    ej%mcold  = 0.0d0
    ej%mcoldz = 0.0d0
    ej%mhot   = 0.0d0
    ej%mhotz  = 0.0d0

#ifdef NO_SN_FEEDBACK
    return
#endif

    if (comp_type == 'disc') then 
       esc_para  = disc_params%esc_para
    else 
       esc_para  = bulge_params%esc_para
    endif 
    vesc       = root2 * comp%speed * esc_para * 0.3

    ! Esn = 10^51 ergs is 10^8/1.989 Msun (kms^-1)^2.  this cancels nicely with eta_sn and 1/v^2.
#ifdef DUAL_IMF
    if (use_second_imf) then 
       feed_fac   = 2.0d0*epsilon*(eta_sn2*dmstar*1.0e8/1.989d0) 
    else
       feed_fac   = 2.0d0*epsilon*(eta_sn*dmstar*1.0e8/1.989d0) 
    end if
#else
    feed_fac   = 2.0d0*epsilon*(eta_sn*dmstar*1.0e8/1.989d0) 
#endif
    
!!$    mcold = feed_fac / vesc / (vhot + vesc)
!!$    mhot  = feed_fac / vhot / (vhot + vesc)
!!$    mej   = mcold + mhot

    ! new attemp
    mej = feed_fac / vesc**2 
    mcold = mej * vhot / (vhot + vesc)
    mhot  = mcold * vesc / vhot ! momentum equality : mcold vcold = mhot vhot ... 

    if (mej > comp%mcold) then
       write(errunit,*) '> rescaling feedback... '
       scale = comp%mcold / mej
       mcold = mcold * scale
       mhot  = mhot  * scale 
    end if
    ej%mcold   = mcold
    ej%mcoldz  = (mcold / comp%mcold) * comp%mcoldz
    ej%mhot    = mhot
    ej%mhotz   = (mhot / comp%mcold) * comp%mcoldz

    return
    
  end subroutine compute_ejecta
!******************************************************************************************************************
  
  subroutine disc_dynamics(gal,macc,h)
    
    ! compute/update dynamical properties of a disc. 

    implicit none

    type(galaxy) :: gal
    type(halo)   :: h
    real(kind=8) :: macc    ! mass of gas to be accreted onto the disc
    real(kind=8) :: r_half_mass,rhalfmgal,core_mass,mbulge_in_disc,mbar_in_disc
    real(kind=8) :: mdm_in_disc,fdm
    real(kind=8) :: newrad, fact
    real(kind=8),parameter   :: speed_crit = 500.0d0 
 
    if (macc + gal%disc%mgal <= 0.0d0) return ! there's no disc here and there ain't gonna be any soon

    ! 1/ update disc size in case of accretion
    ! ----------------------------------------
    if (macc > 0.0d0) then 
       newrad     = h%spin * h%datas%rvir / root2
       newrad     = max(newrad,min_size)   ! limit minimum size to 100pc for extremely low-spin haloes
       if (gal%disc%mgal > 0.0d0) then  ! the disc is not new
!!$          fact          = (macc + gal%disc%mcold) / (macc * gal%disc%rgal**2 + gal%disc%mcold * newrad**2)
          fact          = (macc + gal%disc%mgal) / (macc * gal%disc%rgal**2 + gal%disc%mgal * newrad**2)
          gal%disc%rgal = gal%disc%rgal * newrad * sqrt(fact)
       else                           ! this disc is new
          gal%disc%rgal = newrad
       end if
    end if
    
    ! jeje : the following line has no real physical sense : i remove it ... 
    ! do not allow the half mass radius of the disc to become smaller than that of the bulge
!!$    gal%disc%rgal     = max(gal%disc%rgal,gal%bulge%rgal*bulge_params%s_to_m/disc_params%s_to_m)
    ! do not let half mass radius of disc get bigger than halo radius 
    gal%disc%rgal     = min(gal%disc%rgal,h%rfof/disc_params%s_to_m)

    ! 2/ update dynamics (speed and dynamical time)
    ! ---------------------------------------------
    r_half_mass       = gal%disc%rgal*disc_params%s_to_m
    if (gal%r <= min_size) then
       ! recompute DM core properties assuming it extends at max inside 2*r_1/2 
       gal%core%rho_0  = h%halo_profile%rho_0 
       gal%core%r_0    = h%halo_profile%r_c
       rhalfmgal       = (macc*r_half_mass + (gal%bulge%mgal*gal%bulge%rgal+gal%burst%mgal*gal%burst%rgal) &
            &           *bulge_params%s_to_m)/(total_mass(gal)+macc)
       call compute_core_mass(gal,rhalfmgal,core_mass,'calc')
    endif
    mbulge_in_disc  = (gal%bulge%mgal+gal%burst%mgal)*(r_half_mass/(r_half_mass+gal%bulge%rgal))**2
    mbar_in_disc    = 0.5d0*(gal%disc%mgal+macc) + mbulge_in_disc
    call compute_core_mass(gal,r_half_mass,core_mass,'nope')
    mdm_in_disc     = core_mass
    gal%disc%speed  = sqrt(gravconst*(mdm_in_disc+mbar_in_disc)/r_half_mass)
    ! fraction of dark matter mass remaining in half mass radius of disc 
    fdm             = max(1.d0-(gal%disc%speed/speed_crit)**2,0.0d0)
    ! update DM core to be consistent
    gal%core%fdm    = fdm
    mdm_in_disc     = fdm*mdm_in_disc
    gal%disc%speed  = sqrt(gravconst*(mdm_in_disc+mbar_in_disc)/r_half_mass)
    gal%disc%tdyn   = disc_tdyn(gal%disc)

    return 

  end subroutine disc_dynamics

!******************************************************************************************************************
  subroutine check_for_disc_instability(h,gal)
    
    implicit none
    
    type(halo)               :: h
    type(galaxy)             :: gal
    real(kind=8)             :: alpha_crit,r_half_mass,core_mass,mbulge_in_disc,mdm_in_disc
    real(kind=8)             :: tot_transf,sta_transf,gas_transf,met_transf
    real(kind=8),allocatable :: sfh_tab_trans(:,:)
    integer(kind=4)          :: size_disc    
    real(kind=8)             :: macc
#ifdef RECORD_SFR  
    real(kind=8)             :: sfr_transf
#ifdef DUAL_IMF
    real(kind=8)             :: sfr_transf2
#endif
#endif
#ifdef DUAL_IMF
    real(kind=8)             :: sta_transf2
    real(kind=8),allocatable :: sfh_tab2_trans(:,:)
    integer(kind=4)          :: size_disc2
#endif
    
    alpha_crit = 8.d0*disc_instability_threshhold**2.d0/(1.d0-4.d0*disc_instability_threshhold**2.d0)

    r_half_mass     = gal%disc%rgal*disc_params%s_to_m
    call compute_core_mass(gal,r_half_mass,core_mass,'nope')
    mbulge_in_disc  = (gal%burst%mgal+gal%bulge%mgal)*(r_half_mass/(r_half_mass+gal%bulge%rgal))**2
    ! fraction of dark matter mass remaining in half mass radius of disc 
    mdm_in_disc     = gal%core%fdm*core_mass

    ! if disc is more than meta_stable make it marginally stable ... 
    if (gal%disc%mgal .gt. meta_stable_fudge*alpha_crit*(mbulge_in_disc+mdm_in_disc)) then
       ! total mass to transfer to burst/bulge in order to make the disk stable
       tot_transf = (gal%disc%mgal-alpha_crit*(mbulge_in_disc+mdm_in_disc))/(1.d0+alpha_crit)
       sta_transf = gal%disc%minstar
       gas_transf = gal%disc%mcold
#ifdef DUAL_IMF 
       sta_transf2 = gal%disc%minstar2
       call star_gas_transf(tot_transf,sta_transf,sta_transf2,gas_transf)
#else
       call star_gas_transf(tot_transf,sta_transf,gas_transf)
#endif
       ! remove stars from disc
       if (sta_transf > 0.0d0) then 
#ifdef RECORD_SFR 
          sfr_transf = sta_transf / gal%disc%minstar
#endif
          size_disc        = size(gal%disc%sfh_tab,dim=1)
          allocate(sfh_tab_trans(size_disc,nfile))
          sfh_tab_trans    = 0.0d0
          sfh_tab_trans    = gal%disc%sfh_tab *  sta_transf / gal%disc%minstar
          gal%disc%sfh_tab = gal%disc%sfh_tab - sfh_tab_trans                
          gal%disc%minstar = sum(gal%disc%sfh_tab)
          sta_transf       = sum(sfh_tab_trans)
          gal%disc%transp  = gal%disc%transp  + sta_transf
       else
          size_disc        = 1
          allocate(sfh_tab_trans(size_disc,nfile))
          sfh_tab_trans    = 0.0d0
       end if
#ifdef DUAL_IMF 
       if (sta_transf2 > 0.0d0) then 
#ifdef RECORD_SFR 
          sfr_transf2 = sta_transf2 / gal%disc%minstar2
#endif
          size_disc2        = size(gal%disc%sfh_tab2,dim=1)
          allocate(sfh_tab2_trans(size_disc2,nfile))
          sfh_tab2_trans    = 0.0d0
          sfh_tab2_trans    = gal%disc%sfh_tab2 *  sta_transf2 / gal%disc%minstar2
          gal%disc%sfh_tab2 = gal%disc%sfh_tab2 - sfh_tab2_trans                
          gal%disc%minstar2 = sum(gal%disc%sfh_tab2)
          sta_transf2       = sum(sfh_tab2_trans)
          gal%disc%transp2  = gal%disc%transp2  + sta_transf2
       else
          size_disc2        = 1
          allocate(sfh_tab2_trans(size_disc2,nfile))
          sfh_tab2_trans    = 0.0d0
       end if
#endif
       ! remove gas from disc
       if (gas_transf > 0.0d0) then 
          met_transf      = gas_transf * gal%disc%mcoldz / gal%disc%mcold
          gal%disc%mcold  = max(gal%disc%mcold  - gas_transf,0.0d0)
          gal%disc%mcoldz = max(gal%disc%mcoldz - met_transf,0.0d0)
       else
          met_transf         = 0.0d0
       end if
       ! recompute bulge/burst properties (and add stars removed from the disc to the burst). 
#ifdef DUAL_IMF
       call recompute_bulge(gal,h,sta_transf,sfh_tab_trans,gas_transf,met_transf,size_disc,sta_transf2,sfh_tab2_trans,size_disc2)
       deallocate(sfh_tab2_trans)
#else
       call recompute_bulge(gal,h,sta_transf,sfh_tab_trans,gas_transf,met_transf,size_disc)
#endif
       deallocate(sfh_tab_trans)
    
#ifdef RECORD_SFR
#ifdef DUAL_IMF 
       call transfer_sfr_tabs(gal,sfr_transf,sfr_transf2)
#else 
       call transfer_sfr_tabs(gal,sfr_transf)
#endif
#endif

    endif
    
    macc = 0.0d0
    call disc_dynamics(gal,macc,h)

    return

  end subroutine check_for_disc_instability

!******************************************************************************************************************
#ifdef CLUMPY_SF  
  subroutine buxxx_dynamics(gal,h)

    implicit none

    type(galaxy) :: gal
    type(halo)   :: h
    real(kind=8) :: mtot_bulge,rhalf_bulge,core_mass

    mtot_bulge       = gal%burst%mgal + gal%bulge%mgal
    gal%bulge%rgal   = 0.1 * gal%disc%rgal * disc_params%s_to_m / bulge_params%s_to_m  
    !gal%bulge%rgal   = gal%disc%rgal * disc_params%s_to_m / bulge_params%s_to_m  

    gal%bulge%rgal   = max(gal%bulge%rgal,min_size)  

    rhalf_bulge      = gal%bulge%rgal * bulge_params%s_to_m
    call compute_core_mass(gal,rhalf_bulge,core_mass,'nope')
    core_mass        = gal%core%fdm*min(2.d0*core_mass,gal%core%mass)
    ! 3D vel disp @ half mass radius --> 0.177 factor comes also from Hernquist profile for bulges (ApJ 1990)
    gal%bulge%speed  = sqrt(0.177d0*gravconst*(mtot_bulge+core_mass)/gal%bulge%rgal)
    ! size of burst according to eq 4.5 of GalICS I
    gal%burst%rgal   = max(0.1 * gal%bulge%rgal, min_size)
    ! velocity of burst is the same as that of bulge
    gal%burst%speed  = gal%bulge%speed
    ! compute new dynamical timescales
    gal%bulge%tdyn   = bulge_tdyn(gal%bulge)       
    gal%burst%tdyn   = bulge_tdyn(gal%burst)

    return
    
  end subroutine buxxx_dynamics
#endif
!******************************************************************************************************************
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

end module BOOKKEEP_BARYONS




