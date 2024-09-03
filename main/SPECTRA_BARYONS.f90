module SPECTRA_BARYONS

  use IO_BARYONS
  use BOOKKEEP_BARYONS
  
  public
  
contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!*****************************************************************************************************************
  subroutine compute_spectra(st)

    ! compute and output galaxy spectra of all galaxies in halo h, at timestep st

    implicit none

    real(kind=8)              :: r_shell
    real(kind=8)              :: bulge_spec_aa(lmax),burst_spec_aa(lmax),disc_spec_aa(lmax),fdisc_spec_aa(lmax)   ! spectrum after dust absorption
    real(kind=8)              :: uv_1200_dblg_ext(21)                                                             ! total spectrum after dust absorption
    integer(kind=4)           :: st
    character(MAXPATHSIZE)    :: fname
    logical(kind=4)           :: compute_spectra_flag     ! set to true if we want to compute spectra.
    type(galaxy)              :: gal                      ! current galaxy to compute magnitudes for
    integer(kind=4)           :: j,k                      ! halo/galaxy indexes for main loop
    integer(kind=4)           :: i                        ! dummy loosp variable
    integer(kind=4)           :: nb_tot_gals              ! local copy of tsno(st)%n_total
    real(kind=8)              :: bol(ncomp)               ! bolometric luminosity of gal (in the ncomp components)
    real(kind=8)              :: IR_bol(ncomp)            ! IR bol lum of gal (in the ncomp components)
    real(kind=8)              :: spectrum(lmax,ncomp)     ! stellar spectrum of gal 
    character(80)             :: fmtstring 
    real(kind=8)              :: inc,xxx,yyy    
    real(kind=8)              :: shift, shift2            ! (1+z(ts)), (1+z(ts-1))
    real(kind=8)              :: dz                       ! z(ts-1) - z(ts)
#ifdef MOMAF_INPUTS
    integer(kind=4)           :: igal,ifilt !,icomp
    real(kind=8),allocatable  :: momafmag(:,:,:)
    real(kind=8),allocatable  :: momafdmag(:,:,:)
    real(kind=8),allocatable  :: momag(:),modmag(:)
    character(MAXPATHSIZE)    :: filename
#endif
    integer(kind=4)           :: do_output

#ifdef DESIRED_OUTPUTS
    integer(kind=4)          :: i_st
#endif

#ifdef LYA  
    real(kind=8)             :: logLlya_out,Llya_out
#endif

#ifndef MOMAF_INPUTS
    if (n_ascii_want == 0) then 
       return
    end if
#endif

    ! define format to write magnitudes 
    write(fmtstring,'(a1,i2,a)') '(',nftot,'(f7.3,1x))'

    ! decide if we output ascii files 
    rf_flag = .false.
    do i = 1, n_ascii_want
       if (st == ascii_list(i)) rf_flag = .true.
    end do

    ! initialise stuff for obs-frame mags with cosmic extinction
#ifdef TIBO
    shift              = 1.0d0/tsno(st)%aexp                   ! = (1+z) of timestep st  if simu stops before z=0
#else
    shift              = (tsno(nsteps)%aexp/tsno(st)%aexp)    ! = (1+z) of timestep st
#endif

    if (st > 1) then 
#ifdef TIBO
       shift2             = 1.0d0/tsno(st-1)%aexp                 ! = (1+z) of timestep st-1 if simu stops before z=0
#else
       shift2             = (tsno(nsteps)%aexp/tsno(st-1)%aexp)  ! = (1+z) of timestep st-1
#endif
    else
       shift2 = shift + 1.0d0   ! if first timestep, then just use a dz=1 for computing dm/dz ...
    end if
    dz                 = shift2 - shift                       ! = z(ts-1) - z(ts) (> 0)
    call get_cosmic_extinction(shift2)              ! compute cosmic extinction @ z(ts-1)
    tot_wave2          = restframe_wave *shift2      ! = wavelength @ z(ts-1) in observer frame
    cosmic_extinction2 = cosmic_extinction/shift2   ! do division here to save doing for each galaxy later
    call get_cosmic_extinction(shift)               ! same @ z(ts)
    tot_wave           = restframe_wave *shift
    call ce_effect                                  ! compute mag change in each filter due to cosmic extinction 
    write(ce_unit,fmtstring) amabs_diff(1:nftot)         ! write that out to file ce_effect.dat
    cosmic_extinction  = cosmic_extinction/shift    ! do division here to save doing in all_that_spectra_stuff

#ifdef MOMAF_INPUTS
    compute_spectra_flag = .true.
#else
    compute_spectra_flag = rf_flag
#endif

    if (compute_spectra_flag) then 

       nb_tot_gals        = tsno(st)%n_total                     ! total number of gals in timestep st

       if (rf_flag) then
          write(errunit,*) '> timestep ',st,': will output ascii file for rest-frame (rf) magnitudes'
          call open_magfile(46,'rf_mags_disc',st,nb_tot_gals)
          call open_magfile(47,'rf_mags_face',st,nb_tot_gals)
          call open_magfile(48,'rf_mags_bulge',st,nb_tot_gals)
          call open_magfile(49,'rf_mags_burst',st,nb_tot_gals)
          call open_magfile(96,'rf_noext_mags_disc',st,nb_tot_gals)
          call open_magfile(97,'rf_noext_mags_face',st,nb_tot_gals)
          call open_magfile(98,'rf_noext_mags_bulge',st,nb_tot_gals)
          call open_magfile(99,'rf_noext_mags_burst',st,nb_tot_gals)
       end if
       call open_magfile(70,'bol_lum',st,nb_tot_gals)       
#ifdef MOMAF_INPUTS
       if (nb_tot_gals > 0) then 
          allocate(momafmag(nb_tot_gals,3,nftot),momafdmag(nb_tot_gals,3,nftot))
          momafmag  = 99.0d0
          momafdmag = 0.0d0
          igal      = 1
       end if
#endif

       do_output = 0
       
#ifdef DESIRED_OUTPUTS
       do i_st = 1,ndes_ts
          if (st .ne. desired_ts(i_st)) then
             do_output = 1
          else
             do_output = 0
             go to 1
          end if
       end do
#endif
       
    1 continue
       
       if (do_output == 0) then ! global output, i.e output "all" info

          write(fname,'(a,a,i3.3)') trim(data_dir),'/gal_Nion_phot.',st
          open(444,file=fname,status='unknown')
          write(444,'(i3,1x,i6)') st,tsno(st)%n_total
          
          write(fname,'(a,a,i3.3)') trim(data_dir),'/uv1200_gal.',st
          open(uv1200_gal_unit,file=fname,status='unknown')
          write(uv1200_gal_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
          
          !tibo
#ifndef MIN_IO
          write(fname,'(a,a,i3.3)') trim(data_dir),'/uv1200_dblg_ext.',st
          open(uv1200_dblg_ext_unit,file=fname,status='unknown')
          write(uv1200_dblg_ext_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
          
          write(fname,'(a,a,i3.3)') trim(data_dir),'/uv1200_burst.',st
          open(uv1200_burst_unit,file=fname,status='unknown')
          write(uv1200_burst_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
          
          write(fname,'(a,a,i3.3)') trim(data_dir),'/uv1200_dblg.',st
          open(uv1200_dblg_unit,file=fname,status='unknown')
          write(uv1200_dblg_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
#endif
          
          write(fname,'(a,a,i3.3)') trim(data_dir),'/uv1500_gal.',st
          open(uv1500_gal_unit,file=fname,status='unknown')
          write(uv1500_gal_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ       
          
          write(fname,'(a,a,i3.3)') trim(data_dir),'/uvSED_beta_gal.',st
          open(beta_unit,file=fname,status='unknown')
          write(beta_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ  
          
          write(fname,'(a,a,i3.3)') trim(data_dir),'/age_lum_weighted.',st
          open(age_lum_weighted_unit,file=fname,status='unknown')
          write(age_lum_weighted_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ       
          
       endif

       ! TIBO
#ifdef LYA         
       if (st .gt. 4 .and. st .le. 9) then
          write(fname,'(a,a,i2.2)') trim(data_dir),'/llya4galics.',st
          open(unit=lyaout_unit,file=fname,status='unknown')
       endif
       if (st .gt. 9) then
          write(fname,'(a,a,i3.3)') trim(data_dir),'/llya4galics.',st
          open(unit=lyaout_unit,file=fname,status='unknown')
       endif
#endif       
      ! OBIT
       
      
       do j=1,tsno(st)%nb_of_halos
          do k=1, tsno(st)%liste_halos(j)%datas%nbgal
#ifdef LYA
             logLlya_out = -99.0d0
             Llya_out = 0.0d0
             if (st .gt. 4) then
                read(lyaout_unit,'(e14.6)') logLlya_out
                !Llya_out = 10.0d0**logLlya_out
                Llya_out = logLlya_out
             endif
#endif             
             ! initialize a few things
             gal                   = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)
             spectrum              = 0.0d0
             bol                   = - 1.0d0 - log_L_sun 
             IR_bol                = - 1.0d0 - log_L_sun 
             bulge_spec_aa(:)      = 0.0d0
             burst_spec_aa(:)      = 0.0d0
             disc_spec_aa(:)       = 0.0d0
             fdisc_spec_aa(:)      = 0.0d0
             uv_1200_dblg_ext(:)   = 0.0d0

             call calc_spectrum(gal,spectrum,st) ! compute spectrum of stars

             ! dont bother calculating spectrum if no stars 
!!$#ifdef DUAL_IMF
!!$             if (gal%disc%minstar > 0.0d0 .or. gal%bulge%minstar > 0.0d0 .or. gal%burst%minstar > 0.0d0 .or. &
!!$                  & gal%disc%minstar2 > 0.0d0 .or. gal%bulge%minstar2 > 0.0d0 .or. gal%burst%minstar2 > 0.0d0) then
!!$#else
!!$             if (gal%disc%minstar > 0.0d0 .or. gal%bulge%minstar > 0.0d0 .or. gal%burst%minstar > 0.0d0) then 
!!$#endif                                
                !call calc_spectrum(gal,spectrum,st) ! compute spectrum of stars
                
                ! disk spectrum
                !---------------
!!$#ifdef DUAL_IMF
!!$                if (gal%disc%minstar > 0.0d0 .or. gal%disc%minstar2 > 0.0d0) then 
!!$#else
!!$                if (gal%disc%minstar > 0.0d0) then 
!!$#endif

                   inc = gal%inclination
                   r_shell = gal%disc%rgal
#ifdef LYA                   
                   call all_that_spectra_stuff(r_shell,gal%disc,spectrum(:,1),bol(1),IR_bol(1),& 
                        & inc,'scrn','disc',disc_spec_aa,st,Llya_out)
#else
 call all_that_spectra_stuff(r_shell,gal%disc,spectrum(:,1),bol(1),IR_bol(1),& 
                        & inc,'scrn','disc',disc_spec_aa)
#endif                   
#ifdef MOMAF_INPUTS
                   momafmag(igal,1,:)  = amabs_ce
                   momafdmag(igal,1,:) = amabs_diff / dz
#endif
!!$                else
!!$                   amabs_rf = 99.0d0 ; amabs_rf_noext = 99.0d0
!!$                   bol(1) = - 1.0d0 - log_L_sun ; IR_bol(1) = - 1.0d0 - log_L_sun 
!!$                end if
                if (rf_flag) then                   
                   write(46,fmtstring) amabs_rf
                   write(96,fmtstring) amabs_rf_noext
                end if
                
                ! face-on disc spectrum
                !-----------------------
                if (rf_flag) then  ! dont want face-on stuff for momaf... 
!!$#ifdef DUAL_IMF
!!$                   if (gal%disc%minstar > 0.0d0 .or. gal%disc%minstar2 > 0.0d0) then         
!!$#else
!!$                   if (gal%disc%minstar > 0.0d0) then         
!!$#endif
                   inc = 1.0d0
#ifdef LYA                   
                   call all_that_spectra_stuff(r_shell,gal%disc,spectrum(:,1),xxx,yyy,& 
                        & inc,'scrn','disc',fdisc_spec_aa,st,Llya_out)
#else
                   call all_that_spectra_stuff(r_shell,gal%disc,spectrum(:,1),xxx,yyy,& 
                        & inc,'scrn','disc',fdisc_spec_aa)
#endif
  
!!$                   else
!!$                      amabs_rf = 99.0d0 ; amabs_rf_noext = 99.0d0
!!$                   end if
                   write(47,fmtstring) amabs_rf
                   write(97,fmtstring) amabs_rf_noext
                end if

                ! bulge spectrum
                inc = 0.0d0
!!$#ifdef DUAL_IMF
!!$                if (gal%bulge%minstar > 0.0d0 .or. gal%bulge%minstar2 > 0.0d0) then
!!$#else
!!$                if (gal%bulge%minstar > 0.0d0) then
!!$#endif

#ifdef LYA                
                   call all_that_spectra_stuff(r_shell,gal%bulge,spectrum(:,2),bol(2),IR_bol(2),& 
                        & inc,'scrn','bulge',bulge_spec_aa,st,Llya_out)
#else
                  call all_that_spectra_stuff(r_shell,gal%bulge,spectrum(:,2),bol(2),IR_bol(2),& 
                        & inc,'scrn','bulge',bulge_spec_aa)
#endif
                  
#ifdef MOMAF_INPUTS
                   momafmag(igal,2,:)  = amabs_ce
                   momafdmag(igal,2,:) = amabs_diff / dz
#endif
!!$                else
!!$                   amabs_rf = 99.0d0 ; amabs_rf_noext = 99.0d0
!!$                   bol(2) = - 1.0d0 - log_L_sun ; IR_bol(2) = - 1.0d0 - log_L_sun 
!!$                endif
                if (rf_flag) then
                   write(48,fmtstring) amabs_rf
                   write(98,fmtstring) amabs_rf_noext
                end if

                ! burst spectrum
                !----------------
!!$#ifdef DUAL_IMF
!!$                if (gal%burst%minstar > 0.0d0 .or. gal%burst%minstar2 > 0.0d0) then
!!$#else
!!$                if (gal%burst%minstar > 0.0d0) then
!!$#endif

#ifdef LYA                
                call all_that_spectra_stuff(r_shell,gal%burst,spectrum(:,3),bol(3),IR_bol(3),& 
                     & inc,'scrn','burst',burst_spec_aa,st,Llya_out)
#else
                call all_that_spectra_stuff(r_shell,gal%burst,spectrum(:,3),bol(3),IR_bol(3),& 
                     & inc,'scrn','burst',burst_spec_aa)
#endif
                
#ifdef MOMAF_INPUTS
                   momafmag(igal,3,:)  = amabs_ce
                   momafdmag(igal,3,:) = amabs_diff / dz
#endif
!!$                else
!!$                   amabs_rf = 99.0d0 ; amabs_rf_noext = 99.0d0
!!$                   bol(3) = - 1.0d0 - log_L_sun ; IR_bol(3) = - 1.0d0 - log_L_sun 
!!$                endif
                if (rf_flag) then
                   write(49,fmtstring) amabs_rf
                   write(99,fmtstring) amabs_rf_noext
                end if
                
!!$             else ! if there aint no stars
!!$                
!!$                amabs_rf = 99.0d0 ; amabs_rf_noext = 99.0d0
!!$                bol = - 1.0d0 - log_L_sun ; IR_bol = - 1.0d0 - log_L_sun 
!!$                if (rf_flag) then
!!$                   write(46,fmtstring) amabs_rf
!!$                   write(47,fmtstring) amabs_rf
!!$                   write(48,fmtstring) amabs_rf
!!$                   write(49,fmtstring) amabs_rf
!!$                   write(96,fmtstring) amabs_rf_noext
!!$                   write(97,fmtstring) amabs_rf_noext
!!$                   write(98,fmtstring) amabs_rf_noext
!!$                   write(99,fmtstring) amabs_rf_noext
!!$                end if
!!$             end if
   
                if (do_output == 0) then            
!tibo
#ifndef MIN_IO
                   uv_1200_dblg_ext(:) = disc_spec_aa(iluv_min:iluv_max) + bulge_spec_aa(iluv_min:iluv_max)
#endif
                   write(70,'(2(e14.6,1x))') sum(10.0d0**bol),sum(10.0d0**IR_bol)
                   
                   !tibo
#ifndef MIN_IO
                   write(uv1200_dblg_ext_unit,'(22(e14.6,1x))',advance='yes') sum(uv_1200_dblg_ext)/21.0d0, uv_1200_dblg_ext(:)
#endif
                   
#ifdef MOMAF_INPUTS
                   igal = igal + 1
#endif
                end if
             end do
          end do

          ! TIBO
#ifdef LYA          
          if (st .gt. 4) then
             close(lyaout_unit)
          endif
#endif          
          ! OBIT
          
          if (do_output == 0) then   
             close(444)
             close(uv1200_gal_unit)
             !tibo
#ifndef MIN_IO
             close(uv1200_burst_unit)
             close(uv1200_dblg_unit)
             close(uv1200_dblg_ext_unit)
#endif
             close(uv1500_gal_unit)
             close(age_lum_weighted_unit)
             close(beta_unit)
          end if
       
       ! close output files of the timestep
       do i=46,49
          if (rf_flag) then
             close(i)
             close(i + 50)
          end if
          close(i+10) ; close(i+20)
       end do
       close(70)
      
#ifdef MOMAF_INPUTS
       ! write out momaf inputs at once
       if (nb_tot_gals > 0) then 
          allocate(momag(nb_tot_gals),modmag(nb_tot_gals))
          do ifilt = 1,nftot
             ! disc
             momag  = momafmag(:,1,ifilt)
             modmag = momafdmag(:,1,ifilt)
             write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/disc_',trim(filt(ifilt)%name),'.',st
             filename = trim(filename)//char(0) ! needed for call to C
             call write_momaf_magfile(filename,real(momag,4),real(modmag,4),nb_tot_gals)
             ! bulge
             momag  = momafmag(:,2,ifilt)
             modmag = momafdmag(:,2,ifilt)
             write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/bulge_',trim(filt(ifilt)%name),'.',st
             filename = trim(filename)//char(0)
             call write_momaf_magfile(filename,real(momag,4),real(modmag,4),nb_tot_gals)
             ! burst
             momag  = momafmag(:,3,ifilt)
             modmag = momafdmag(:,3,ifilt)
             write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/burst_',trim(filt(ifilt)%name),'.',st
             filename = trim(filename)//char(0)
             call write_momaf_magfile(filename,real(momag,4),real(modmag,4),nb_tot_gals)
          end do
          deallocate(momafmag,momafdmag,momag,modmag)
       else
          do ifilt = 1, nftot
             write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/disc_',trim(filt(ifilt)%name),'.',st
             filename = trim(filename)//char(0) ! "char(0)" needed for call to C
             call write_momaf_empty_magfile(filename)
             write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/bulge_',trim(filt(ifilt)%name),'.',st
             filename = trim(filename)//char(0) ! needed for call to C
             call write_momaf_empty_magfile(filename)
             write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/burst_',trim(filt(ifilt)%name),'.',st
             filename = trim(filename)//char(0) ! needed for call to C
             call write_momaf_empty_magfile(filename)
          end do
       end if
#endif
 
    end if

    return

  contains

    !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    subroutine get_cosmic_extinction(zplus1)

      implicit none

      real(kind=8)    :: zplus1,z
      integer(kind=4) :: i

      z = zplus1 -1.0d0  !redshift of the source
      do i = 1, lmax+fins
         cosmic_extinction(i) = igm_transmission(z, zplus1 * restframe_wave(i)*10000.0d0)  ! zplus1 * rest = total. 
      end do

      return

    end subroutine get_cosmic_extinction

    !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  end subroutine compute_spectra

!*****************************************************************************************************************
  real function igm_transmission(z_source,lambda_obs)
    
    implicit none
    
    real(kind=8)     :: z_source,lambda_obs
    real(kind=8)     :: A(1:4)
    real(kind=8)     :: lambda(1:5), red_lambda(1:5)
    real(kind=8)     :: xc, xem, igm_absorption
    
    ! define constants from Madau et al 1996 :
    A(1)             = 0.0036d0    ! Ly_alpha coefficient
    A(2)             = 1.7e-3    ! Ly_beta ...
    A(3)             = 1.2e-3    ! Ly_gamma ...
    A(4)             = 9.3e-4    ! Ly_delta ...
    
    lambda(1)        = 1216.d0 ! Ly_alpha in angstroms
    lambda(2)        = 1026.d0 ! Ly_beta
    lambda(3)        = 973.d0  ! Ly_gamma
    lambda(4)        = 950.d0  ! Ly_delta
    lambda(5)        = 912.d0  ! Lyman limit
    
    red_lambda(1:5)  = (1.d0+z_source) * lambda  ! redshifted lambdas
    
    ! variables for continuum absorption :
    xc               = lambda_obs/lambda(5)
    xem              = 1.0d0 + z_source
    
    ! We consider that between 912A and 950A, there is only Ly_delta
    ! absorption, because lacking coefficients in Madau et al 1996.
    
    ! initialisation : if no absorption, igm_abs = 0.0 :
    igm_absorption   = 0.0d0
    
    ! Ly_alpha only :
    if ((lambda_obs <= red_lambda(1)) .and. (lambda_obs > red_lambda(2))) then
       igm_absorption = A(1) * (lambda_obs/lambda(1))**(3.46d0)
    end if
    
    ! Ly_alpha and Ly_beta :
    if ((lambda_obs <= red_lambda(2)) .and. (lambda_obs > red_lambda(3))) then
       igm_absorption = A(1) * (lambda_obs/lambda(1))**(3.46d0) + &
                      & A(2) * (lambda_obs/lambda(2))**(3.46d0)
    end if
    
    !Ly_alpha, Ly_beta, and Ly_gamma :
    if ((lambda_obs <= red_lambda(3)) .and. (lambda_obs > red_lambda(4))) then
       igm_absorption = A(1) * (lambda_obs/lambda(1))**(3.46d0) + &
                      & A(2) * (lambda_obs/lambda(2))**(3.46d0) + &
                      & A(3) * (lambda_obs/lambda(3))**(3.46d0)
    end if
    
    !Ly_alpha, Ly_beta, Ly_gamma, and Ly_delta :
    if ((lambda_obs <= red_lambda(4)) .and. (lambda_obs > red_lambda(5))) then
       igm_absorption = A(1) * (lambda_obs/lambda(1))**(3.46d0) + &
                      & A(2) * (lambda_obs/lambda(2))**(3.46d0) + &
                      & A(3) * (lambda_obs/lambda(3))**(3.46d0) + &
                      & A(4) * (lambda_obs/lambda(4))**(3.46d0)
    end if
    
    ! All these lines plus the photoelectric effect (shortwards 912 A)
    if (lambda_obs <= red_lambda(5)) then 
       igm_absorption = A(1) * (lambda_obs/lambda(1))**(3.46d0)      + &
                      & A(2) * (lambda_obs/lambda(2))**(3.46d0)      + &
                      & A(3) * (lambda_obs/lambda(3))**(3.46d0)      + &
                      & A(4) * (lambda_obs/lambda(4))**(3.46d0)      + &
                      & 0.25 * xc**3 * (xem**0.46d0 - xc**0.46d0)      + &
                      & 9.4 * xc**1.5 * (xem**0.18d0 - xc**0.18d0)     - &
                      & 0.7 * xc**3 * (xc**(-1.32d0) - xem**(-1.32d0)) - &
                      & 0.023 * (xem**1.68d0 - xc**1.68d0)
    end if
    
    igm_transmission = exp(-igm_absorption)
    
    return
    
  end function igm_transmission

!*****************************************************************************************************************
  subroutine ce_effect

  ! we compute the difference in magnitude due to cosmic extinction in each filter, for a flat spectrum.

    implicit none         
    
    real(kind=8) :: amabs(nftot)
    real(kind=8) :: copied_spectrum(lmax+fins)  

    copied_spectrum(1:lmax+fins) = 1.0d0                                    ! spectrum before CE
    call convolution(copied_spectrum,tot_wave,lmax+fins,amabs)      ! mags before CE 
    amabs_diff                   = amabs
    copied_spectrum(1:lmax+fins) = cosmic_extinction(1:lmax+fins)         ! spectrum after CE
    call convolution(copied_spectrum,tot_wave,lmax+fins,amabs)      ! mags after CE
    amabs_diff                   = amabs_diff - amabs                     ! difference as global variable amabs_diff.  

    return

  end subroutine ce_effect

!*****************************************************************************************************************
  subroutine convolution(spectrum,lambda,n,flu)

  ! computation of fluxes through the filters to get observed magnitudes and fluxes

    implicit none

    integer(kind=4),intent(in) :: n
    real(kind=8),intent(in)    :: spectrum(n),lambda(n)
    real(kind=8),intent(out)   :: flu(nftot)
    integer(kind=4)            :: l,nf,lf,i_up,i_down    
    real(kind=8)               :: r1
    real(kind=8),allocatable   :: filter_wave(:),filter_trans(:)
    integer(kind=4)            :: filter_lftot

    flu = 0.0d0 

    do nf = 1,nftot

       filter_lftot = filt(nf)%lftot
       allocate(filter_wave(filter_lftot),filter_trans(filter_lftot))
       filter_wave  = filt(nf)%wave
       filter_trans = filt(nf)%trans

       call locate(lambda,n,filter_wave(1),i_down)            ! lowest wavelength of filter 
       i_up = i_down
       call hunt(lambda,n,filter_wave(filter_lftot),i_up)     ! highest wavelength of filter 
       lf = 1                                                 ! good first guess....
       i_down = max(i_down,1) ;  i_up = min(i_up+1,n-1)       ! just for safety's sake........ 
       if (i_up-i_down >= filter_lftot) then                  ! spectrum has better resolution than filter
          do l=i_down,i_up
             call hunt(filter_wave(1:filter_lftot),filter_lftot,lambda(l),lf)
             lf = max(lf,1)
             if (lf < filter_lftot) then 
                r1      = filter_trans(lf)+(filter_trans(lf+1) - filter_trans(lf))* &
                     & (lambda(l)-filter_wave(lf)) / (filter_wave(lf+1) - filter_wave(lf))                
                flu(nf) = flu(nf) + (spectrum(l)+spectrum(l+1)) * r1 * (lambda(l+1)-lambda(l))/2.d0
             end if
          end do
       else                                                   ! filter has better resolution than spectrum 
          do lf=1,filter_lftot -1 
             call hunt(lambda,n,filter_wave(lf),l)
             l = max(l,1)
             if (l < n) then                 
                r1      = spectrum(l) + (spectrum(l+1) - spectrum(l))*(filter_wave(lf) - lambda(l)) &
                     & / (lambda(l+1) - lambda(l)) 
                flu(nf) = flu(nf) + (filter_trans(lf) + filter_trans(lf+1)) * r1 * & 
                     & (filter_wave(lf+1) - filter_wave(lf))/2.d0
             end if
          end do
       end if
       
       deallocate(filter_wave,filter_trans)

    end do

    do nf=1,nftot
       if (flu(nf) > 0.0d0) then
          flu(nf) = filt(nf)%cfil - 2.5d0*log10(flu(nf)) 
       else
          flu(nf) = 99.0d0
       endif
    enddo
                      
    return

  end subroutine convolution

!*****************************************************************************************************************
  subroutine calc_spectrum(g,spectrum,st)  
  
    implicit none    

    integer(kind=4)                    :: st
    character(MAXPATHSIZE)             :: fname
    type(galaxy)                       :: g
    real(kind=8),dimension(lmax,ncomp) :: spectrum
    integer(kind=4)                    :: i,j,nagebins,num_mets
    character(MAXPATHSIZE)             :: file
    real(kind=8),allocatable           :: Nphot_tab(:,:)
    real(kind=8)                       :: age_lc_disc,age_lc_burst,age_lc_bulge,lc_disc,lc_burst,lc_bulge            
    real(kind=8),dimension(21)         :: uv1500_gal
    real(kind=8),dimension(67)         :: uvSED_beta
!tibo
    integer(kind=4)                    :: do_output
#ifdef DESIRED_OUTPUTS
    integer(kind=4)                    :: i_st
#endif


#ifdef TH4
    write(file,'(4a)') trim(inputpath),'/Nion_photons_tab_',trim(imfname),'.bin'  ! modif_imf : TH4 is contains a digit...
#else
    write(file,'(4a)') trim(inputpath),'/Nion_photons_tab_',imfname,'.bin'
#endif

    open(333,file=file,form='unformatted')
    rewind(333)
    read(333) nagebins,num_mets
    allocate(Nphot_tab(nagebins,num_mets))
    read(333) Nphot_tab(:,:)   
    close(333)

    
    g%disc%Nion_phot = 0.0d0
    if (g%disc%minstar > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%disc%sfh_tab,dim=1)
             spectrum(:,1)    = spectrum(:,1) + sburst(:,j,i) * g%disc%sfh_tab(j,i)
             g%disc%Nion_phot = g%disc%Nion_phot + 10.0d0**Nphot_tab(j,i) * g%disc%sfh_tab(j,i) * 10.0d0**11
          end do
       end do
    else
       spectrum(:,1)    = 0.0d0
    end if

!--------------------- mean age of the galaxy ---------------------------------

    age_lc_disc = 0.0d0
    lc_disc     = 0.0d0
    if (g%disc%minstar > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%disc%sfh_tab,dim=1)
             !age_lc_disc = age_lc_disc + g%disc%sfh_tab(j,i) * 10.0d0**11 * timetab(j) * sum(sburst(iluv_min:iluv_max,j,i))/21.0d0
             !lc_disc     = lc_disc + g%disc%sfh_tab(j,i) * 10.0d0**11 * sum(sburst(iluv_min:iluv_max,j,i))/21.0d0
             age_lc_disc = age_lc_disc + g%disc%sfh_tab(j,i) * 10.0d0**11 * timetab(j) * g%disc%minstar
             lc_disc     = lc_disc + g%disc%sfh_tab(j,i) * 10.0d0**11 * g%disc%minstar
          end do
       end do
    endif

    age_lc_burst = 0.0d0
    lc_burst     = 0.0d0
    if (g%burst%minstar > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%burst%sfh_tab,dim=1)
             !age_lc_burst = age_lc_burst + g%burst%sfh_tab(j,i) * 10.0d0**11 * timetab(j) * sum(sburst(iluv_min:iluv_max,j,i))/21.0d0
             !lc_burst     = lc_burst + g%burst%sfh_tab(j,i) * 10.0d0**11 * sum(sburst(iluv_min:iluv_max,j,i))/21.0d0
             age_lc_burst = age_lc_burst + g%burst%sfh_tab(j,i) * 10.0d0**11 * timetab(j) * g%burst%minstar
             lc_burst     = lc_burst + g%burst%sfh_tab(j,i) * 10.0d0**11 * g%burst%minstar
          end do
       end do
    endif

    age_lc_bulge = 0.0d0
    lc_bulge     = 0.0d0
    if (g%bulge%minstar > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%bulge%sfh_tab,dim=1)
             !age_lc_bulge = age_lc_bulge + g%bulge%sfh_tab(j,i) * 10.0d0**11 * timetab(j) * sum(sburst(iluv_min:iluv_max,j,i))/21.0d0
             !lc_bulge     = lc_bulge + g%bulge%sfh_tab(j,i) * 10.0d0**11 * sum(sburst(iluv_min:iluv_max,j,i))/21.0d0
             age_lc_bulge = age_lc_bulge + g%bulge%sfh_tab(j,i) * 10.0d0**11 * timetab(j) * g%bulge%minstar
             lc_bulge     = lc_bulge + g%bulge%sfh_tab(j,i) * 10.0d0**11 * g%bulge%minstar
          end do
       end do
    endif

             
    if (g%bulge%minstar .ne. 0.0d0 .or. g%burst%minstar .ne. 0.0d0 .or. g%disc%minstar .ne. 0.0d0) then
       g%age_lum_weighted = (age_lc_disc + age_lc_bulge + age_lc_burst) / (lc_disc + lc_burst + lc_bulge)
    else
       g%age_lum_weighted = -99.0d0
    end if

!------------------------------------------------------------------------------

    g%bulge%Nion_phot  = 0.0d0
    if (g%bulge%minstar > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%bulge%sfh_tab,dim=1)
             spectrum(:,2)     = spectrum(:,2) + sburst(:,j,i) * g%bulge%sfh_tab(j,i)
             g%bulge%Nion_phot = g%bulge%Nion_phot + 10.0d0**Nphot_tab(j,i) * g%bulge%sfh_tab(j,i) * 10.0d0**11
          end do
       end do
    else
       spectrum(:,2)      = 0.0d0
    end if

    g%burst%Nion_phot  = 0.0d0
    if (g%burst%minstar > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%burst%sfh_tab,dim=1)
             spectrum(:,3)     = spectrum(:,3) + sburst(:,j,i) * g%burst%sfh_tab(j,i)
             g%burst%Nion_phot = g%burst%Nion_phot + 10.0d0**Nphot_tab(j,i) * g%burst%sfh_tab(j,i) * 10.0d0**11
          end do
       end do
    else
       spectrum(:,3)      = 0.0d0
    end if

    
    g%uv_1200_gal(:) = spectrum(iluv_min:iluv_max,1) + spectrum(iluv_min:iluv_max,2) + spectrum(iluv_min:iluv_max,3)
    g%uv_1200_burst(:) = spectrum(iluv_min:iluv_max,3)
    g%uv_1200_dblg(:)  = spectrum(iluv_min:iluv_max,1) + spectrum(iluv_min:iluv_max,2) 

   uv1500_gal(:) = spectrum(iburstuv_min:iburstuv_max,1) + spectrum(iburstuv_min:iburstuv_max,2) + spectrum(iburstuv_min:iburstuv_max,3)

   uvSED_beta(:) = spectrum(ibeta_min:ibeta_max,1) + spectrum(ibeta_min:ibeta_max,2) + spectrum(ibeta_min:ibeta_max,3)

#ifdef DUAL_IMF
! add IMF2 stell-pop if present
    if (g%disc%minstar2 > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%disc%sfh_tab2,dim=1)
             spectrum(:,1) = spectrum(:,1) + sburst2(:,j,i) * g%disc%sfh_tab2(j,i)
             g%disc%Nion_phot = g%disc%Nion_phot + 10.0d0**Nphot_tab(j,i) * g%disc%sfh_tab2(j,i) * 10.0d0**11
          end do
       end do
    else
       spectrum(:,1)    = 0.0d0
       g%disc%Nion_phot = 0.0d0
    end if
    
    if (g%bulge%minstar2 > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%bulge%sfh_tab2,dim=1)
             spectrum(:,2) = spectrum(:,2) + sburst2(:,j,i) * g%bulge%sfh_tab2(j,i)
             g%bulge%Nion_phot = g%bulge%Nion_phot + 10.0d0**Nphot_tab(j,i) * g%bulge%sfh_tab2(j,i) * 10.0d0**11
          end do
       end do
    else
       spectrum(:,2)      = 0.0d0
       g%bulge%Nion_phot  = 0.0d0
    end if
    if (g%burst%minstar2 > 0.0d0) then
       do i = 1, nfile
          do j = 1, size(g%burst%sfh_tab2,dim=1)
             spectrum(:,3) = spectrum(:,3) + sburst2(:,j,i) * g%burst%sfh_tab2(j,i)
             g%burst%Nion_phot = g%burst%Nion_phot + 10.0d0**Nphot_tab(j,i) * g%burst%sfh_tab2(j,i) * 10.0d0**11
          end do
       end do
    else
       spectrum(:,3)      = 0.0d0
       g%burst%Nion_phot  = 0.0d0
    end if
    
#endif

! do a modif here for desired output...

    do_output = 0

#ifdef DESIRED_OUTPUTS
       do i_st = 1,ndes_ts
          if (st .ne. desired_ts(i_st)) then
             do_output = 1
          else
             do_output = 0
             go to 2
          end if
       end do
#endif

    2 continue

       if (do_output == 0) then ! global output, i.e output "all" info   

          write(444,'(e14.6,1x,e14.6,1x,e14.6)',advance='yes') g%disc%Nion_phot,g%burst%Nion_phot,g%bulge%Nion_phot    
          
          write(uv1200_gal_unit,'(22(e14.6,1x))',advance='yes') sum(g%uv_1200_gal(:)) / 21.0d0, g%uv_1200_gal(:)
          !tibo
#ifndef MIN_IO
          write(uv1200_burst_unit,'(22(e14.6,1x))',advance='yes') sum(g%uv_1200_burst(:)) / 21.0d0, g%uv_1200_burst(:)
          write(uv1200_dblg_unit,'(22(e14.6,1x))',advance='yes') sum(g%uv_1200_dblg(:)) / 21.0d0, g%uv_1200_dblg(:)
#endif
          
          write(uv1500_gal_unit,'(21(e14.6,1x))',advance='yes') uv1500_gal(:)
          write(beta_unit,'(67(e14.6,1x))',advance='yes') uvSED_beta(:)

          write(age_lum_weighted_unit,'(e14.6)',advance='yes') g%age_lum_weighted
  
       end if

#ifdef NO_ION       
       !! Tibo: set ionizing flux to 0 to be consistent with Lya emission
       spectrum(1:121,1) = 0.0d0
       spectrum(1:121,2) = 0.0d0
       spectrum(1:121,3) = 0.0d0
       
       ! Next bin is 915 A => interpolation
       spectrum(122,1) = spectrum(122,1) * (915.d0 - 912.d0) / (915.d0 - 905.d0)
       spectrum(122,2) = spectrum(122,2) * (915.d0 - 912.d0) / (915.d0 - 905.d0)
       spectrum(122,3) = spectrum(122,3) * (915.d0 - 912.d0) / (915.d0 - 905.d0)
#endif
       
    return
    
  end subroutine calc_spectrum

  !*****************************************************************************************************************
#ifdef LYA  
  subroutine all_that_spectra_stuff(r_shell,comp,spec,bol,IR_bol,inc,name,comp_type,new_spec,st,Llya_out)
#else
  subroutine all_that_spectra_stuff(r_shell,comp,spec,bol,IR_bol,inc,name,comp_type,new_spec)
#endif
    

    implicit none


    type(gal_comp)    :: comp
    real(kind=8)      :: spec(lmax),inc                    ! stellar spectrum and cos(inclination)
    character(4)      :: comp_type                         ! either disc, bulg, or burs
    character*4       :: name                              ! name of dust model ('Dwek', 'Slab')
    real(kind=8)      :: bol,IR_bol                        ! total and IR bolometric luminosity. 
    real(kind=8)      :: local_wave(lmax+fins)
    real(kind=8)      :: amabs(nftot)
    real(kind=8)      :: new_spec(lmax)                    ! stellar spectrum after dust absorption 
    real(kind=8)      :: final_spectrum(lmax+fins)         ! full spectrum, including dust emission. 
#ifdef MOMAF_INPUTS
    real(kind=8)      :: copied_spectrum(lmax+fins) 
#endif
    real(kind=8)      :: obsc_rad,r_shell

    ! TIBO
#ifdef LYA      
    real(kind=8)      :: Llya_out,Llya_out_2add
    integer(kind=4)   :: st
#endif    
    ! OBIT
    
    if (comp_type == 'disc') then
       obsc_rad = disc_params%obsc_rad
    else
       obsc_rad = bulge_params%obsc_rad
    end if
    
    if (rf_flag) then 
       new_spec       = spec 
       bol            = 0.0d0  
       IR_bol         = 0.0d0
       final_spectrum = 0.0d0
       call build_spectrum(final_spectrum,new_spec,IR_bol)
       final_spectrum = (final_spectrum/ 5.321d0 * 3.83d0)
       local_wave     = restframe_wave
       call convolution(final_spectrum,local_wave,lmax+fins,amabs)
       amabs_rf_noext = amabs
    end if

    new_spec        = spec       ! copy to internal variable in case we need this later
    call absorption(r_shell,comp,new_spec,IR_bol,bol,obsc_rad,name,inc) 
    final_spectrum  = 0.0d0
    ! build spectrum including dust emission 
    call build_spectrum(final_spectrum,new_spec,IR_bol)


    ! TIBO
!!$#ifdef LYA     
!!$    if (st .gt. 4) then
!!$       if (comp_type == 'disc') then
!!$          Llya_out_2add = -77.0d0
!!$          Llya_out_2add = Llya_out
!!$       end if
!!$                      
!!$       if (Llya_out_2add .gt. 0.0d0) then
!!$          final_spectrum(153) = final_spectrum(153) + Llya_out_2add
!!$       endif
!!$    end if
!!$#endif

#ifdef LYA     
    if (st .gt. 4) then
       if (comp_type == 'disc') then
          if (Llya_out .gt. 0.0d0) then
             Llya_out_2add = 10.0d0**Llya_out
             final_spectrum(153) = final_spectrum(153) + Llya_out_2add

      ! print*,Llya_out_2add,Llya_out,st
       ! These conversions have been done when writing the file : write_llya4galics.pro                   
       !logLlya_out_2add = logLlya_out - 1.0d0 ! erg/s -> erg/s/A (dlambda = 10 A)
       !logLlya_out_2add = logLlya_out_2add - log_L_sun - 11.0d0 + 4.0d0                                                                  
       !! Add Llya_out luminosity density to SED at Lambda = 1215-1225 AA    
      
          !print*,Llya_out_2add,Llya_out,final_spectrum(153),final_spectrum(153) + Llya_out_2add
          !final_spectrum(153) = final_spectrum(153) + Llya_out_2add
          !print*,Llya_out_2add,Llya_out,final_spectrum(153)
          endif
       end if
    end if
#endif
    
    ! OBIT

    
    
    final_spectrum  = (final_spectrum / 5.321d0 * 3.83d0)     ! 5.321 is L_sun in V band, 3.83 is L_bol_sun

    if (rf_flag) then
       local_wave = restframe_wave
       call convolution(final_spectrum,local_wave,lmax+fins,amabs)
       amabs_rf   = amabs
    end if

#ifdef MOMAF_INPUTS
    copied_spectrum = final_spectrum * cosmic_extinction
    local_wave      = tot_wave
    call convolution(copied_spectrum,local_wave,lmax+fins,amabs)
    amabs_ce        = amabs
    copied_spectrum = final_spectrum * cosmic_extinction2
    local_wave      = tot_wave2
    call convolution(copied_spectrum,local_wave,lmax+fins,amabs)
    amabs_diff      = amabs - amabs_ce
#endif

    return 

  end subroutine all_that_spectra_stuff

!*****************************************************************************************************************
  subroutine absorption(r_shell,gal,spectrum,alir,bol_lum,nh_par,model,cth)

! name is rather explicit: computes how much the stellar population synthesis spectra is affected by dust
    
    implicit none  
          
    type(gal_comp),intent(in)          :: gal
    real(kind=8),intent(out)           :: alir,bol_lum   ! bol_lum will be total bolometric lum, alir = IR lum.  
    real(kind=8),intent(in)            :: cth,nh_par     ! cth is cos(inclination ang), nh_par is the obscuration radius
    character(4),intent(in)            :: model
    real(kind=8),intent(inout)         :: spectrum(lmax)
    real(kind=8)                       :: Z_g,nh,al,nh_tibo,r_shell
    integer(kind=4)                    :: i,j
    real(kind=8),dimension(lmax)       :: extinction,ext_spe,ext_spe_int,tau,tau_tibo
    real(kind=4),parameter             :: mu_gas = 1.22

    if (gal%mcold <= 0.0d0) then  ! no absorption bc no cold gas

       if (abs(gal%mcoldz) > rel_prec) stop 'Pb with %coldz in absorption'
       ! Now, integrate the SED to get the total bolometric luminosity  
       bol_lum = 0.0d0
       alir    = 0.0d0
       do i=2,lmax
          bol_lum = bol_lum + (spectrum(i) + spectrum(i-1)) * dalamb(i)
       enddo
       alir      = log10(max(alir,0.1d0))    ! in log_L_sun: 0.1 is minimum value (anything but 0)
       bol_lum   = log10(max(bol_lum,0.1d0)) ! idem

    else

       if (gal%mcoldz <= 0.0d0) then  ! no absorption bc no metals in cold gas

          ! Now, integrate the SED to get the total bolometric luminosity
          bol_lum = 0.0d0
          alir    = 0.0d0
          do i=2,lmax
             bol_lum = bol_lum + (spectrum(i)+spectrum(i-1))*dalamb(i)
          enddo
          alir      = log10(max(alir,0.1d0))              ! in log_L_sun
          bol_lum   = log10(max(bol_lum,0.1d0))           ! idem

       
       else                         ! absorption (normal case)

          nh = gal%mcold / (pi * (gal%rgal*nh_par)**2 ) ! compute the HI column density.  
          nh = log10(nh) + 18.95d0                      ! this is 10^11/1.4 M_sun / Mpc^2 --> 
                                                        ! convert to atoms/cm^2 
          ! the commented bits in following equations are bc HI absorption has not been tested  
          Z_g = log10((gal%mcoldz / gal%mcold)/Z_sun)
!begin jeje 
          extinction(:) = 10**(tau_hires(:,1) + nh + tau_hires(:,2)*Z_g) !+ 10**(tau_hires(:,3) +nh) 
          extinction    = extinction * (1.0d0 + global_redshift)**(-0.5)
!end jeje


          ! two different models, Dwek and Varosi 1996 for spheroid, Slab for disc.  
          if (model == 'Dwek') then

             tau = extinction 
             j   = 10           ! first guess....
             do i=1,lmax
                if (tau(i) >= tau_tab(ntau)) then 
                   ext_spe(i) = 0.75d0/tau(i)
                else if (tau(i) < tau_tab(1)) then 
                   ext_spe(i) = 1.0d0 + tau(i)*series(1)
                else
                   call hunt(tau_tab,ntau,tau(i),j)
                   ext_spe(i) = xinterp(tau_tab(j),ext_tab(j),tau_tab(j+1),ext_tab(j+1),tau(i))
                end if
             end do
             ext_spe     = (ext_spe/(1.0d0 - albedo + albedo * ext_spe)) ! energy which is put in the IR
             ext_spe_int = 1.0d0 - ext_spe                               ! energy taken from unextincted spectrum

          else if (model == 'scrn') then
             if (r_shell == 0.0d0) then
                r_shell = gal%rgal
             end if
             nh_tibo     = gal%mcold / 4.0d0 / pi / mu_gas / r_shell**2 ! mu_gas = 1.22... ok
             nh_tibo     = log10(nh_tibo) + 19.15d0 ! 18.95  
             tau_tibo(:) = 10**(tau_hires(:,1) + nh_tibo + tau_hires(:,2)*Z_g)
             tau_tibo    = tau_tibo * (1.0d0 + global_redshift)**(-0.5)
             ext_spe     = exp(-tau_tibo)
             ext_spe_int = 1.0d0 - ext_spe  
          else if (model == 'Slab') then
             
             tau = extinction*sqrt(1.0d0-albedo)
             j   = 10           ! first guess....         
             do i=1,lmax
                al = tau(i)/cth
                if (al >= tau_tab(ntau)) then 
                   ext_spe(i) = 1.0d0 / al
                else if (al < tau_tab(1)) then 
                   ext_spe(i) = 1.0d0 - 0.5d0*al
                else
                   call hunt(tau_tab,ntau,al,j)
                   ext_spe(i) = xinterp(tau_tab(j),slab_tab(j),tau_tab(j+1),slab_tab(j+1),al)
                end if
             end do
             ext_spe_int = 1.0d0 - ext_spe           

             ! integrated spectrum:
             do i=1,lmax
                if (tau(i) >= tau_tab(ntau)) then 
                   ext_spe_int(i) = 1.0d0 - (0.5d0/tau(i))
                else if (tau(i) < tau_tab(1)) then 
                   ext_spe_int(i) = 0.5d0*tau(i)*(1.5d0 - log(tau(i)) +0.577215d0)
                else
                   call hunt(tau_tab,ntau,tau(i),j)
                   ext_spe_int(i) = xinterp(tau_tab(j),slab_int_tab(j),tau_tab(j+1),slab_int_tab(j+1),tau(i))
                end if
             end do

          else if (model == 'QSO ') then

             tau = extinction
             j   = 10           ! first guess....         
!             do i=1,lmax  
             ! apply some model more ext_spe, the single, densest line of sight.  
!             end do
             ext_spe_int = 1.0d0 - ext_spe           
       
             ! integrated spectrum:
!             do i=1,lmax
             ! apply some model more ext_spe_int, the total extinction over all l.o.s
!             end do
       
          else

             stop '> error: model not implemented in absorption!'

          end if
    
          ! Now, integrate the SED to get the total bolometric luminosity, and the amount absorbed.  
          bol_lum = 0.0d0
          alir    = 0.0d0
    
          do i=2,lmax
             alir    = alir + (ext_spe_int(i)*spectrum(i)+ext_spe_int(i-1)*spectrum(i-1))*dalamb(i)
             bol_lum = bol_lum + (spectrum(i)+spectrum(i-1))*dalamb(i)
          enddo

          spectrum  = spectrum * ext_spe      ! new spectrum after dust extinction 
          alir      = log10(max(alir,0.1d0))    ! in log_L_sun
          bol_lum   = log10(max(bol_lum,0.1d0)) ! idem.

       endif

    endif
      
    return

  end subroutine absorption

!*****************************************************************************************************************
  subroutine build_spectrum(final_spectrum,anew_spec,alir)

! once again pretty explicit name: this sub takes the stellar population spectrum after it has been 
! corrected by dust absorption and patches it with the dust emission spectrum (conserving energy)
    
    implicit none

    real(kind=8),intent(in)    :: anew_spec(lmax),alir    
    real(kind=8),intent(out)   :: final_spectrum(fins+lmax)
    integer(kind=4)            :: i,j,k,l,jj
    real(kind=8)               :: red(fins)

    ! compute dust emission spectrum corresponding to correct bolometric luminosity via interpolation
    ! of table containing template spectra    
    if (alir <= lum(1)) then
       red(:) = alir -  (lum(1)) + (spec_ir(:,1))
    else if (alir >= lum(nbre)) then
       red(:) = alir -  (lum(nbre)) + (spec_ir(:,nbre))
    else
       call locate (lum,nbre,alir,jj)
       do i=1,fins
          red(i) = xinterp(lum(jj),spec_ir(i,jj),lum(jj+1),spec_ir(i,jj+1),alir)
       enddo
    end if 
   
    ! now patch it to the absorbed stellar population synthesis spectrum
    l=1 ;  j=1  ;  k=1
    do while (k <= fins+lmax)
       if (l <= lmax) then  
          if (alamb(l) < alam_ir(j)) then 
             if (j == 1) then
                final_spectrum(k) = anew_spec(l)
             else 
                final_spectrum(k) = anew_spec(l) + 10**(red(j-1) + interp_arr(k)*(red(j)-red(j-1)))
             endif
             l = l+1
          else
             if (alam_ir(j) < alamb(lmax)) then
                if (anew_spec(l-1) > 0.0d0) then 
                   final_spectrum(k) = 10**red(j) + anew_spec(l-1)*(anew_spec(l)/anew_spec(l-1))**interp_arr(k)
                else
                   if (interp_arr(k) == 1.0d0) then 
                      final_spectrum(k) = 10**red(j) + anew_spec(l)
                   else
                      final_spectrum(k) = 10**red(j) 
                   end if
                end if
             else
                final_spectrum(k) = 10**red(j)
             endif
             j = j+1
          end if
       else
          if (alam_ir(j) < alamb(lmax)) then
             if (anew_spec(l-1) > 0.0d0) then 
                final_spectrum(k) = 10**red(j) + anew_spec(l-1)*(anew_spec(l)/anew_spec(l-1))**interp_arr(k)
             else
                if (interp_arr(k) == 1.0d0) then 
                   final_spectrum(k) = 10**red(j) + anew_spec(l)
                else
                   final_spectrum(k) = 10**red(j) 
                end if
             end if
          else
             final_spectrum(k) = 10**red(j)
          endif
          j = j+1
       endif
       k = l+j-1
    end do

    return

  end subroutine build_spectrum

!*****************************************************************************************************************
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
end module SPECTRA_BARYONS
