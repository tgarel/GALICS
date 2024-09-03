module IO_BARYONS
  
  ! this module contains subroutines for input/output of data files.  
  
  ! line 684 : comment/uncomment 512 or 1024

  use INIT_BARYONS
  
  public


contains
  
  !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  !*****************************************************************************************************************
  subroutine init_put_baryons
    
    ! calls several subs that read in info from the input directory. 
    
    implicit none
    
    call read_baryons_datfile
    call read_filter_names   
    call read_cooling_curves 
    call read_eta            
    call read_ejecta         
    call read_spectra        
    call prepare_spectra
    call read_sfh_list
    
    return 
    
  end subroutine init_put_baryons
  
  !*****************************************************************************************************************
  subroutine  read_baryons_datfile
    
    ! read the file baryons.dat which contains all the user-adjustable parameters
    ! see GLOB_DEFS.F for a definition of these parameters
    
    implicit none                       
    
    character(MAXPATHSIZE)  :: file,line,name,value  
    integer(kind=4)         :: i
    
    write(errunit,*) '> Using data directory: '
    write(errunit,*) '> --------------------- '    
    write(errunit,*) '> ',trim(adjustl(data_dir))
    write(errunit,*)
    write(errunit,*) '> Parameter values set in GLOB_DEFS.F: '
    write(errunit,*) '> ------------------------------------ ' 
    write(errunit,*) '> Halos density profile : ',trim(adjustl(profile))
    
#ifdef CLUMPY_SF
    frac_acc_to_burst = 0.2d0    
#endif
    omega_0        = 0.333d0                                ! first, set up default values for LCDM 256^3 simu  
    hubble         = 0.667d0
    omega_b        = 0.045d0                                ! parameters are defined in module GLOB_DEFS.F
    Lbox_phys      = 150.0d0
    alphapar       = 0.033d0
    epsilon        = 0.100d0
    nu             = 0.000d0
    psi            = 0.017d0
    chi            = 3.333d0
    delay_fountain = 1.0
    duration_fountain = 1.0
    imfname        = 'salp'
#ifdef DUAL_IMF
    imfname2       = 'TH8'
#if (DUAL_IMF == 1)
    tdyn_threshold  = 16.28d0
    z_threshold = 3.0d0
#elif (DUAL_IMF == 2)
    z_threshold = 3.0d0
    m_threshold = 1.0d0
#elif (DUAL_IMF == 3)
    write(errunit,*) 'using option 3 for IMF ... which is really for testing only ... '
#else 
    write(errunit,*) 'Please set DUAL_IMF=x to a proper value ... '
    stop
#endif
#endif
    nsteps_do      = 0
    z_reionisation = 10.0d0
    alpha_reheat   = 0.1d0
    
    logMmin = -0.5d0 
    logMmax = 1.5d0
    ColdAccretionZ = 0.0

    ! initialization of random number generator 
    iseed         = -100 
    output_fits_spectra = .false.
    
#ifdef MOMAF_INPUTS
    momaf_snap_dir = ''
#endif
    
#ifdef RENUMBERING
    global_t_igm = 1.0e4
#endif
    
    write(file,'(a,a)') trim(data_dir),'/baryons.dat'  
    call test_stop(file)                                 ! test for files existence.  
    open(unit=baryons_unit,status='old',file=file)
    
    write(errunit,*)
    write(errunit,*) '> Input parameters from baryons.dat file:'
    write(errunit,*) '> ---------------------------------------'
    do
       read (baryons_unit,'(a)',end=2) line                        ! read each line as a string ...
       i     = scan(line,'=')                            ! line either begins with a '#' - ie its a comment
       if (i == 0 .or. line(1:1) == '#') cycle           ! or has parameter = value ! comment       
       name  = trim(adjustl(line(:i-1)))              
       value = trim(adjustl(line(i+1:)))
       i     = scan(value,'!')                           ! check for a comment ('!') at end of line 
       if (i /= 0) value = trim(adjustl(value(:i-1)))
       write(errunit,*) '> ',name(1:21),' : ',trim(value) ! write parameters to log file 
       
       select case (trim(name))
       case ('omega_0','Omega_0')                         ;  read(value,*) omega_0
       case ('hubble','h','hub')                          ;  read(value,*) hubble
       case ('omega_b','Omega_b')                         ;  read(value,*) omega_b
       case ('alphapar', 'alpha')                         ;  read(value,*) alphapar 
       case ('epsilon')                                   ;  read(value,*) epsilon
       case ('nu')                                        ;  read(value,*) nu
       case ('psi')                                       ;  read(value,*) psi
       case ('chi')                                       ;  read(value,*) chi
       case ('delay_fountain','fountain_delay')           ;  read(value,*) delay_fountain
       case ('duration_fountain','fountain_duration')     ;  read(value,*) duration_fountain
       case ('meta_stable','meta_stable_fudge')           ;  read(value,*) meta_stable_fudge
       case ('disc_instability_threshhold','alpha_disc')  ;  read(value,*) disc_instability_threshhold
       case ('treefile')                                  ;  treefile   = trim(value)
       case ('inputpath')                                 ;  inputpath  = trim(value)
       case ('filterpath')                                ;  filterpath = trim(value)          
       case ('spec_dir')                                  ;  spec_dir   = trim(value)
       case ('sfh_dir')                                   ;  sfh_dir    = trim(value)
       case ('imf','IMF')                                 ;  read(value,*) imfname
#ifdef CLUMPY_SF
       case ('frac_acc_to_burst')                         ;  read(value,*) frac_acc_to_burst
#endif
#ifdef DUAL_IMF
       case ('second_imf','imf2','IMF2','imfname2')             ;  imfname2   = trim(value)
#if (DUAL_IMF == 1)
       case ('tdyn_threshold','t_threshold','tdyn-threshold')   
          read(value,*) tdyn_threshold   
          tdyn_threshold = tdyn_threshold/1000.0d0 ! conversion to the right units, tdyn is entered in Myrs and the times in Galics are in Gyrs    
       case ('z_threshold','redshift_threshold','z-threshold')  ;  read(value,*) z_threshold   
#elif (DUAL_IMF == 2)
       case ('m_threshold','mass_threshold','m-threshold')      ;  read(value,*) m_threshold   
       case ('z_threshold','redshift_threshold','z-threshold')  ;  read(value,*) z_threshold   
#endif
#endif
       case ('nsteps','nsteps_do')                        ;  read(value,*) nsteps_do       
       case ('output_fits','output_fits_spectra')         ;  read(value,*) output_fits_spectra       
       case ('z_reionisation','z_reionization')           ;  read(value,*) z_reionisation       
       case ('alpha_reheat')                              ;  read(value,*) alpha_reheat
       case ('Lbox','Lbox_phys')                          ;  read(value,*) Lbox_phys 
#ifdef RENUMBERING
       case ('T_IGM','t_igm', 'T_igm') ;  read(value,*) global_t_igm
#else
       case ('T_IGM','t_igm', 'T_igm') ;  write(errunit,*) 'needs RENUMBERING compilation option' ; stop
#endif
#ifdef MOMAF_INPUTS
       case ('momaf_dir', 'momaf_snap_dir')               ; momaf_snap_dir = trim(value)
#endif

       case ('logMmin')
          read(value,*) logMmin
       case ('logMmax')
          read(value,*) logMmax
       case('ColdAccretionZ','coldAccretionZ','coldaccretionz')
          read(value,*) ColdAccretionZ

       case default                                       ;  write(errunit,*) 'dont recognise parameter:',trim(adjustl(name)),':'
       end select
       
    end do
2   close (baryons_unit)
    
#ifdef MOMAF_INPUTS 
    if (momaf_snap_dir .eq. '') stop 'Please define momaf_snap_dir in baryons.dat'
#endif
    
    return
    
  end subroutine read_baryons_datfile
  
  !*****************************************************************************************************************
  subroutine  read_filter_names
    
    ! read the filter names file.
    
    implicit none                       
    
    character(MAXPATHSIZE)  :: file
    integer(kind=4)         :: i,ierr
    
    write(file,'(a,a)') trim(data_dir),'/filters.dat'  
    call test_stop(file)                                 ! test for files existence.  
    open(unit=filters_unit,status='old',file=file)
    
    write(errunit,*)
    write(errunit,*) '> Input for spectra and metals: '
    write(errunit,*) '> ----------------------------- '
    
    read(filters_unit,*) nftot 
    allocate(filt(0:nftot),stat=ierr)
    if (ierr /= 0) stop '> problem allocating filters'
    do i = 1,nftot
       read (filters_unit,'(a)') filt(i)%name 
    end do
    close (filters_unit)
    
    return
    
  end subroutine read_filter_names
  
  !*****************************************************************************************************************
  subroutine read_cooling_curves
    
    ! read the 'cooling_curves_new.dat' file which contains metallicity-dependent cooling rates
    ! from Sutherland & Dopita (1992)
    
    implicit none              
    
    integer(kind=4)         :: ierr,i  
    character(MAXPATHSIZE)  :: file
    
    write(file,'(a,a)') trim(inputpath),'/cooling_curves_new.dat'
    open (param_unit,file=file,status='old')
    read(param_unit,*)                                           ! connect on 1st line.   
    read(param_unit,*) ncool,nmets
    allocate(templ(ncool),coolcurlmet(nmets,ncool),coolcurl(ncool),tabmetcool(nmets),stat=ierr)
    if (ierr /= 0 ) stop '> error in opening cooling curve file'
    read(param_unit,*) tabmetcool(1:nmets)                       ! log of the metallicity
    do i=1,ncool
       read (param_unit,*) templ(i),coolcurlmet(1:nmets,i)
    end do
    coolcurlmet = coolcurlmet+23.0d0    
    close(param_unit)
    
    return
    
  end subroutine read_cooling_curves
  
  !*****************************************************************************************************************
  subroutine read_eta
    
    ! read the eta_xxxx.dat file which contains the number of SN per solar mass of star formation. As this is
    ! only IMF dependent, xxxx is the IMF i.e. 'salp', 'scal', 'kenn' etc ...
    
    implicit none 
    
    character(MAXPATHSIZE) :: file     
 
#ifdef TH4
    write(file,'(4a)') trim(inputpath),'/eta_',trim(imfname),'.dat'  ! modif_imf : TH4 is contains a digit...
#else
    write(file,'(4a)') trim(inputpath),'/eta_',imfname,'.dat'
#endif
    write(errunit,*) '> ',trim(adjustl(file))
    print*,imfname
    open(param_unit,file=file,status='old')
    read(param_unit,*) eta_SN                                    ! sn per mass of stars formed.
    close(param_unit)
    
#ifdef DUAL_IMF
    write(file,'(4a)') trim(inputpath),'/eta_',trim(imfname2),'.dat'
    write(errunit,*) '> ',trim(adjustl(file))
    open(param_unit,file=file,status='old')
    read(param_unit,*) eta_sn2                                ! sn per mass of stars formed for second IMF
    close(param_unit)   
#endif
    
    return
    
  end subroutine read_eta
  
  !*****************************************************************************************************************
  subroutine read_ejecta 
    
    ! read the ejected_gas_xxxx.dat file which contains gas (+ metals) ejected as a function of time by an 
    ! instantaneous burst of 1 M_sun of stars as a function of initial metallicity.
    
    implicit none 
    
    character(MAXPATHSIZE)    :: file   
    real(kind=8),allocatable  :: gaztemp(:,:),ztemp(:,:)
    integer(kind=4)           :: i
    real(kind=8)              :: age
#ifdef DUAL_IMF 
    integer(kind=4)           :: nrej2,nfile2,nelements2
    real(kind=4)              :: dummy
#endif
    
#ifdef TH4
    write(file,'(4a)') trim(inputpath),'/ejected_gas_',trim(imfname),'_new.dat' ! modif_imf : TH4 is contains a digit...
# else
    write(file,'(4a)') trim(inputpath),'/ejected_gas_',imfname,'_new.dat'
#endif

    write(errunit,*) '> ',trim(adjustl(file))
    open(param_unit,file=file,status='old')
    read(param_unit,*) 
    read(param_unit,*) nrej, nfile, nelements
    nelements = 1                                        ! we set this to 1 as we only care about Z, not C, O, Fe
    allocate(ztab(nrej,nfile),gaztab(nrej,nfile),timetab(nrej),tabmetspec(nfile))
    allocate(delta_timetab(nrej))
    read(param_unit,*) tabmetspec(1:nfile)
    do i = 1,nrej
       read(param_unit,*) timetab(i), gaztab(i,1:nfile), ztab(i,1:nfile)
    end do
    close(param_unit)
    
    
    ! output timetab 
    write(file,'(a,a)') trim(data_dir), '/time_tab.dat'
    open(param_unit,file=file,status='unknown')
    write(param_unit,'(i8)') nrej
    do i = 1,nrej
       write(param_unit,'(e14.6)') timetab(i)
    end do
    close(param_unit)
    
    
    ! define age_index_transp
    age              = age_for_transp
    call locate(timetab,nrej,age,age_index_transp)
    age_index_transp = age_index_transp + 1
    write(errunit,*) '> Will start burst transport after ',timetab(age_index_transp),' Gyr'
    
    allocate(gaztemp(nrej,nfile),ztemp(nrej,nfile))
    
    do i=1,nrej-1
       gaztemp(i,:)     = gaztab(i,:) + gaztab(i+1,:)
       ztemp(i,:)       = ztab(i,:) + ztab(i+1,:) 
       delta_timetab(i) = timetab(i+1) - timetab(i)
    end do
    gaztemp(nrej,:)     =  gaztemp(nrej-1,:) 
    ztemp(nrej,:)       =  ztemp(nrej-1,:)       
    
    gaztab = 0.5d0 * gaztemp
    ztab   = 0.5d0 * ztemp
    
    deallocate(gaztemp,ztemp)
    
#ifdef DUAL_IMF
    ! do the same for top heavy IMF if required
    ! NB : this assumes same time-bins as the first imf ...  (nrej,nfile)... i put an explicit check...
    write(file,'(4a)') trim(inputpath),'/ejected_gas_',trim(imfname2),'_new.dat'
    write(errunit,*) '> ',trim(adjustl(file))
    open(param_unit,file=file,status='old')
    read(param_unit,*) 
    read(param_unit,*) nrej2,nfile2,nelements2
    if (nrej2 .ne. nrej .or. nfile2 .ne. nfile) then 
       write(errunit,*) 'ERROR : nrej/nfile differ for the two IMFs'
       write(errunit,*) 'nrej,nfile    =',nrej,nfile
       write(errunit,*) 'nrej2,nfile2  =',nrej2,nfile2
       stop
    end if
    allocate(ztab2(nrej,nfile),gaztab2(nrej,nfile))
    read(param_unit,*) tabmetspec(1:nfile)
    do i = 1,nrej
       read(param_unit,*) dummy, gaztab2(i,1:nfile), ztab2(i,1:nfile)
    end do
    close(param_unit)
    allocate(gaztemp(nrej,nfile),ztemp(nrej,nfile))
    do i=1,nrej-1
       gaztemp(i,:)     = gaztab2(i,:) + gaztab2(i+1,:)
       ztemp(i,:)       = ztab2(i,:)   + ztab2(i+1,:) 
    end do
    gaztemp(nrej,:)     =  gaztemp(nrej-1,:) 
    ztemp(nrej,:)       =  ztemp(nrej-1,:)       
    gaztab2 = 0.5d0 * gaztemp
    ztab2   = 0.5d0 * ztemp
    deallocate(gaztemp,ztemp)
#endif
    
    return
    
  end subroutine read_ejecta
  
  !*****************************************************************************************************************
  subroutine read_spectra
    
    ! read in spectra_xxxx.dat which contains the evolution of the SED of an instantaneous burst of stars as a 
    ! function of wavelength, time and initial metallicity (i.e. is the result of a simple stellar population 
    ! synthesis model ... here STARDUST 99)
    
    implicit none  
    
    character(MAXPATHSIZE)   :: file    
    integer(kind=4)          :: nfile,i,l
    real(kind=4),allocatable :: temp_sburst(:,:,:)
#ifdef DUAL_IMF
    integer(kind=4)          :: lmax2, nfile2, nrej2
#endif
    
#ifdef TH4
    write(file,'(4a)') trim(inputpath),'/spectra_',trim(imfname),'_new.dat' ! modif_imf : TH4 is contains a digit...
#else
    write(file,'(4a)') trim(inputpath),'/spectra_',imfname,'_new.dat'
#endif

    write(errunit,*) '> ',trim(adjustl(file))
    call test_stop(file)
    open(spec_unit,file=file,status='old')
    read(spec_unit,*) 
    read(spec_unit,*) lmax, nfile, nrej
    allocate(alamb(lmax),sburst(lmax,nrej,nfile),dalamb(lmax),temp_sburst(lmax, nrej, nfile))
    dalamb = 0.0d0
    read(spec_unit,*) alamb(1:lmax)
    do i = 1,nrej
       do l = 1,lmax
          read(spec_unit,*) temp_sburst(l,i,1:nfile)
          sburst(l,i,1:nfile) = real(temp_sburst(l,i,1:nfile),8)/10.0d0**log_L_sun  ! normalize to L_sun in order to avoid overflows
       end do
    end do
    dalamb(1) = 0.d0
    do l = 2,lmax                                     ! dalamb is difference between wavelengths
       dalamb(l) = (alamb(l)-alamb(l-1))*0.5d0        ! - useful to store as dlambda is needed for many 
    end do                                            ! integrations later.  *0.5 for trapezium rule....
    close(spec_unit)
    
#ifdef DUAL_IMF
    write(file,'(4a)') trim(inputpath),'/spectra_',trim(imfname2),'_new.dat'       
    write(errunit,*) '> ',trim(adjustl(file))
    call test_stop(file)
    open(spec_unit,file=file,status='old')
    read(spec_unit,*) 
    read(spec_unit,*) lmax2, nfile2, nrej2
    if (lmax2 /= lmax .or. nfile2 /= nfile .or. nrej2 /= nrej) then
       write(errunit,*) 'argggg'
       stop
    end if
    allocate(sburst2(lmax,nrej,nfile))
    if(imfname2.eq.'TH8')then
       do i = 1, lmax
          read(spec_unit,*) 
       end do
    else
       read(spec_unit,*) 
    endif
    do i = 1,nrej
       do l = 1,lmax
          read (spec_unit,*) temp_sburst(l,i,1:nfile)
          sburst2(l,i,1:nfile)   = real(temp_sburst(l,i,1:nfile),8)/10.0d0**log_L_sun  ! normalize to L_sun in order to avoid overflows
       end do
    end do
    close(spec_unit)
#endif
    
    deallocate(temp_sburst)
    
    return
    
  end subroutine read_spectra
  
  !*****************************************************************************************************************
  subroutine test_stop(file)
    
    ! test for file existence
    
    implicit none
    
    character(*)    :: file
    logical(kind=4) :: exist

    inquire(file=file,exist=exist)
    if (.not. exist) then
       write(errunit,*) '> fatal error: ',trim(adjustl(file)),'does not exist.'  
       stop
    endif
    
    return
    
  end subroutine test_stop



#ifdef ASCII_TREE
  !*****************************************************************************************************************
  subroutine read_tree_ascii
    
    ! reads the treefile (generally called "tree.dat") which a list of DM halos properties, and links to 
    ! fathers and sons (for details on how it is obtained see build_tree.F)
    
    implicit none
    
    type (halo),pointer         :: h
    integer(kind=4)             :: st,j,unitfile,ierr
#ifndef OLD_TREE_FMT
    integer(kind=4),allocatable :: dum(:)
#endif
    
    unitfile = tree_unit
    write(errunit,*)
    write(errunit,*) '> Input for dark matter halos (ascii): '  
    write(errunit,*) '> ------------------------------------- '
    call test_stop(treefile)
    open(unit=unitfile,file=treefile,status='old',form='formatted')
    
    read(unitfile,*) nsteps  
!    print*, nsteps
    if (nsteps_do == 0) then 
       nsteps_do = nsteps
    else
       !nsteps_do has been set in baryons.dat
    end if
    write(errunit,*) '> total number of time steps             : ',nsteps
    write(errunit,*) '> number of time steps that will be done : ',nsteps_do
    
    allocate(tsno(nsteps),stat=ierr)     ! tsno is the structure that contains the halo list ...
    call init_timestep
    if (ierr /= 0) then 
       write(errunit,*) '> not enough memory to allocate tsno'
       stop 
    end if
    
#ifdef OLD_TREE_FMT
    read(unitfile,*) tsno(1:nsteps)%nb_of_halos      
#else
    allocate(dum(nsteps))
    read(unitfile,*) tsno(1:nsteps)%nb_of_halos, dum(1:nsteps)
!    print*,  tsno(1:nsteps)%nb_of_halos, dum(1:nsteps)
    deallocate(dum)
#endif
    read(unitfile,*) tsno(1:nsteps)%aexp
    read(unitfile,*) tsno(1:nsteps)%omega_t
    read(unitfile,*) tsno(1:nsteps)%age_univ
print*,    tsno(1:nsteps)%age_univ

    tsno(1)%age_univm2   = max(tsno(1)%age_univ - 0.002d0,0.0d0) ! age_univ - 2Myr
    tsno(1)%age_univm10  = max(tsno(1)%age_univ - 0.01d0,0.0d0) ! age_univ - 10Myr
    tsno(1)%age_univm20  = max(tsno(1)%age_univ - 0.02d0,0.0d0) ! age_univ - 20Myr
    tsno(1)%age_univm50  = max(tsno(1)%age_univ - 0.05d0,0.0d0) ! age_univ - 50Myr
    tsno(1)%age_univm100 = max(tsno(1)%age_univ - 0.1d0,0.0d0)  ! age_univ - 100Myr or not
    do st = 2,nsteps 
       tsno(st)%age_univm2  = max(tsno(st)%age_univ - 0.002d0,tsno(st-1)%age_univ) ! age_univ - 2Myr
       tsno(st)%age_univm10  = max(tsno(st)%age_univ - 0.01d0,tsno(st-1)%age_univ) ! age_univ - 10Myr
       tsno(st)%age_univm20  = max(tsno(st)%age_univ - 0.02d0,tsno(st-1)%age_univ) ! age_univ - 20Myr
       tsno(st)%age_univm50  = max(tsno(st)%age_univ - 0.05d0,tsno(st-1)%age_univ) ! age_univ - 50Myr
       tsno(st)%age_univm100 = max(tsno(st)%age_univ - 0.1d0,tsno(st-1)%age_univ)  ! age_univ - 100Myr or not
    end do

#ifdef BIG_RUN
    write(errunit,*) '> normal simulation, no contamination to be read '
#endif
#ifndef BIG_RUN
    write(errunit,*) '> resimulation, will read halo contamination '
#endif
    
    do st=1,nsteps_do
       allocate(tsno(st)%liste_halos(0:tsno(st)%nb_of_halos))     
#ifdef RENUMBERING 
       allocate(tsno(st)%list_halos_number(0:tsno(st)%nb_of_halos))
#endif
       print*, st, tsno(st)%nb_of_halos
       do j = 1,tsno(st)%nb_of_halos
          h => tsno(st)%liste_halos(j)   ! h is an alias for current halo
          call init_halo(h)
          ! read halo information
          call read_halo_ascii(h,unitfile)
#ifdef RENUMBERING
          tsno(st)%list_halos_number(j) = tsno(st)%liste_halos(j)%my_number
#endif
       end do
#ifdef RENUMBERING
    tsno(st)%list_halos_number(0) = 0
#endif
!    print*, "st", st, "list" , tsno(st)%list_halos_number(:)
    end do


    
    close(unitfile) 
    
    write(errunit,*) 
    write(errunit,*) '> Starting to put baryons in DM halos: '
    write(errunit,*) '> ------------------------------------ '
    
    return
    
  end subroutine read_tree_ascii
  
  !*****************************************************************************************************************
#endif
  
  !*****************************************************************************************************************
  subroutine read_tree
    
    ! reads the treefile (generally called "tree.dat") which a list of DM halos properties, and links to 
    ! fathers and sons (for details on how it is obtained see build_tree.F)
    
    implicit none
    
    type (halo),pointer         :: h
    integer(kind=4)             :: st,j,unitfile,ierr
    real(kind=4),allocatable    :: values(:)
#ifndef OLD_TREE_FMT
    integer(kind=4),allocatable :: dum(:)
#endif
    
    unitfile = tree_unit
    write(errunit,*)
    write(errunit,*) '> Input for dark matter halos: '  
    write(errunit,*) '> ---------------------------- '
    call test_stop(treefile)
    open(unit=unitfile,file=treefile,status='old',form='unformatted')
    
    read(unitfile) nsteps  
    if (nsteps_do == 0) then 
       nsteps_do = nsteps
    else
       !nsteps_do has been set in baryons.dat
    end if
    write(errunit,*) '> total number of time steps             : ',nsteps
    write(errunit,*) '> number of time steps that will be done : ',nsteps_do
    
    allocate(tsno(nsteps),stat=ierr)     ! tsno is the structure that contains the halo list ...
    call init_timestep
    if (ierr /= 0) then 
       write(errunit,*) '> not enough memory to allocate tsno'
       stop 
    end if
    
#ifdef OLD_TREE_FMT
    read(unitfile) tsno(1:nsteps)%nb_of_halos      
#else
    allocate(dum(nsteps))
    read(unitfile) tsno(1:nsteps)%nb_of_halos, dum(1:nsteps)
    deallocate(dum)
#endif
    allocate(values(nsteps))
    read(unitfile) values(1:nsteps)
    tsno(1:nsteps)%aexp = values
    read(unitfile) values(1:nsteps)
    tsno(1:nsteps)%omega_t = values
    read(unitfile) values(1:nsteps)
    tsno(1:nsteps)%age_univ = values
    deallocate(values)
    
    tsno(1)%age_univm10  = max(tsno(1)%age_univ - 0.01d0,0.0d0) ! age_univ - 10Myr
    tsno(1)%age_univm100 = max(tsno(1)%age_univ - 0.1d0,0.0d0)  ! age_univ - 100Myr or not
    do st = 2,nsteps 
       tsno(st)%age_univm10  = max(tsno(st)%age_univ - 0.01d0,tsno(st-1)%age_univ) ! age_univ - 10Myr
       tsno(st)%age_univm100 = max(tsno(st)%age_univ - 0.1d0,tsno(st-1)%age_univ)  ! age_univ - 100Myr or not
    end do

#ifdef BIG_RUN
    write(errunit,*) '> normal simulation, no contamination to be read '
#else
    write(errunit,*) '> resimulation, will read halo contamination '
#endif
    
    do st=1,nsteps_do
       allocate(tsno(st)%liste_halos(0:tsno(st)%nb_of_halos))     
#ifdef RENUMBERING 
       allocate(tsno(st)%list_halos_number(0:tsno(st)%nb_of_halos))
#endif
       do j = 1,tsno(st)%nb_of_halos
          h => tsno(st)%liste_halos(j)   ! h is an alias for current halo
          call init_halo(h)
          ! read halo information
          call read_halo(h,unitfile)
#ifdef RENUMBERING
          tsno(st)%list_halos_number(j) = tsno(st)%liste_halos(j)%my_number
#endif
       end do
    end do
    
    close(unitfile) 
    
    write(errunit,*) 
    write(errunit,*) '> Starting to put baryons in DM halos: '
    write(errunit,*) '> ------------------------------------ '
    
    return
    
  end subroutine read_tree
  
  !*****************************************************************************************************************
  subroutine read_halo(h,unitfile)
    
    ! reads properties and links to fathers and sons of a DM halo 
    
    implicit none
    
    integer(kind=4)             :: unitfile, k
#ifndef OLD_TREE_FMT
    integer(kind=4)             :: level, hosthalo, hostsub, nbsub, nextsub
#endif
    type (halo)                 :: h 
    integer(kind=4),allocatable :: tabint(:)
    real(kind=4),allocatable    :: tabreal(:)
    real(kind=4)                :: value(4)
    
    read(unitfile) h%my_number
#ifdef DEFINE_IDS 
    read(unitfile) h%BushID
#endif
    read(unitfile) h%my_timestep   

!-------------- 1024 -----------------------------------
    h%my_form_aexp = 0.0 
!-------------- 512/256 ------------------------------------
  ! read(unitfile) value(1)
  ! h%my_form_aexp = value(1)
!-------------------------------------------------------

#ifndef OLD_TREE_FMT
    read(unitfile) level, hosthalo, hostsub, nbsub, nextsub
#endif
    read(unitfile) value(1)
    h%mfof = value(1)
#ifdef READ_HALO_ACCRETION
    read(unitfile) h%macc  ! macc is written as kind=8 : no need for "value" thing
#endif
    read(unitfile) value(1:3)
    h%p%x = value(1)
    h%p%y = value(2)
    h%p%z = value(3)
    read(unitfile) value(1:3)
    h%v%x = value(1)
    h%v%y = value(2)
    h%v%z = value(3)
    read(unitfile) value(1:3)
    h%L%x = value(1)
    h%L%y = value(2)
    h%L%z = value(3)
    read(unitfile) value(1:4)
    h%rfof = value(1)
    h%sh%a = value(2)
    h%sh%b = value(3)
    h%sh%c = value(4)
    read(unitfile) value(1:3)
    h%ek = value(1)
    h%ep = value(2)
    h%et = value(3)
    read(unitfile) value(1)
    h%spin = value(1)
    
    read(unitfile) h%my_fathers%nb_fathers
    if (h%my_fathers%nb_fathers > 0) then
       allocate(tabint(h%my_fathers%nb_fathers),tabreal(h%my_fathers%nb_fathers), & 
            & h%my_fathers%list_fathers(h%my_fathers%nb_fathers),                 & 
            & h%my_fathers%mass_fathers(h%my_fathers%nb_fathers))
       read(unitfile) tabint(1:h%my_fathers%nb_fathers)
       do k=1,h%my_fathers%nb_fathers
          h%my_fathers%list_fathers(k) = tabint(k)
       end do
       read(unitfile) tabreal(1:h%my_fathers%nb_fathers)
       do k=1,h%my_fathers%nb_fathers
          h%my_fathers%mass_fathers(k) = tabreal(k)
       end do
       deallocate(tabint,tabreal)
    endif
    

    read(unitfile) h%my_sons%nb_sons
    if (h%my_sons%nb_sons > 0) then
       allocate(tabint(h%my_sons%nb_sons),h%my_sons%list_sons(h%my_sons%nb_sons))
       read(unitfile) tabint(1:h%my_sons%nb_sons) 
       do k=1,h%my_sons%nb_sons
          h%my_sons%list_sons(k) = tabint(k)
       end do
       deallocate(tabint)
    endif
    
    read(unitfile) value(1:4)
    h%datas%rvir = value(1)
    h%datas%mvir = value(2)
    h%datas%tvir = value(3)
    h%datas%cvel = value(4)
    read(unitfile) value(1:2)
    h%halo_profile%rho_0 = value(1)
    h%halo_profile%r_c = value(2)
#ifndef BIG_RUN
    read(unitfile) h%ncont
#endif
    
    return
    
  end subroutine read_halo


#ifdef ASCII_TREE
  !*****************************************************************************************************************
  subroutine read_halo_ascii(h,unitfile)
    
    ! reads properties and links to fathers and sons of a DM halo 
    
    implicit none
    
    integer(kind=4)             :: unitfile,k,level, hosthalo, hostsub, nbsub, nextsub
    type (halo)                 :: h 
    integer(kind=4),allocatable :: tabint(:)
    real(kind=4),allocatable    :: tabreal(:)
    
    read(unitfile,*) h%my_number
!    print*,  "my bumber", h%my_number
#ifdef DEFINE_IDS
    read(unitfile,*) h%BushID
#endif
    read(unitfile,*) h%my_timestep   
!    print*, h%my_timestep   
!    read(unitfile,*) h%my_form_aexp
#ifndef OLD_TREE_FMT
    read(unitfile,*) level, hosthalo, hostsub, nbsub, nextsub
#endif
    read(unitfile,*) h%mfof
#ifdef READ_HALO_ACCRETION
    print(*,*) "can not read mass accretion in ascii file"
    stop
#endif

!    print*, "FOF", h%mfof
    read(unitfile,*) h%p%x,h%p%y,h%p%z
    read(unitfile,*) h%v%x,h%v%y,h%v%z
    read(unitfile,*) h%L%x,h%L%y,h%L%z 
    read(unitfile,*) h%rfof, h%sh%a, h%sh%b, h%sh%c
    read(unitfile,*) h%ek,h%ep,h%et
    read(unitfile,*) h%spin


    read(unitfile,*) h%my_fathers%nb_fathers
    if (h%my_fathers%nb_fathers > 0) then
       allocate(tabint(h%my_fathers%nb_fathers),tabreal(h%my_fathers%nb_fathers), & 
            & h%my_fathers%list_fathers(h%my_fathers%nb_fathers),               & 
            & h%my_fathers%mass_fathers(h%my_fathers%nb_fathers))
       read(unitfile,*) tabint(1:h%my_fathers%nb_fathers)
       do k=1,h%my_fathers%nb_fathers
          h%my_fathers%list_fathers(k) = tabint(k)
       end do
!       print*, "n papas", h%my_fathers%nb_fathers
!       print*, "lista papas", h%my_fathers%list_fathers
       read(unitfile,*) tabreal(1:h%my_fathers%nb_fathers)
       do k=1,h%my_fathers%nb_fathers
          h%my_fathers%mass_fathers(k) = tabreal(k)
       end do
       deallocate(tabint,tabreal)
    endif
    
    read(unitfile,*) h%my_sons%nb_sons
    if (h%my_sons%nb_sons > 0) then
       allocate(tabint(h%my_sons%nb_sons),h%my_sons%list_sons(h%my_sons%nb_sons))
       read(unitfile,*) tabint(1:h%my_sons%nb_sons) 
       do k=1,h%my_sons%nb_sons
          h%my_sons%list_sons(k) = tabint(k)
       end do
       deallocate(tabint)
    endif



!    print*, "number sons", h%my_sons%nb_sons
!    print*, "sons",           h%my_sons%list_sons
    read(unitfile,*) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    read(unitfile,*) h%halo_profile%rho_0,h%halo_profile%r_c
!    print*,  h%halo_profile%rho_0,h%halo_profile%r_c


#ifndef BIG_RUN
    read(unitfile,*) h%ncont
#endif
    
    return
    
  end subroutine read_halo_ascii
#endif
  
  !*****************************************************************************************************************
  subroutine write_halo(h,unitfile)
    
    ! write properties and links to fathers and sons of a DM halo into a file for debugging 
    
    implicit none
    
    integer(kind=4)             :: unitfile
    type (halo)                 :: h 
    
    write(unitfile,*) '------------------------------------------------------------------------------'  
    write(unitfile,*) 'halo num',h%my_number
    write(unitfile,*) 'halo ts  ',h%my_timestep   
    write(unitfile,*) 'halo aexp ',h%my_form_aexp
    write(unitfile,*) 'halo fof mass and radius',h%mfof,h%rfof
    write(unitfile,*) 'halo pos ',h%p%x,h%p%y,h%p%z
    write(unitfile,*) 'halo speed ',h%v%x,h%v%y,h%v%z
    write(unitfile,*) 'halo angmom ',h%L%x,h%L%y,h%L%z 
    write(unitfile,*) 'halo shape ',h%sh%a, h%sh%b, h%sh%c
    write(unitfile,*) 'halo ener ',h%ek,h%ep,h%et
    write(unitfile,*) 'halo spin ',h%spin
    write(unitfile,*) 'halo nbfas ',h%my_fathers%nb_fathers
    if (h%my_fathers%nb_fathers /= 0) then
       write(unitfile,*) 'halo falist ',h%my_fathers%list_fathers(1:h%my_fathers%nb_fathers)
       write(unitfile,*) 'halo mfalist ',h%my_fathers%mass_fathers(1:h%my_fathers%nb_fathers)
    endif
    write(unitfile,*) 'halo frag ',h%tree%frag  
    write(unitfile,*) 'halo nbdad ',h%tree%ndads  
    if (h%tree%ndads /= 0) then
       write(unitfile,*) 'halo dadlist ',h%tree%dads(1:h%tree%ndads)
    endif
    write(unitfile,*) 'halo nbfils ',h%my_sons%nb_sons
    if (h%my_sons%nb_sons /= 0) then
       write(unitfile,*) 'halo filslist ',h%my_sons%list_sons(1:h%my_sons%nb_sons) 
    endif
    write(unitfile,*) 'halo nbsons ',h%tree%nsons  
    if (h%tree%nsons /= 0) then
       write(unitfile,*) 'halo sonlist ',h%tree%sons(1:h%tree%nsons)
    endif
    write(unitfile,*) 'halo virprops ',h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    write(unitfile,*) 'halo prof ',h%halo_profile%rho_0,h%halo_profile%r_c
#ifndef BIG_RUN
    write(unitfile,*) 'halo cont ',h%ncont
#endif
    
    return
    
  end subroutine write_halo
  
  !*****************************************************************************************************************
  subroutine output_ts_info(st)
    
    implicit none
    
    character(MAXPATHSIZE)  :: fname
    integer(kind=4)         :: st

    if (st == 1) then
       write(fname,'(a,a)') trim(data_dir),'/halo-gas_count.dat'
       open(unit=hg_cnt_unit,file=fname,status='unknown')
       write(hg_cnt_unit,'(i3)') nsteps 
       
       write(fname,'(a,a)') trim(data_dir),'/gal_count.dat'
       open(unit=gcnt_unit,file=fname,status='unknown')  
       write(gcnt_unit,'(i3)') nsteps 
       
       write(fname,'(a,a)') trim(data_dir),'/SFR.dat'
       open(unit=sfr_unit,file=fname,status='unknown')
       write(sfr_unit,'(i3)') nsteps 
       
       write(fname,'(a,a)') trim(data_dir),'/merging.dat'
       open(unit=merging_unit,file=fname,status='unknown')
       write(merging_unit,'(i3)') nsteps 
       
       write(fname,'(a,a)') trim(data_dir),'/bookkeep.dat'
       open(unit=bookkeep_unit,file=fname,status='unknown')
       write(bookkeep_unit,'(i3)') nsteps 
       
       write(fname,'(a,a)') trim(data_dir),'/tau.dat'
       open(unit=tau_unit,file=fname,status='unknown')
       write(tau_unit,'(i3)') nsteps
       
       write(fname,'(a,a)') trim(data_dir),'/T_igm.dat'
       open(unit=tigm_unit,file=fname,status='unknown')
       write(tigm_unit,'(i3)') nsteps
    else
       write(fname,'(a,a)') trim(data_dir),'/halo-gas_count.dat'
       open(unit=hg_cnt_unit,file=fname,status='unknown',position='append')
       
       write(fname,'(a,a)') trim(data_dir),'/gal_count.dat'
       open(unit=gcnt_unit,file=fname,status='unknown',position='append')  
       
       write(fname,'(a,a)') trim(data_dir),'/SFR.dat'
       open(unit=sfr_unit,file=fname,status='unknown',position='append')
       
       write(fname,'(a,a)') trim(data_dir),'/merging.dat'
       open(unit=merging_unit,file=fname,status='unknown',position='append')
       
       write(fname,'(a,a)') trim(data_dir),'/bookkeep.dat'
       open(unit=bookkeep_unit,file=fname,status='unknown',position='append')
       
       write(fname,'(a,a)') trim(data_dir),'/tau.dat'
       open(unit=tau_unit,file=fname,status='unknown',position='append')
       
       write(fname,'(a,a)') trim(data_dir),'/T_igm.dat'
       open(unit=tigm_unit,file=fname,status='unknown',position='append')
    end if
    
    write(hg_cnt_unit,'(i3,1x,i6,1x,e14.6)') st,tsno(st)%nb_of_halos,tsno(st)%aexp 
    write(gcnt_unit,'(i3,1x,i6,2(1x,e14.6))') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
    write(sfr_unit,'(i3,9(1x,e14.6))') st,tsno(st)%global_SF,tsno(st)%mass_bary,tsno(st)%mass_cold, & 
         & tsno(st)%mass_star,tsno(st)%mass_metals,tsno(st)%mass_halo,tsno(st)%mass_gal_gas, &
         & tsno(st)%mass_hot_gas, tsno(st)%mass_hot_mets
    write(merging_unit,'(i3,2(1x,i6),1x,i6)') st,tsno(st)%n_ss,tsno(st)%n_dynfric,tsno(st)%n_unst
    write(bookkeep_unit,'(i3,8(1x,i6))') st,tsno(st)%n_total,tsno(st)%n_form,(tsno(st)%n_ss+tsno(st)%n_dynfric), &
         & tsno(st)%n_lost,tsno(st)%n_biglost,tsno(st)%nb_maj_mergs,tsno(st)%nb_min_mergs,tsno(st)%nb_halo_mergs
    write(tau_unit,'(i3,1x,e14.6)') st,tsno(st)%average_tau  
#ifdef TIBO
    write(tigm_unit,'(i3,3(1x,e14.6))') st,1.d0/tsno(st)%aexp-1.0d0,tsno(st)%t_igm,tsno(st)%age_univ
#else
    write(tigm_unit,'(i3,3(1x,e14.6))') st,tsno(nsteps)%aexp/tsno(st)%aexp-1.0d0,tsno(st)%t_igm,tsno(st)%age_univ
#endif    

    close(hg_cnt_unit)
    close(gcnt_unit)
    close(sfr_unit)
    close(merging_unit)
    close(bookkeep_unit)
    close(tau_unit)
    close(tigm_unit)
    
    return
    
  end subroutine output_ts_info
  
  !*****************************************************************************************************************
  subroutine output_galaxy_info(st)
    
    implicit none
    
    character(MAXPATHSIZE)   :: fname
    integer(kind=4)          :: st,j,k,i,ig   
    !tibo
    integer(kind=4)          :: do_output
    real(kind=8)             :: m,mz
    !tibo
    real(kind=8)             :: minrad,newrad
#ifdef DESIRED_OUTPUTS
    integer(kind=4)          :: i_st
#endif

! tibo
    write(fname,'(a,a,i3.3)') trim(data_dir),'/numbers.',st !will contain ts,hno,gno for each gal
    open(numbers_unit,file=fname,status='unknown')

#ifdef DEFINE_IDS
    write(fname,'(a,a,i3.3)') trim(data_dir),'/GalaxyIDs.',st !will contain ts,hno,gno for each gal (and be erase by routine define_ids)                                                                                  
    open(ids_unit,file=fname,status='unknown')
#endif
    
    do_output = 0

#ifdef DESIRED_OUTPUTS
    !allocate(desired_ts(ndes_ts))
    do i_st = 1,ndes_ts
       if (st .ne. desired_ts(i_st)) then
          do_output = 1
       else
           do_output = 0
          ! write(errunit,'(a)') 'it does happen, do_output is zero...''
           go to 3
        end if
    end do
    !deallocate(desired_ts)
#endif
    
3   continue

    if (do_output == 0) then ! global output, i.e output "all" info
       write(errunit,'(a3,i3)') 'Do output at ts = ',st
       write(fname,'(a,a,i3.3)') trim(data_dir),'/gal_results.',st
       open(unit=galres_unit,file=fname,status='unknown')
!tibo
#ifndef MIN_IO
       write(fname,'(a,a,i3.3)') trim(data_dir),'/cores.',st
       open(cores_unit,file=fname,status='unknown')
#endif
       write(fname,'(a,a,i3.3)') trim(data_dir),'/gal_sfr.',st
       open(galsfr_unit,file=fname,status='unknown')
       
!#ifdef DEFINE_IDS
!       write(fname,'(a,a,i3.3)') trim(data_dir),'/GalaxyIDs.',st !will contain ts,hno,gno for each gal (and be erase by routine define_ids)
!       open(ids_unit,file=fname,status='unknown')
!#endif
       
       
       write(galres_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
       !tibo
#ifndef MIN_IO
       write(cores_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
#endif
       write(galsfr_unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') st,tsno(st)%n_total,tsno(st)%aexp,tsno(st)%age_univ
       
       do j=1,tsno(st)%nb_of_halos
          
          call det_inclination(tsno(st)%liste_halos(j))
          call det_gal_positions(tsno(st)%liste_halos(j))
          
          do k=1,tsno(st)%liste_halos(j)%datas%nbgal         
             
             !tibo : set lambdaRvir /sqrt(2) as galaxy disc radius if central
             call get_central_galaxy(tsno(st)%liste_halos(j),ig)
             if (k == ig) then
                newrad     = tsno(st)%liste_halos(j)%spin * tsno(st)%liste_halos(j)%datas%rvir / root2
                newrad     = max(newrad,min_size)   ! limit minimum size to 100pc for extremely low-spin haloes
                !if (tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%mgal > 0.0d0) then
                tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%r_output = newrad
                !else ! this disc is new
                !tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%r_output =  tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%rgal
                !end if
             else      ! this is a satellite                     
                tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%r_output = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%rgal
             end if
             
             !end tibo
             
             call advanced_out_gal(galres_unit,tsno(st)%liste_halos(j)%datas%liste_galaxies(k), &
                  & tsno(st)%liste_halos(j)%my_number, tsno(st)%liste_halos(j)%datas%cvel,tsno(st)%liste_halos(j)%datas%mvir,tsno(st)%liste_halos(j)%datas%rvir,tsno(st)%liste_halos(j)%mfof)
             
             
             ! write out SFRs :
             write(galsfr_unit,'(9e14.6)') tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr1, & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr1,   & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr1,   &
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr10,   & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr10,  & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr10,  & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr100,  & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr100, & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr100
             ! call out_QSO(41,tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%QSO)
             !tibo
#ifndef MIN_IO
             write(cores_unit,'(2e14.6)') tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%r, & 
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%core%mass
#endif
             ! output numbers of gals :
             write(numbers_unit,'(3i8)') st,tsno(st)%liste_halos(j)%my_number,                 &
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%my_number
#ifdef DEFINE_IDS
             !          write(ids_unit,'(i8,1x,i16)') tsno(st)%liste_halos(j)%BushID, tsno(st)%liste_halos(j)%HaloID
             write(ids_unit,*) tsno(st)%liste_halos(j)%BushID, tsno(st)%liste_halos(j)%HaloID
#endif
             
          enddo
       end do
       
       close(galres_unit) ; close(galsfr_unit) ; close(numbers_unit) 
       
       
       !tibo
#ifndef MIN_IO
       close(cores_unit) 
#endif
       
#ifdef DEFINE_IDS
       close(ids_unit)
#endif

    else ! no global output, only numbers and IDs...
       do j=1,tsno(st)%nb_of_halos
          do k=1,tsno(st)%liste_halos(j)%datas%nbgal 
             ! output numbers of gals :
             write(numbers_unit,'(3i8)') st,tsno(st)%liste_halos(j)%my_number,                 &
                  & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%my_number
#ifdef DEFINE_IDS
             !          write(ids_unit,'(i8,1x,i16)') tsno(st)%liste_halos(j)%BushID, tsno(st)%liste_halos(j)%HaloID
             write(ids_unit,*) tsno(st)%liste_halos(j)%BushID, tsno(st)%liste_halos(j)%HaloID
#endif
             
          enddo
       end do
       close(numbers_unit) 
#ifdef DEFINE_IDS
       close(ids_unit)
#endif
    end if

    return
    
  end subroutine output_galaxy_info
  
  !*****************************************************************************************************************
  subroutine advanced_out_gal(unit,gal,halono,cvel,mvir,rvir,mfof)
    
    ! output data for each component of a galaxy: bulge, burst & disc
    
    implicit none
    
    real(kind=8)    :: cvel,mvir,rvir,mfof,mwind_gal,mwind_gal_z
    integer(kind=4) :: unit,halono
    type(galaxy)    :: gal                    
    
!--------------------- wind for each galaxy -------------------------------

    mwind_gal = 0.0d0
    mwind_gal = accretion_remaining_mass(gal%wind,gal%tgal)
    mwind_gal_z = accretion_remaining_mass_z(gal%wind,gal%tgal)
!--------------------------------------------------------------------------

    write(unit,'(12(e14.6,1x))',advance='no')                                                             &
         &    mwind_gal, mwind_gal_z,    gal%disc%mgal,       gal%disc%mcold,      gal%disc%sfr1,     gal%disc%minstar,       &
         & gal%disc%mcoldz,     gal%disc%rgal,      gal%disc%r_output,      gal%disc%tdyn,     gal%disc%speed,         &
         & gal%disc%transp 
    
    if (gal%bulge%mgal > 0.0d0 .or. gal%burst%mgal > 0.0d0) then 
       write(unit,fmt=11,advance='no')                                                          &
            & gal%disturb, gal%r, halono, cvel, gal%inclination, gal%nb_merg, gal%tbirth,mvir,rvir,mfof
       write(unit,fmt=12,advance='no')                                                          &
            & gal%bulge%mgal,      gal%bulge%mcold,     gal%bulge%sfr1 ,      gal%bulge%minstar, &
            & gal%bulge%mcoldz,    gal%bulge%rgal,      gal%bulge%tdyn,      gal%bulge%speed,   &
            & gal%bulge%transp
       write(unit,fmt=12,advance='yes')                                                         &
            & gal%burst%mgal,      gal%burst%mcold,     gal%burst%sfr1 ,      gal%burst%minstar, &
            & gal%burst%mcoldz,    gal%burst%rgal,      gal%burst%tdyn,      gal%burst%speed,   &
            & gal%burst%transp 
    else
       write(unit,fmt=11,advance='yes')                                                         &
            & gal%disturb, gal%r, halono, cvel, gal%inclination, gal%nb_merg, gal%tbirth,mvir,rvir,mfof
    end if
    
11  format (2(e14.6,1x),i7,1x,2(e14.6,1x),i6,1x,4(e14.6,1x)) 
12  format (9(e14.6,1x))
    
    return
    
  end subroutine advanced_out_gal
  
  !*****************************************************************************************************************
  subroutine output_halo_info(st)
    
    ! Output the halo (baryonic) info......              
    
    implicit none
    
    character(MAXPATHSIZE)  :: fname
    integer(kind=4)         :: st,j
    real(kind=8)            :: mtemp,fthtemp
    type(halo)              :: h
    
    write(fname,'(a,a,i3.3)') trim(data_dir),'/halo-gas_results.',st
    open(unit=hgres_unit,file=fname,status='unknown')

!tibo
#ifndef MIN_IO
    write(fname,'(a,a,i3.3)') trim(data_dir),'/halo_gas-profiles.',st
    open(unit=hgprof_unit,file=fname,status='unknown')
    write(fname,'(a,a,i3.3)') trim(data_dir),'/halo-reservoir.',st
    open(unit=hreserv_unit,file=fname,status='unknown')
    write(fname,'(a,a,i3.3)') trim(data_dir),'/halo_dm_results.',st
    open(unit=hdmres_unit,file=fname,status='unknown')
    write(hdmres_unit,'(i3,1x,i6,1x,e14.6)') st,tsno(st)%nb_of_halos,tsno(st)%aexp 
    write(hgprof_unit,'(i3,1x,i6,1x,e14.6)') st,tsno(st)%nb_of_halos,tsno(st)%aexp 
    write(hreserv_unit,'(i3,1x,i6,1x,e14.6)') st,tsno(st)%nb_of_halos,tsno(st)%aexp 
#endif
    write(hgres_unit,'(i3,1x,i6,1x,e14.6)') st,tsno(st)%nb_of_halos,tsno(st)%aexp 

    
    do j=1,tsno(st)%nb_of_halos
       
       h = tsno(st)%liste_halos(j)
       ! tibo: convert halo positions into comoving Mpc/h: NEEDED for lightcones of halos in MOMAF...
       ! call det_halo_positions(h)

       fthtemp = h%my_form_aexp      ! fthtemp = formation expansion factor
       mtemp   = (omega_b/omega_0) * h%mfof
       call halo_out(hgres_unit,h,mtemp,fthtemp)
!tibo
#ifndef MIN_IO
       call dm_out(hdmres_unit,h)
#endif
       
    end do
    
    close(hgres_unit)
!tibo
#ifndef MIN_IO
    close(hgprof_unit)
    close(hreserv_unit)
    close(hdmres_unit)
#endif
    
    return
    
  end subroutine output_halo_info
  
  !*****************************************************************************************************************
#ifndef BIG_RUN
  subroutine contamination_out(st)
    
    ! output data for contaminated halos and galaxies in case of resimulation
    
    implicit none
    
    integer(kind=4)         :: st,j,k
    character(MAXPATHSIZE)  :: fname

!tibo
#ifndef MIN_IO    
    write(fname,'(a,a,i3.3)') trim(data_dir),'/gal_contamination.',st
    open(contam_unit,file=fname,status='unknown')
    do j=1,tsno(st)%nb_of_halos
       if (tsno(st)%liste_halos(j)%datas%nbgal /= 0) then           
          do k=1,tsno(st)%liste_halos(j)%datas%nbgal
             write(contam_unit,'(i4)')  tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%contam
          end do
       end if
    end do
    close(contam_unit)
    
    write(fname,'(a,a,i3.3)') trim(data_dir),'/halo_contamination.',st
    open(contam_unit,file=fname,status='unknown')
    do j=1,tsno(st)%nb_of_halos
       write(contam_unit,'(i4)')  tsno(st)%liste_halos(j)%ncont
    end do
    close(contam_unit)
#endif
    
    return
    
  end subroutine contamination_out
#endif
  !*****************************************************************************************************************
  subroutine halo_out(unit,h,mtemp,fthtemp)
    
    ! output data for DM and baryons in halos
    
    implicit none
    
    integer(kind=4) :: unit
    type(halo)      :: h
    real(kind=8)    :: mtemp,fthtemp             
    
    real(kind=8)    :: metsinstars ! mass of metals that went into stars and is still there 
    real(kind=8)    :: mstar,m 
    integer(kind=4) :: i,j
    
    metsinstars = 0.0d0
    mstar       = 0.0d0
    do i = 1,h%datas%nbgal ! sum over all galaxies in the halo 
       if (h%datas%liste_galaxies(i)%disc%minstar > 0.0d0) then
          do j = 1,nfile ! sum over all initial metallicities
             m           = sum(h%datas%liste_galaxies(i)%disc%sfh_tab(:,j))
             mstar       = mstar       + m
             metsinstars = metsinstars + m * tabmetspec(j)
          end do
       end if
       if (h%datas%liste_galaxies(i)%bulge%minstar > 0.0d0) then
          do j = 1,nfile ! sum over all initial metallicities
             m           = sum(h%datas%liste_galaxies(i)%bulge%sfh_tab(:,j))
             mstar       = mstar       + m
             metsinstars = metsinstars + m * tabmetspec(j)
          end do
       end if
       if (h%datas%liste_galaxies(i)%burst%minstar > 0.0d0) then
          do j = 1,nfile ! sum over all initial metallicities
             m           = sum(h%datas%liste_galaxies(i)%burst%sfh_tab(:,j))
             mstar       = mstar       + m
             metsinstars = metsinstars + m * tabmetspec(j)
          end do
       end if
    end do

!jeje 
mtemp = h%datas%mgaz
     
#ifdef RENUMBERING
    write(unit,'(i8,1x,11(e14.6,1x),i6,1x,4(e14.6,1x))')  h%my_number,                  & 
         & h%mfof          , h%datas%mvir       , h%datas%mcoldgaz     , mtemp        , & 
         & h%datas%mhotgaz , h%datas%mhotz      , h%datas%mcoldz       ,                & 
         & h%rfof          , h%datas%rvir       , h%datas%tvir         , fthtemp      , & 
         & h%datas%nbgal   , h%halo_profile%r_c , h%halo_profile%rho_0 , &
         & mstar           , metsinstars
#else    
    write(unit,'(11(e14.6,1x),i6,1x,4(e14.6,1x))')                                      & 
         & h%mfof          , h%datas%mvir       , h%datas%mcoldgaz     , mtemp        , & 
         & h%datas%mhotgaz , h%datas%mhotz      , h%datas%mcoldz       ,                & 
         & h%rfof             , h%datas%rvir         ,                & 
         & h%datas%tvir    , fthtemp            ,                                       & 
         & h%datas%nbgal   , h%halo_profile%r_c , h%halo_profile%rho_0 , &
         & mstar           , metsinstars
#endif
    
    return
    
  end subroutine halo_out
  
  !*****************************************************************************************************************
  subroutine dm_out(dmunit,h)
    ! output DM props of the halos
    
    implicit none
    
    integer(kind=4) :: dmunit
    type(halo)      :: h
    
    write(dmunit,'(13(e14.6,1x))') &
         & h%p%x,h%p%y,h%p%z, &
         & h%v%x,h%v%y,h%v%z, &
         & h%L%x,h%L%y,h%L%z, &
         & h%ek, h%ep, h%et, h%spin
    
    return
    
  end subroutine dm_out

  !*****************************************************************************************************************
  subroutine output_wave(localwave,n,word)
    
    implicit none
    
    integer(kind=4)         :: i,n
    character(10)           :: word
    real(kind=8)            :: localwave(1:n) 
    character(MAXPATHSIZE)  :: fname
    
    write(fname,'(a,a,a,a)') trim(data_dir),'/',trim(word),'_wavelength.dat'
    open(unit=wave_unit,file=fname,form='formatted',status='unknown')
    do i=1,n
       write(wave_unit,*) localwave(i)
    end do
    close(wave_unit)
    
    return
    
  end subroutine output_wave
  
  !*****************************************************************************************************************
  subroutine read_filters             
    
    ! read in the filters from filter library, and calculate normalizations 
    
    implicit none                      
    
    integer(kind=4)          :: nf,lf,ierr
    real(kind=8)             :: total,total_V
    real(kind=8),allocatable :: temp(:)
    character(MAXPATHSIZE)   :: file
    logical(kind=4)          :: exist

    total_V = 0.0d0
    
    ! first, read the JOHNSON_V.dat file to calibrate AB mags.      
    filt(0)%name = 'JOHNSON_V'
    do nf=0,nftot 
       
       write(file,'(a,a,a,a)') trim(filterpath),'/',trim(filt(nf)%name),'.dat'
       
       if (nf ==0) then 
          inquire(file=file,exist=exist)
          if (.not.exist) stop 'need the file JOHNSON_V.dat to calibrate magnitudes'
       end if
       
       open(unit=filters_unit,status='old',file=file)
       read(filters_unit,*) filt(nf)%lftot,filt(nf)%cfil,filt(nf)%aire
       allocate(filt(nf)%wave(filt(nf)%lftot),filt(nf)%trans(filt(nf)%lftot),stat=ierr)
       do lf=1,filt(nf)%lftot
          read(filters_unit,*) filt(nf)%wave(lf),filt(nf)%trans(lf) 
       enddo
       ! files are normalized such that maximum transmission = 1.0.  ie width is fwhm.    
       filt(nf)%trans = filt(nf)%trans / maxval(filt(nf)%trans) 
       ! if cfil and width are not specified, calculate them.
       if (filt(nf)%cfil == 0.0d0 .or. nf == 0) then           
          total                  = 0.0d0 
          filt(nf)%aire          = 0.0d0
          allocate(temp(1:filt(nf)%lftot))
          temp(1:filt(nf)%lftot) = filt(nf)%trans(1:filt(nf)%lftot)/filt(nf)%wave(1:filt(nf)%lftot)**2 
          total                  = trap(filt(nf)%lftot,temp(1:filt(nf)%lftot),filt(nf)%wave(1:filt(nf)%lftot))
          deallocate(temp)
          filt(nf)%aire          = trap(filt(nf)%lftot,filt(nf)%trans(1:filt(nf)%lftot), &
               &       filt(nf)%wave(1:filt(nf)%lftot))
          total                  = filt(nf)%aire/total
          if (nf == 0) then
             total_V = total
          else
             filt(nf)%cfil = filt(0)%cfil/(total/total_V)
          end if
       end if
       filt(nf)%wave(1:filt(nf)%lftot) = filt(nf)%wave(1:filt(nf)%lftot)/10000.0d0  ! 10000 does Angstroms -> microns 
       filt(nf)%aire                   = filt(nf)%aire/10000.0d0  
       close(filters_unit) 
       
    enddo
    
    ! these redefinitions help save a few calculations in the convolution code.   
    
    filt(1:nftot)%cfil = filt(1:nftot)%aire * filt(1:nftot)%cfil / filt(0)%cfil 
    filt(1:nftot)%cfil = 4.83d0 + 2.5d0 * log10(filt(1:nftot)%cfil) 
    
    return 
    
  contains
    
    !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real function trap(n,y,x)              
      
      ! simple array trapezium rule:
      
      implicit none
      
      integer(kind=4)         :: i,n
      real(kind=8),intent(in) :: x(n),y(n)
      
      trap = 0.0d0
      do i = 1,n-1
         trap = trap + (y(i) + y(i+1))* (x(i+1) - x(i))
      end do
      trap = trap*0.5d0
      
      return
      
    end function trap
    !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  end subroutine read_filters
  
  !*****************************************************************************************************************
  subroutine read_dust_info           
    
    ! read in the data from nebform.dat (extinction curves) and spec_IR_new.dat (dust emission spectra templates) 
    
    implicit none                      
    
    integer(kind=4) :: i,l
    
    open(unit=param_unit,file=trim(inputpath)//'/nebform.dat',status='old')
    read(param_unit,'(i2)') find
    !write(errunit,*) find
    allocate(alambda(find),exti(find,3))
    do i=1,find
       read(param_unit,'(4(e15.6,1x))') alambda(i),(exti(i,l),l=1,3)
    end do
    close(param_unit)

!tibo: info: find = 38 = # number of A_l/A_v(lambda) = low_res tabulation

#ifdef DALE_HELOU
    open(unit=param_unit,file=trim(inputpath)//'/spectra_DaleHelou.dat',status='old')
    write(errunit,*) trim(inputpath)//'/spectra_DaleHelou.dat'
#else
    open(unit=param_unit,file=trim(inputpath)//'/spec_IR_new.dat',status='old')    
    write(errunit,*) trim(inputpath)//'/spec_IR_new.dat'
#endif    
    read(param_unit,*)                          ! 1st line is a comment
    read(param_unit,*) fins,nbre
    !write(errunit,*) fins,nbre
    allocate(lum(nbre),alam_ir(fins),spec_ir(fins,nbre))
    read(param_unit,*) lum(1:nbre)
    lum = log10(lum)
    do i=1,fins
       read(param_unit,*) alam_ir(i),spec_ir(i,1:nbre)
    end do
    close(param_unit)
    spec_ir = log10(spec_ir)             ! taking log here saves a lot of work later
    
    return
    
  end subroutine read_dust_info
  
  !*****************************************************************************************************************
  subroutine prepare_spectra
    
    implicit none
    
    integer(kind=4)          :: i
    logical(kind=4)          :: ascii_existlist
    character(MAXPATHSIZE)   :: fname
    
    
    ! check wether there is a list of timesteps for which to output restframe mags
    write(fname,'(a,a)') trim(data_dir),'/ascii_ts_list.dat'
    inquire(file=fname,exist=ascii_existlist)
    if (ascii_existlist) then 
       open(unit=param_unit,status='old',file=fname)
       read(param_unit,*)  ! header
       read(param_unit,*) n_ascii_want
       allocate(ascii_list(n_ascii_want))
       do i = 1,n_ascii_want
          read(param_unit,*) ascii_list(i)
       end do
       close(param_unit)
    else
       n_ascii_want = 0
    endif    
    
    call output_wave(alamb,lmax,'base      ') ! output alamb to file base_wavelength.dat
    call read_dust_info                       ! define the global variables : find, alambda, exti, fins, nbre, 
    ! lum, alam_ir, spec_ir
    call read_filters                         ! define global variable filt
    allocate(interp_arr(fins+lmax),wavelength(fins+lmax),restframe_wave(lmax+fins)) ! allocate global variables 
    call setup_spectra_saves                                                        ! and define them in setup_spectra_saves
    ! allocate global vars cosmic_extinction(2) 
    allocate(cosmic_extinction(lmax+fins),cosmic_extinction2(lmax+fins))
    allocate(tot_wave(lmax+fins),tot_wave2(lmax+fins))
    allocate(amabs_diff(nftot),amabs_rf_noext(nftot),amabs_rf(nftot))
#ifdef MOMAF_INPUTS
    allocate(amabs_ce(nftot))
#endif
    write(fname,'(a,a)') trim(data_dir),'/ce_effect.dat'
    open(unit=ce_unit,file=fname,status='unknown') 
    write(ce_unit,'(i3)') nsteps

    return
    
  end subroutine prepare_spectra
  
  !*****************************************************************************************************************
  subroutine setup_spectra_saves   
    
    ! call this subroutine once to set up global variables, look-up tables etc.
    
    implicit none     
    
    integer(kind=4)          :: i,j,k,l
    real(kind=8)             :: ip1fac,res
    real(kind=8),allocatable :: tau_lores(:,:)
    
    allocate(tau_hires(lmax,3),albedo(lmax),tau_lores(find,3))

!tibo: info : lmax = 1221 = # of bins of stellar spectra = hi_res tabulation

    ! First, combine the STELLAR and DUST wavelength grids......... 
    l = 1 ;  j = 1  ;  k = 1
    do while (k <= fins+lmax)
       if (l <= lmax) then
          if (alamb(l) < alam_ir(j)) then
             wavelength(k) = alamb(l)
             if (j /=1 ) interp_arr(k) = log10(wavelength(k)/alam_ir(j-1)) / &
                  & log10(alam_ir(j)/alam_ir(j-1))
             l             = l+1
          else
             wavelength(k) = alam_ir(j)    
             if (alam_ir(j) < alamb(lmax)) interp_arr(k) = log10(wavelength(k)/alamb(l-1)) / &
                  & log10(alamb(l) / alamb(l-1))
             j             = j+1
          end if
       else
          wavelength(k) = alam_ir(j)    
          if (alam_ir(j) < alamb(lmax)) interp_arr(k) = log10(wavelength(k)/alamb(l-1)) / &
               & log10(alamb(l) / alamb(l-1))
          j             = j+1
       end if
       k = l+j-1
    end do
    
    call output_wave(wavelength,fins+lmax,'corrected ')  ! output this combined wavelength array 
    restframe_wave = wavelength
    
    ! Second, make look-up table for the optical depth: 
    call compute_optical_depth(tau_lores)

    ! Third, interpolate to get a hires grid in wavelength:  
    do i=1,lmax         
       
       call locate(alambda,find,alamb(i),j)
       if ((j < find).and.(j > 0)) then
          ! here the interpolation is logarithmic because tau_lam is given as E(lam-V)/E(B-V) which 
          ! are differences of extinction in magnitudes i.e. propto log(lambda)
          do k=1,3
             tau_hires(i,k) = xinterp(log10(alambda(j)),tau_lores(j,k),log10(alambda(j+1)),tau_lores(j+1,k), &
                  &         log10(alamb(i)))
          end do
          if (exti(j+1,2) > 0.0d0) then
             albedo(i)      = 10.0d0**xinterp(log10(alambda(j)),log10(exti(j,2)),log10(alambda(j+1)), & 
                  &               log10(exti(j+1,2)),log10(alamb(i)))             
          else
             albedo(i)      = 0.0d0             
          endif
       else
          do k=1,3
             tau_hires(i,k) = tau_lores(1,k)
          end do
          albedo(i) = exti(1,2)
       endif
       
    enddo
    
    ! Lastly, tabulate tau and the corresponding extinction:
    do i=1,ntau
       tau_tab(i) = 10**(real(i,8)/20.0d0 -4.d0)
    end do
    
    ! series expansion: x = 3*Sigma_1^infinity (2tau)^n (-1)^n / (n+1)! (n+3)
    ip1fac = 2
    do i=1,nseries        
       series(i) = 3.0d0 * (-2.0d0)**i / ((i+3) * (ip1fac))
       ip1fac    = ip1fac*(i+1)
    end do
    
    do i=1,ntau
       if (tau_tab(i) > 10.0d0) then                                                   ! Somewhere we really should 
          ext_tab(i) = (3.0d0/(4.0d0*tau_tab(i))) * (1.0d0 - (1.0d0/(2.0d0 * tau_tab(i)**2)))  ! keep track of all the 
       else                                                                          ! calculations that have gone 
          if (tau_tab(i) > 0.1d0) then                                                 ! into these expansions.  
             ext_tab(i) =  (3.0d0/(4.0d0*tau_tab(i))) * (1.0d0 - (1.0d0/(2.0d0 * tau_tab(i)**2)) + &     
                  & ( (1.0d0/tau_tab(i)) + 1.0d0/(2.0d0 * tau_tab(i)**2))*exp(-2.0d0*tau_tab(i)))
          else
             ! order in the ifs is important below bc we're adding more and more terms in the series
             ! depending on how large the value of tau_tab(i)
             ext_tab(i) = 1.0d0 + tau_tab(i)*series(1)
             if (tau_tab(i) > 0.001d0) ext_tab(i) = ext_tab(i) + series(2)*tau_tab(i)**2
             if (tau_tab(i) > 0.003d0) ext_tab(i) = ext_tab(i) + series(3)*tau_tab(i)**3
             if (tau_tab(i) > 0.01d0)  ext_tab(i) = ext_tab(i) + series(4)*tau_tab(i)**4
             if (tau_tab(i) > 0.03d0)  ext_tab(i) = ext_tab(i) + series(5)*tau_tab(i)**5
          endif
       endif
    enddo
    
    do i = 1,ntau
       if (tau_tab(i) < 75.0d0) then 
          call expint(1,tau_tab(i),res)
          slab_int_tab(i) = 1.0d0 - (0.5d0/tau_tab(i))*(1.0d0 + (tau_tab(i)-1.0d0)*exp(-tau_tab(i)) - &
               &                         tau_tab(i)**2 * res)
          slab_tab(i)     = (1.0d0 - exp(-tau_tab(i)))/tau_tab(i)
       else
          slab_int_tab(i) = 1.0d0 - (0.5d0/tau_tab(i))
          slab_tab(i)     = 1.0d0 / tau_tab(i)
       end if
    end do
    
    deallocate(tau_lores)
    
    return
    
  contains
    
    !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    subroutine compute_optical_depth(tau_lam) 
      
      ! this subroutine sets up the coarse bins of tau.
      
      implicit none
      
      real(kind=8),intent(out) :: tau_lam(find,3)
      integer(kind=4)          :: i
      
      do i=1,find
         !tau_lam(i,1) = log10(3.25d0*exti(i,1)/6.8d0) -21.d0 
         !tibo:
         tau_lam(i,1) = log10(exti(i,1)) -21.d0 
         tau_lam(i,2) = (1.6d0 + exti(i,3))
         tau_lam(i,3) = log10(6.3d0*(1/0.0900d0*alamb(i))**3) -18.d0 ! HI absortion hasn't been tested yet!
      enddo
      
      return
      
    end subroutine compute_optical_depth
    !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  end subroutine setup_spectra_saves
  
  !*****************************************************************************************************************
  subroutine open_magfile(unit,name,st,nbtot_gal)
    
    ! open ascii file to write final magnitudes of galaxies
    
    implicit none
    
    integer(kind=4), intent(in) :: unit
    integer(kind=4)             :: st,nbtot_gal
    character(MAXPATHSIZE)      :: fname
    character(*)                :: name
    
    write(fname,'(a,a,a,a1,i3.3)') trim(data_dir),'/',trim(name),'.',st
    open(unit=unit,file=fname,form='formatted',status='unknown')   
    write(unit,*) nftot,nbtot_gal 
    
    return
    
  end subroutine open_magfile
    
  !*****************************************************************************************************************
  subroutine output_gal_tree(st)
    
    implicit none
    
    character(MAXPATHSIZE)    :: fname
    integer(kind=4)           :: st,j,k
    
    write(fname,'(a,a,i3.3)') trim(data_dir),'/gal_daughters.',st
    open(unit=gal_tree_unit,file=fname,status='unknown')
    do j=1,tsno(st)%nb_of_halos
       do k=1,tsno(st)%liste_halos(j)%datas%nbgal
          write(gal_tree_unit,'(2i8)')  tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%my_girls%hno, & 
               & tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%my_girls%gno
       end do
    end do
    close(gal_tree_unit)
    
    return
    
  end subroutine output_gal_tree
  
  !*****************************************************************************************************************
  subroutine read_sfh_list
    
    ! read the list of galaxies for which we want to output the sfh_tab's. This list is in a file named 'sfh_gal_list.dat',
    ! which is sorted with ts,hno,gno.
    
    implicit none
    
    character(MAXPATHSIZE)  :: fname
    logical(kind=4)         :: existlist
    integer(kind=4)         :: i
    
    write(fname,'(a,a)') trim(data_dir),'/sfh_gal_list.dat'
    inquire(file=fname,exist=existlist)
    if (existlist) then 
       open(unit=param_unit,status='old',file=fname)
       read(param_unit,*)  ! header
       read(param_unit,*) n_sfh_gal
       allocate(sfh_gal_list(n_sfh_gal,3))
       do i = 1, n_sfh_gal
          read(param_unit,*) sfh_gal_list(i,1),sfh_gal_list(i,2),sfh_gal_list(i,3)
       end do
       close(param_unit)
    else
       n_sfh_gal = 0
    endif
    
    
    ind_sfh_gal = 1
    
    return
    
  end subroutine read_sfh_list
  
  !*****************************************************************************************************************
  subroutine output_sfh_tab(st)
    
    implicit none
    
    integer(kind=4)         :: st
    integer(kind=4)         :: hno,gno
    character(MAXPATHSIZE)  :: filename
    integer(kind=4)         :: n_age,n_mets
    integer(kind=4)         :: icomp,j,k,ierr
    
    if (st == sfh_gal_list(ind_sfh_gal,1)) then
       
       do hno = 1, tsno(st)%nb_of_halos
          do gno = 1, tsno(st)%liste_halos(hno)%datas%nbgal
             
             ! test wether galaxy is in the list to be outputed
             if (st == sfh_gal_list(ind_sfh_gal,1) .and. hno == sfh_gal_list(ind_sfh_gal,2) & 
                  & .and. gno == sfh_gal_list(ind_sfh_gal,3)) then 
                ! output one file per existing component
                do icomp = 1,3
                   if (icomp == 1) then
                      if (tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%disc%minstar > 0.0d0) then
                         write(filename,'(a,a,''sfh_disc.'',i2.2,i5.5,i3.3,''.dat'')') trim(data_dir), trim(sfh_dir), &
                              & st,hno,gno
                         open(param_unit,file=filename,status='new',iostat=ierr)  
                         n_age  = size(tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%disc%sfh_tab(:,:),dim=1)
                         n_mets = size(tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%disc%sfh_tab(:,:),dim=2)
                         write(param_unit,'(2(i8,1x))') n_age, n_mets
                         do j = 1,n_mets
                            do k = 1,n_age
                               write(param_unit,'(e14.6)') tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%disc%sfh_tab(k,j)
                            end do
                         end do
                         close(param_unit)
                      end if
                   else if (icomp == 2) then
                      if (tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%bulge%minstar > 0.0d0) then
                         write(filename,'(a,a,''sfh_bulge.'',i2.2,i5.5,i3.3,''.dat'')')  trim(data_dir), trim(sfh_dir), &
                              & st,hno,gno
                         open(unit=param_unit,file=filename,status='unknown')                      
                         n_age  = size(tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%bulge%sfh_tab(:,:),dim=1)
                         n_mets = size(tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%bulge%sfh_tab(:,:),dim=2)
                         write(param_unit,'(2(i8,1x))') n_age, n_mets
                         do j = 1,n_mets
                            do k = 1,n_age
                               write(param_unit,'(e14.6)') tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%bulge%sfh_tab(k,j)
                            end do
                         end do
                         close(param_unit)
                      end if
                   else if (icomp == 3) then
                      if (tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%burst%minstar > 0.0d0) then
                         write(filename,'(a,a,''sfh_burst.'',i2.2,i5.5,i3.3,''.dat'')') trim(data_dir),  trim(sfh_dir), &
                              & st,hno,gno
                         open(unit=param_unit,file=filename,status='unknown')                      
                         n_age  = size(tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%burst%sfh_tab(:,:),dim=1)
                         n_mets = size(tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%burst%sfh_tab(:,:),dim=2)
                         write(param_unit,'(2(i8,1x))') n_age, n_mets
                         do j = 1,n_mets
                            do k = 1,n_age
                               write(param_unit,'(e14.6)') tsno(st)%liste_halos(hno)%datas%liste_galaxies(gno)%burst%sfh_tab(k,j)
                            end do
                         end do
                         close(param_unit)
                      end if
                   endif
                end do
                
                ! increment index of the list 
                ind_sfh_gal = ind_sfh_gal + 1
                if (ind_sfh_gal > n_sfh_gal) stop
             end if
             
          end do
       end do
       
    else
       
       ! check that no galaxy was missed 
       if (st > sfh_gal_list(ind_sfh_gal,1)) then
          write(errunit,*) '> Problem in output_sfh_tab : skipped galaxies... or galaxy missing...' 
          write(errunit,*) st,ind_sfh_gal
          write(errunit,*) sfh_gal_list(ind_sfh_gal,1),sfh_gal_list(ind_sfh_gal,2),sfh_gal_list(ind_sfh_gal,3)
          stop
       end if
       
    end if
    
    return
    
  end subroutine output_sfh_tab
  
  !*****************************************************************************************************************
  subroutine write_snaps_redshifts
    
    implicit none

    integer(kind=4)         :: ts
    real(kind=8)            :: z
    character(MAXPATHSIZE)  :: filename
    
    write(filename,'(a,a)') trim(data_dir),'/snaps_redshifts.dat'
    open(unit=param_unit,file=filename,status='unknown',form='formatted')
    do ts = 1,nsteps
#ifdef TIBO
       z = 1.d0 / tsno(ts)%aexp - 1.d0
#else
       z = tsno(nsteps)%aexp / tsno(ts)%aexp - 1.d0
#endif
       write(param_unit,*) z
    end do
    close(param_unit)

#ifdef MOMAF_INPUTS
    write(filename,'(a,a)') trim(momaf_snap_dir),'/snaps_redshifts.dat'
    open(unit=param_unit,file=filename,status='unknown',form='formatted')
    do ts = 1,nsteps
#ifdef TIBO
       z = 1.d0 / tsno(ts)%aexp - 1.d0
#else
       z = tsno(nsteps)%aexp / tsno(ts)%aexp - 1.d0
#endif
       write(param_unit,*) z
    end do
    close(param_unit)
#endif

    return

  end subroutine write_snaps_redshifts

  !*****************************************************************************************************************
#ifdef MOMAF_INPUTS
  subroutine write_momaf_snap(st)
    
    implicit none

    integer(kind=4)             :: st,ih,ig,igal
    integer(kind=4)             :: nb_tot_gals
    integer(kind=8),allocatable :: galid(:),haloid(:)
    real(kind=8),allocatable    :: pos(:,:),vel(:,:),inc(:)
    character(MAXPATHSIZE)      :: filename

    nb_tot_gals        = tsno(st)%n_total
    if (nb_tot_gals > 0) then 
       allocate(pos(3,nb_tot_gals),vel(3,nb_tot_gals),galid(nb_tot_gals),haloid(nb_tot_gals),inc(nb_tot_gals))
       igal = 1
       do ih = 1, tsno(st)%nb_of_halos
          do ig = 1,tsno(st)%liste_halos(ih)%datas%nbgal
             pos(1,igal)  = tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%p%x
             pos(2,igal)  = tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%p%y
             pos(3,igal)  = tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%p%z
             vel(1,igal)  = tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%v%x
             vel(2,igal)  = tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%v%y
             vel(3,igal)  = tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%v%z
             inc(igal)    = tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%inclination
             galid(igal)  = int(igal,8)
             haloid(igal) = int(ih,8)
             igal         = igal + 1
          end do
       end do
       write(filename,'(a,a,i3.3)') trim(momaf_snap_dir),'/IDsPosVel.',st
       filename = trim(filename)//char(0)
       call write_momaf_propfile(filename,nb_tot_gals,galid,haloid,real(inc,4),real(pos,4),real(vel,4))
       deallocate(pos,vel,galid,haloid,inc)
    else
       write(filename,'(a,a,i3.3)') trim(momaf_snap_dir),'/IDsPosVel.',st
       filename = trim(filename)//char(0)
       call write_momaf_empty_propfile(filename)
    end if

    return

  end subroutine write_momaf_snap

  !*****************************************************************************************************************

  subroutine write_momaf_empty_fileset(st)

    implicit none

    integer(kind=4)         :: st,ifilt
    character(MAXPATHSIZE)  :: filename

    write(filename,'(a,a,i3.3)') trim(momaf_snap_dir),'/IDsPosVel.',st
    filename = trim(filename)//char(0)
    call write_momaf_empty_propfile(filename)
    
    do ifilt = 1, nftot
       write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/disc_',trim(filt(ifilt)%name),'.',st
       filename = trim(filename)//char(0)
       call write_momaf_empty_magfile(filename)
       write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/bulge_',trim(filt(ifilt)%name),'.',st
       filename = trim(filename)//char(0)
       call write_momaf_empty_magfile(filename)
       write(filename,'(a,a,a,a,i3.3)') trim(momaf_snap_dir),'/burst_',trim(filt(ifilt)%name),'.',st
       filename = trim(filename)//char(0)
       call write_momaf_empty_magfile(filename)
    end do

    return

  end subroutine write_momaf_empty_fileset

#endif

  !*****************************************************************************************************************
#ifdef RECORD_SFR
  subroutine output_sfr_tabs(st)
    
    implicit none 

    integer(kind=4)        :: st,j,k,n1,n2
    real(kind=8)           :: sfrtot
    character(MAXPATHSIZE) :: filename

    ! write out disc component 
    write(filename,'(a,a,i3.3)') trim(data_dir),'/disc_sfrtabs.',st
    open(unit=sfrtab_unit,file=filename,form='unformatted',status='unknown')
    write(sfrtab_unit) tsno(st)%n_total,nrej,nfile  ! nb of galaxies, nb of agebins, and nb of metallicities
    write(sfrtab_unit) timetab(1:nrej)
    write(sfrtab_unit) tabmetspec(1:nfile)
    do j=1,tsno(st)%nb_of_halos
       do k=1,tsno(st)%liste_halos(j)%datas%nbgal
          sfrtot = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%totsfr
          if (sfrtot > 0.0d0) then 
             n1 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr_tab,dim=1)
             n2 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr_tab,dim=2)
             write(sfrtab_unit) sfrtot, n1, n2                
             write(sfrtab_unit) tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr_tab(:,:)
          else 
             n1 = 0
             n2 = 0
             write(sfrtab_unit) sfrtot, n1, n2
          end if
       end do
    end do
    close(sfrtab_unit)
#ifdef DUAL_IMF 
    write(filename,'(a,a,i3.3)') trim(data_dir),'/disc_sfrtabs_imf2.',st
    open(unit=sfrtab_unit,file=filename,form='unformatted',status='unknown')
    write(sfrtab_unit) tsno(st)%n_total,nrej,nfile  ! nb of galaxies, nb of agebins, and nb of metallicities
    write(sfrtab_unit) timetab(1:nrej)
    write(sfrtab_unit) tabmetspec(1:nfile)
    do j=1,tsno(st)%nb_of_halos
       do k=1,tsno(st)%liste_halos(j)%datas%nbgal
          sfrtot = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%totsfr2
          if (sfrtot > 0.0d0) then 
             n1 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr_tab2,dim=1)
             n2 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr_tab2,dim=2)
             write(sfrtab_unit) sfrtot, n1, n2
             write(sfrtab_unit) tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%disc%sfr_tab2(:,:)          
          else
             n1 = 0 
             n2 = 0
             write(sfrtab_unit) sfrtot, n1, n2
          end if
       end do
    end do
    close(sfrtab_unit)
#endif

    ! write out bulge component 
    write(filename,'(a,a,i3.3)') trim(data_dir),'/bulge_sfrtabs.',st
    open(unit=sfrtab_unit,file=filename,form='unformatted',status='unknown')
    write(sfrtab_unit) tsno(st)%n_total,nrej,nfile  ! nb of galaxies, nb of agebins, and nb of metallicities
    write(sfrtab_unit) timetab(1:nrej)
    write(sfrtab_unit) tabmetspec(1:nfile)
    do j=1,tsno(st)%nb_of_halos
       do k=1,tsno(st)%liste_halos(j)%datas%nbgal
          sfrtot = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%totsfr
          if ( sfrtot > 0.0d0) then 
             n1 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr_tab,dim=1)
             n2 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr_tab,dim=2)
             write(sfrtab_unit) sfrtot, n1, n2
             write(sfrtab_unit) tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr_tab(:,:)
          else 
             n1 = 0 
             n2 = 0 
             write(sfrtab_unit) sfrtot, n1, n2
          end if
       end do
    end do
    close(sfrtab_unit)
#ifdef DUAL_IMF 
    write(filename,'(a,a,i3.3)') trim(data_dir),'/bulge_sfrtabs_imf2.',st
    open(unit=sfrtab_unit,file=filename,form='unformatted',status='unknown')
    write(sfrtab_unit) tsno(st)%n_total,nrej,nfile  ! nb of galaxies, nb of agebins, and nb of metallicities
    write(sfrtab_unit) timetab(1:nrej)
    write(sfrtab_unit) tabmetspec(1:nfile)
    do j=1,tsno(st)%nb_of_halos
       do k=1,tsno(st)%liste_halos(j)%datas%nbgal
          sfrtot = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%totsfr2
          if (sfrtot > 0.0d0) then 
             n1 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr_tab2,dim=1)
             n2 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr_tab2,dim=2)
             write(sfrtab_unit) sfrtot, n1, n2
             write(sfrtab_unit) tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%bulge%sfr_tab2(:,:)          
          else 
             n1 = 0
             n2 = 0 
             write(sfrtab_unit) sfrtot, n1, n2
          end if
       end do
    end do
    close(sfrtab_unit)
#endif

    ! write out burst component 
    write(filename,'(a,a,i3.3)') trim(data_dir),'/burst_sfrtabs.',st
    open(unit=sfrtab_unit,file=filename,form='unformatted',status='unknown')
    write(sfrtab_unit) tsno(st)%n_total,nrej,nfile  ! nb of galaxies, nb of agebins, and nb of metallicities
    write(sfrtab_unit) timetab(1:nrej)
    write(sfrtab_unit) tabmetspec(1:nfile)
    do j=1,tsno(st)%nb_of_halos
       do k=1,tsno(st)%liste_halos(j)%datas%nbgal
          sfrtot = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%totsfr
          if (sfrtot > 0.0d0) then 
             n1 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr_tab,dim=1)
             n2 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr_tab,dim=2)
             write(sfrtab_unit) sfrtot, n1, n2 
             write(sfrtab_unit) tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr_tab(:,:)
          else
             n1 = 0
             n2 = 0 
             write(sfrtab_unit) sfrtot, n1, n2 
          end if
       end do
    end do
    close(sfrtab_unit)
#ifdef DUAL_IMF 
    write(filename,'(a,a,i3.3)') trim(data_dir),'/burst_sfrtabs_imf2.',st
    open(unit=sfrtab_unit,file=filename,form='unformatted',status='unknown')
    write(sfrtab_unit) tsno(st)%n_total,nrej,nfile  ! nb of galaxies, nb of agebins, and nb of metallicities
    write(sfrtab_unit) timetab(1:nrej)      ! age bins 
    write(sfrtab_unit) tabmetspec(1:nfile)  ! metallicity bins 
    do j=1,tsno(st)%nb_of_halos
       do k=1,tsno(st)%liste_halos(j)%datas%nbgal
          sfrtot = tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%totsfr2
          if (sfrtot > 0.0d0) then 
             n1 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr_tab2,dim=1)
             n2 = size(tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr_tab2,dim=2)
             write(sfrtab_unit) sfrtot, n1, n2
             write(sfrtab_unit) tsno(st)%liste_halos(j)%datas%liste_galaxies(k)%burst%sfr_tab2(:,:)          
          else 
             n1 = 0
             n2 = 0
             write(sfrtab_unit) sfrtot, n1, n2
          end if
       end do
    end do
    close(sfrtab_unit)
#endif

    return

  end subroutine output_sfr_tabs
#endif
  !*****************************************************************************************************************

  subroutine output_acc_stories(st)

    implicit none

    character(MAXPATHSIZE) :: filename
    integer(kind=4)        :: ih,ig,st
   
!!$    write(filename,'(a,a)') trim(data_dir),'/acc_story.dat'
!!$    open(unit=acc_unit,file=filename,form='unformatted',status='unknown')
!!$    write(acc_unit) MaxAccretionEvents
!!$    write(acc_unit) (((ig+0.5)*accretion_deltat),ig=1,MaxAccretionEvents)
!!$    do ih = 1,tsno(st)%nb_of_halos
!!$       do ig=1,tsno(st)%liste_halos(ih)%datas%nbgal
!!$          write(acc_unit) tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%acc%accretion_rate
!!$       end do
!!$    end do
!!$    close(acc_unit)
!!$
!!$    write(filename,'(a,a)') trim(data_dir),'/wind_story.dat'
!!$    open(unit=acc_unit,file=filename,form='unformatted',status='unknown')
!!$    write(acc_unit) MaxAccretionEvents
!!$    write(acc_unit) (((ig+0.5)*accretion_deltat),ig=1,MaxAccretionEvents)
!!$    do ih = 1,tsno(st)%nb_of_halos
!!$       do ig=1,tsno(st)%liste_halos(ih)%datas%nbgal
!!$          write(acc_unit) tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)%wind%accretion_rate
!!$       end do
!!$    end do
!!$    close(acc_unit)
    
    return

  end subroutine output_acc_stories

  !*****************************************************************************************************************
#ifdef DEBUG_OUTPUTS
  subroutine debug_outputs(st)

    implicit none
    
    integer(kind=4)        :: st,ih,ig,ng,i,i1,i10,i5
    real(kind=8)           :: mcum,t, mstar,mism,mwind,macc,mhot
    character(MAXPATHSIZE) :: filename
    type(halo),pointer     :: h
    type(galaxy),pointer   :: g

    ! count the total number of central galaxies 
    ng = 0
    do ih=1,tsno(st)%nb_of_halos
       if (tsno(st)%liste_halos(ih)%datas%nbgal > 0) ng = ng + 1
    end do
    ! dump their properties and those of associated haloes

!tibo
#ifndef MIN_IO
    write(filename,'(a,a,i3.3,a)') trim(data_dir),'/centraldiscprops',st,'.dat'
    open(unit=debug_unit,file=filename,status='unknown',form='formatted')
    write(debug_unit,'(i8)') ng
    do ih=1,tsno(st)%nb_of_halos
       if (tsno(st)%liste_halos(ih)%datas%nbgal > 0) then 
          h => tsno(st)%liste_halos(ih)
          call get_central_galaxy(h,ig)
          g => tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)
          write(debug_unit,'(11(e14.6,1x),i8,1x,3(e14.6,1x))') g%disc%mgal,g%disc%rgal,g%disc%mcold,g%disc%tdyn,g%disc%minstar,g%disc%mcoldz,g%disc%speed, &
               & g%disc%sfr100, h%datas%rvir,h%datas%mvir,h%datas%cvel,h%datas%nbgal,h%rfof,h%mfof,h%spin
       end if
    end do
    close(debug_unit)

    ! props of central galaxies 
    write(filename,'(a,a,i3.3,a)') trim(data_dir),'/centralprops',st,'.dat'
    open(unit=debug_unit,file=filename,status='unknown',form='formatted')
    write(debug_unit,'(i8)') ng
    do ih=1,tsno(st)%nb_of_halos
       if (tsno(st)%liste_halos(ih)%datas%nbgal > 0) then 
          h => tsno(st)%liste_halos(ih)
          call get_central_galaxy(h,ig)
          g => tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)
          mstar = max(tsno(st)%age_univm100, tsno(st-1)%age_univ)
          call accretion_mass(g%acc,mstar,tsno(st)%age_univ,macc,mism)
          write(debug_unit,'(10(e14.6,1x))') h%datas%rvir,h%datas%mvir,h%datas%cvel,h%rfof,h%mfof, & 
               & g%disc%sfr100+g%bulge%sfr100+g%burst%sfr100,g%disc%sfr10+g%bulge%sfr10+g%burst%sfr10, total_stellar_mass(g),total_gal_gas_mass(g), & 
               & macc
       end if
    end do
    close(debug_unit)

    ! dump baryonic decompositition for each halo
    write(filename,'(a,a,i3.3,a)') trim(data_dir),'/halo_baryondec',st,'.dat'
    open(unit=debug_unit,file=filename,status='unknown',form='formatted')
    do ih = 1,tsno(st)%nb_of_halos
       if (tsno(st)%liste_halos(ih)%datas%nbgal > 0) then 
          h => tsno(st)%liste_halos(ih)
          mstar = 0.0d0
          macc  = 0.0d0
          mwind = 0.0d0
          mhot  = 0.0d0
          mism  = 0.0d0
          do ig = 1,tsno(st)%liste_halos(ih)%datas%nbgal
             g => tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)
             mstar = mstar + total_stellar_mass(g)
             mism  = mism  + total_gal_gas_mass(g)
             mwind = mwind + accretion_remaining_mass(g%wind,g%tgal)
             macc  = macc  + accretion_remaining_mass(g%acc,g%tgal)
          end do
          mhot = h%datas%mhotgaz
          write(debug_unit,'(i8,1x,9(e14.6,1x))') h%datas%nbgal,cold_gas_in_halo(h),mhot,h%datas%mcoldgaz,mstar,mism,mwind,macc,h%mfof,h%datas%mvir
       end if
    end do
    close(debug_unit)

#endif

!!$    ! output times of last accretion events (from filaments, not wind), for central galaxies
!!$    write(filename,'(a,a,i3.3,a)') trim(data_dir),'/lastaccretiontimes',st,'.dat'
!!$    open(unit=debug_unit,file=filename,status='unknown',form='formatted')
!!$    write(debug_unit,'(i8)') ng
!!$    do ih=1,tsno(st)%nb_of_halos
!!$       if (tsno(st)%liste_halos(ih)%datas%nbgal > 0) then 
!!$          h => tsno(st)%liste_halos(ih)
!!$          call get_central_galaxy(h,ig)
!!$          g => tsno(st)%liste_halos(ih)%datas%liste_galaxies(ig)
!!$          i1   = 0
!!$          i5   = 0
!!$          i10  = 0
!!$          mcum = 0
!!$          do i=MaxAccretionEvents,1,-1
!!$             t = tsno(st)%age_univ
!!$             if (i * accretion_deltat > t) cycle
!!$             if (g%acc%accretion_rate(i) > 0.0d0) then 
!!$                if (i1 == 0) i1 = i
!!$                mcum = mcum + g%acc%accretion_rate(i) * accretion_deltat
!!$                if (mcum > 0.01 * omega_b / omega_0  *h%mfof .and. i5 == 0) then 
!!$                   i5 = i
!!$                end if
!!$                if (mcum > 0.05 * omega_b / omega_0  *h%mfof) then 
!!$                   i10 = i
!!$                   exit
!!$                end if
!!$             end if
!!$          end do
!!$          write(debug_unit,'(11(e14.6,1x))') h%mfof,i1*accretion_deltat,i5*accretion_deltat,i10*accretion_deltat,g%disc%mgal,g%disc%mcold,g%disc%minstar,g%disc%rgal,g%disc%speed,h%datas%mvir,h%datas%rvir,h%datas%cvel
!!$       end if
!!$    end do
!!$    close(debug_unit)


    

    return

  end subroutine debug_outputs

  !*****************************************************************************************************************
  
  subroutine output_gbk(ts)

    implicit none
    
    integer(kind=4) :: ts,unit,ig,ih,j
    real(kind=8)    :: dum,m
    type(galaxy),pointer    :: g
    character(MAXPATHSIZE) :: filename
    

    if (ts < 2) return
    unit = 12
!tibo
!#ifndef MIN_IO
    write(filename,'(a,a,i3.3)') trim(data_dir),'/gbk.',ts
    open(unit=unit,file=filename,status='unknown',form='formatted')
    write(unit,'(i3,1x,i6,1x,e14.6,1x,e14.6)') ts,tsno(ts)%n_total,tsno(ts)%aexp,tsno(ts)%age_univ
    do ih = 1,tsno(ts)%nb_of_halos
       if (tsno(ts)%liste_halos(ih)%datas%nbgal > 0) then 
          !tibo: output halo accretion for each galaxy (0 if satellite, all the halo accretion if central)
          do ig = 1,tsno(ts)%liste_halos(ih)%datas%nbgal
             !call get_central_galaxy(tsno(ts)%liste_halos(ih),ig)
             g => tsno(ts)%liste_halos(ih)%datas%liste_galaxies(ig)
             !if (ig .ne. j) then
             !   write(unit,'(13(e14.6,1x))') 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0
             !else
             call accretion_mass(g%acc,tsno(ts-1)%age_univ,tsno(ts)%age_univ,m,dum)
             call inc_gbk_macc_stream(g%gbk,m)
             call accretion_mass(g%wind,tsno(ts-1)%age_univ,tsno(ts)%age_univ,m,dum)
             call inc_gbk_macc_fountain(g%gbk,m)
             write(unit,'(13(e14.6,1x))') tsno(ts)%liste_halos(ih)%macc*1d11,tsno(ts)%liste_halos(ih)%mfof*1d11,&
                     cold_gas_in_halo_streams(tsno(ts)%liste_halos(ih))*1d11, &
                     tsno(ts)%liste_halos(ih)%datas%mhotgaz*1d11,tsno(ts)%liste_halos(ih)%datas%mhotz*1d11,&
                     g%gbk%macc_stream*1d11,g%gbk%macc_fountain*1d11,g%gbk%mejcold*1d11,g%gbk%mejhot*1d11,g%gbk%mejout*1d11,&
                     g%gbk%mstarform*1d11,g%gbk%mstarform_bulge*1d11,tsno(ts)%age_univ-tsno(ts-1)%age_univ
             
             !end if
          end do
       end if
    end do
    
    close(unit)
!#endif    

    return

  end subroutine output_gbk

  !*****************************************************************************************************************


#endif
  !*****************************************************************************************************************
  
  subroutine output_halo_accinfo

    implicit none
    
    character(MAXPATHSIZE) :: filename
    integer(kind=4)        :: ts,ih,unit
    real(kind=8)           :: delta_t 

!tibo
#ifndef MIN_IO    
    unit = 12
    do ts = 2,nsteps_do
       delta_t = tsno(ts)%age_univ - tsno(ts-1)%age_univ ! in Gyr
       
       write(filename,'(a,a,i3.3)') trim(data_dir),'/HaloAcc.',ts
       open(unit=unit,file=filename,status='unknown',form='formatted')
       write(unit,*) tsno(ts)%nb_of_halos
       do ih = 1,tsno(ts)%nb_of_halos
          write(unit,'(2(e14.6,1x))') tsno(ts)%liste_halos(ih)%macc/delta_t*100.d0,tsno(ts)%liste_halos(ih)%mfof
       end do
       close(unit)
    end do
#endif

    return

  end subroutine output_halo_accinfo


  !***************************************************************************************************************** 
  subroutine output_halo_props
    
    implicit none
    
    character(MAXPATHSIZE) :: filename
    integer(kind=4)        :: ts,ih,unit
    real(kind=8)           :: delta_t 
    
    unit = 61
    do ts = 2,nsteps_do
     !  if ((ts .lt. 23) .or. (ts .ge. 93)) then
       delta_t = tsno(ts)%age_univ - tsno(ts-1)%age_univ ! in Gyr                                   
       write(filename,'(a,a,i3.3)') trim(data_dir),'/Halo_props.',ts
       open(unit=unit,file=filename,status='unknown',form='formatted')
       write(unit,*) tsno(ts)%nb_of_halos
       do ih = 1,tsno(ts)%nb_of_halos
          !if (ih .eq. 1) then
             !write(errunit,*) 'Mvir, rfof, spin, hid', tsno(ts)%liste_halos(ih)%datas%mvir,tsno(ts)%liste_halos(ih)%rfof,tsno(ts)%liste_halos(ih)%spin,tsno(ts)%liste_halos(ih)%HaloID
          !end if
          write(unit,'(6(e14.6,1x),i20,1x)') tsno(ts)%liste_halos(ih)%macc/delta_t*100.d0,tsno(ts)%liste_halos(ih)%mfof,tsno(ts)%liste_halos(ih)%datas%mvir,tsno(ts)%liste_halos(ih)%rfof,tsno(ts)%liste_halos(ih)%datas%rvir,tsno(ts)%liste_halos(ih)%spin,tsno(ts)%liste_halos(ih)%HaloID
       end do
       close(unit)
       !end if
    end do
    
    return

  end subroutine output_halo_props

  !*****************************************************************************************************************
  !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
end module IO_BARYONS

