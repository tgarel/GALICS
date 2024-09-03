module cooling
  
  public

  !*****************************************************************************************************************
  ! variables pertaining to cooling files
  !*****************************************************************************************************************
  integer(kind=4)              :: ncool                            ! number of temperatures for cooling curves
  integer(kind=4)              :: nmets                            ! number of metallicities for cooling curves
  real(kind=4),allocatable     :: templ(:)                         ! values of temperature for cooling curves
  real(kind=4),allocatable     :: coolcurlmet(:,:)                 ! cooling curves as a function of metallicity
  real(kind=4),allocatable     :: coolcurl(:)                      ! cooling curve
  real(kind=4),allocatable     :: tabmetcool(:)                    ! metallicity array for cooling curves
  !*****************************************************************************************************************

contains

  !*****************************************************************************************************************
  subroutine read_cooling_curves
    
    ! read the 'cooling_curves_new.dat' file which contains metallicity-dependent cooling rates
    ! from Sutherland & Dopita (1992)
    
    implicit none              
    
    integer(kind=4) :: ierr,i  
    character(200)  :: file
    
    write(file,'(a,a)') trim(inputpath),'/cooling_curves_new.dat'
    call test_stop(file)
    open(cool_unit,file=file,status='old')
    read(cool_unit,*)                                           ! connect on 1st line.   
    read(cool_unit,*) ncool,nmets
    allocate(templ(ncool),coolcurlmet(nmets,ncool),coolcurl(ncool),tabmetcool(nmets),stat=ierr)
    if (ierr /= 0 ) stop '> error in opening cooling curve file'
    templ       = 0.0
    coolcurlmet = 0.0
    coolcurl    = 0.0
    tabmetcool  = 0.0
    read(cool_unit,*) tabmetcool(1:nmets)                       ! log of the metallicity
    do i=1,ncool
       read (cool_unit,*) templ(i),coolcurlmet(1:nmets,i)
    end do
    coolcurlmet = coolcurlmet+23.0  ! convert to 10^-23 erg cm^3/s  
    close(cool_unit)
    
    return
    
  end subroutine read_cooling_curves
  
  !*****************************************************************************************************************
  
  subroutine compute_cooling(rvir,tvir,mhot,mhotz,t1,t2,mcool,mcoolz)
    
    ! compute mass of gas (and metals) which cools from t1 to t2
    ! in a halo defined by rvir, mhot etc... 

    lambda = get_lambda(tvir,mhotz,mhot) ! cooling rate from SD93 in units of 10^-23 erg.cm^3/s
    mcool = !some function of lambda... 
    
    if (get_mcool > mhotgas) then 
       mcool = mhotgas
    end if
    
  end subroutine compute_cooling
    
  !*****************************************************************************************************************


end module cooling
