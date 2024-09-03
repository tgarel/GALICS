module GLOB_DEFS

  use accretion_type
  use gal_bookkeep_type

  public

  !-------------------------------------------------------------------------------------------------------------  
  ! input-output units
  !-------------------------------------------------------------------------------------------------------------  
  integer(kind=4),parameter :: errunit        = 6            ! output unit for messages on screen (6 is stdio)
  integer(kind=4),parameter :: logerrunit     = 133          ! output unit for error messages in file fort.logerrunit 
  integer(kind=4),parameter :: baryons_unit   = 18
  integer(kind=4),parameter :: filters_unit   = 19
  integer(kind=4),parameter :: param_unit     = 20
  integer(kind=4),parameter :: spec_unit      = 21
  integer(kind=4),parameter :: tree_unit      = 22
  integer(kind=4),parameter :: hg_cnt_unit    = 23
  integer(kind=4),parameter :: gcnt_unit      = 24
  integer(kind=4),parameter :: sfr_unit       = 25
  integer(kind=4),parameter :: merging_unit   = 26
  integer(kind=4),parameter :: bookkeep_unit  = 27
  integer(kind=4),parameter :: tau_unit       = 28
  integer(kind=4),parameter :: tigm_unit      = 29
  integer(kind=4),parameter :: ce_unit        = 30
  integer(kind=4),parameter :: galres_unit    = 31
  integer(kind=4),parameter :: cores_unit     = 32
  integer(kind=4),parameter :: galsfr_unit    = 33
  integer(kind=4),parameter :: numbers_unit   = 34
  integer(kind=4),parameter :: hgres_unit     = 35
  integer(kind=4),parameter :: hgprof_unit    = 36
  integer(kind=4),parameter :: hreserv_unit   = 37
  integer(kind=4),parameter :: hdmres_unit    = 38
#ifndef BIG_RUN 
  integer(kind=4),parameter :: contam_unit    = 39
#endif
  integer(kind=4),parameter :: wave_unit      = 40
  integer(kind=4),parameter :: gal_tree_unit  = 41
  integer(kind=4),parameter :: asc_tree_unit  = 42
  integer(kind=4),parameter :: ids_unit       = 43
#ifdef RECORD_SFR
  integer(kind=4),parameter :: sfrtab_unit    = 44
#endif
  integer(kind=4),parameter :: acc_unit       = 45
#ifdef DEBUG_OUTPUTS
  integer(kind=4),parameter :: debug_unit     = 46
#endif
  integer(kind=4),parameter :: eject_unit     = 47
  integer(kind=4),parameter :: mfs_unit       = 48
  integer(kind=4),parameter :: burst_out_unit = 49
  integer(kind=4),parameter :: wind_unit      = 50
  integer(kind=4),parameter :: shell_unit     = 51
  integer(kind=4),parameter :: uv1200_dblg_unit = 52
  integer(kind=4),parameter :: uv1200_gal_unit  = 53
  integer(kind=4),parameter :: uv1200_dblg_ext_unit = 54
  integer(kind=4),parameter :: uv1200_burst_unit = 55
  integer(kind=4),parameter :: uv1500_gal_unit = 56
  integer(kind=4),parameter :: ej_burst_unit     = 57
  integer(kind=4),parameter :: age_lum_weighted_unit = 58
  integer(kind=4),parameter :: beta_unit = 59
#ifdef DESIRED_OUTPUTS
  integer(kind=4),parameter :: desired_unit = 60
#endif
  ! TIBO
#ifdef LYA  
  integer(kind=4),parameter :: lyaout_unit = 62
#endif
  ! OBIT
  
  !-------------------------------------------------------------------------------------------------------------  
  ! useful data types to deal with halos and galaxies
  !-------------------------------------------------------------------------------------------------------------  
  type vector                      ! vector coordinates
     real(kind=8)                  :: x
     real(kind=8)                  :: y
     real(kind=8)                  :: z
  end type vector

  type father            ! list of progenitors of a halo
     integer(kind=4)               :: nb_fathers         ! total number of progenitors
     integer(kind=4),pointer       :: list_fathers(:)    ! list of their indices
     real(kind=8),pointer          :: mass_fathers(:)    ! array containing % of father's mass given to son
  end type father

  type tree_info         ! simpler tree information for halos : does not include the background as a progenitor
     integer(kind=4)               :: frag               ! does halo originates from previous fragmentation (0=no, 1=yes)  
     integer(kind=4)               :: ndads              ! total number of progenitors : 0 (background), 1..n
     integer(kind=4),pointer       :: dads(:)            ! list of their indices
     integer(kind=4)               :: nsons              ! total number of sons        : 0 (end of branch) or 1
     integer(kind=4),pointer       :: sons(:)            ! list of their indices
  end type tree_info

  type shape             ! shape (triaxiality) of a halo
     real(kind=8)                  :: a                  ! 1st principal axis size
     real(kind=8)                  :: b                  ! 2nd    "       "    "
     real(kind=8)                  :: c                  ! 3rd    "       "    "
  end type shape

  type gal_comp           ! component of a galaxy (i.e. disc, bulge or burst; NB: burst == small bulge) 

     real(kind=8)                  :: r_output           
     real(kind=8)                  :: mgal               ! total mass (stars and cold gas: minstar+mcold)
     real(kind=8)                  :: rgal               ! exp scale (disc), Hernquist's r_B (bulge)
     real(kind=8)                  :: mcold              ! mass of cold gas (includes metals)
     real(kind=8)                  :: tdyn               ! dyn time (.5 turn @ r1/2 [disc], r1/2 --> 0 [bulge])
     real(kind=8)                  :: minstar            ! total mass in stars
     real(kind=8)                  :: rstrip             ! RAM pressure stripping radius for gas
     real(kind=8)                  :: sfr1               ! instantaneous star formation rate (over 1Myr)
     real(kind=8)                  :: sfr10              ! sfr averaged over 10Myr
     real(kind=8)                  :: sfr100             ! sfr averaged over 100Myr
     real(kind=8)                  :: mcoldz             ! mass of metals in cold phase
     real(kind=8)                  :: speed              ! rot. speed (disc), 3D vel disp (bulge)  
     real(kind=8)                  :: transp             ! mass of stars transported to other components
     real(kind=8)                  :: Nion_phot          ! Number of ionizing photons emitted per second
     real(kind=8)                  :: contF_1216         ! UV flux summed on 201 bins (between 1105 and 1305 A) around lya

#ifdef DUAL_IMF
     real(kind=8)                  :: minstar2           ! total mass of stars of second IMF (==sum(sfh_tab2))
     real(kind=8)                  :: sfr21              ! ...
     real(kind=8)                  :: sfr210             ! sfr averaged over 10 Myr
     real(kind=8)                  :: sfr2100            ! sfr averaged over 100 Myr
     real(kind=8)                  :: transp2            ! ... 
#endif
     real(kind=8),pointer          :: sfh_tab(:,:)       ! star formation history (dim1 = age, dim2 = mets)
#ifdef DUAL_IMF 
     real(kind=8),pointer          :: sfh_tab2(:,:)      ! star formation history for second stellar population
#endif
#ifdef RECORD_SFR
     real(kind=8)                  :: totsfr             ! total mass of stars formed (== sum(sfr_tab))
     real(kind=8),pointer          :: sfr_tab(:,:)       ! real star formation history (i.e. mass of stars actually
#ifdef DUAL_IMF 
     real(kind=8),pointer          :: sfr_tab2(:,:)      ! formed, not remaining after evolution...)
     real(kind=8)                  :: totsfr2
#endif
#endif
  end type gal_comp

  type quasar            ! quasar
     real(kind=8)                  :: mass               ! mass of QSO , 10^11 M_sun
     real(kind=8)                  :: accretion_rate     ! accretion rate, 10^11 M_sun / Gyr
  end type quasar

  type DM_core           ! DM core attached to a satellite galaxy used to compute dynamical friction more accurately
     real(kind=8)                  :: mass               ! mass of core
     real(kind=8)                  :: rho_0              ! central density
     real(kind=8)                  :: r_0                ! core radius  
     real(kind=8)                  :: fdm                ! fraction of dark matter still in the core 
  end type DM_core

  type gal_prog
     integer(kind=4)               :: nb_prog            ! galaxy's number of progs during a timestep
     integer(kind=4),pointer       :: hno(:)             ! list of halo numbers of progs of gal
     integer(kind=4),pointer       :: gno(:)             ! list of gal numbers of progs of gal
  end type gal_prog

  type gal_girls
     integer(kind=4)               :: hno                ! halo number of gal's descendent
     integer(kind=4)               :: gno                ! gal number of gal's descendent
  end type gal_girls

  type shell
     real(kind=8)                  :: radius_subm1
     real(kind=8)                  :: radius
     real(kind=8)                  :: speed
     real(kind=8)                  :: age
     real(kind=8)                  :: zmet
     real(kind=8)                  :: nh
     real(kind=8)                  :: tau_dust
     real(kind=8)                  :: Rfrag
     real(kind=8)                  :: Rmax
     real(kind=8)                  :: Mshell
     real(kind=8)                  :: Mshellz
     logical(kind=4)               :: onoff
     real(kind=8)                  :: tstart
     real(kind=8)                  :: c_factor
     real(kind=8)                  :: mcold_gal
     real(kind=8)                  :: rgal
     integer(kind=4)               :: ii
     integer(kind=4)               :: oo
     integer(kind=4)               :: a
     real(kind=8)                  :: v_esc_halo
     integer(kind=4)               :: vv
     integer(kind=4)               :: c
     real(kind=8)                  :: Rfinal

  end type shell
  
  type igm
     real(kind=8)                  :: mcold
     real(kind=8)                  :: mcoldz
  end type igm

  type galaxy            ! contains all properties necessary to define a galaxy
     integer(kind=4)               :: my_number          ! index of galaxy in list
     integer(kind=4)               :: hno                ! number of parent halo
     type(gal_prog)                :: my_progs           ! progs IDs of galaxy
     type(gal_girls)               :: my_girls           ! daughter ID of galaxy
     integer(kind=4)               :: nb_merg            ! number of mergers it has undergone previously
#ifndef BIG_RUN
     integer(kind=4)               :: contam             ! flag contaminated galaxies (resimulations only)
#endif
     integer(kind=4)               :: inst_bulg          ! flag galaxy disc which were unstable during timestep 
     real(kind=8)                  :: tgal               ! galaxy clock time at timestep st (st,st-1,midstep)
     real(kind=8)                  :: tbirth             ! time at which galaxy was first spotted   
     real(kind=8)                  :: tmerge             ! expected merging time.
     real(kind=8)                  :: disturb            ! disturbance of morphological type (recent merger)
     real(kind=8)                  :: r                  ! galaxy's orbital radius in halo
     real(kind=8)                  :: inclination        ! inclination of disc
     type(vector)                  :: p                  ! galaxy position 
     type(vector)                  :: v                  ! galaxy velocity through the Universe !
     type (DM_core)                :: core               ! galaxy's attached DM core 
     type (gal_comp)               :: disc               ! galaxy's disc component
     type (gal_comp)               :: bulge              ! galaxy's bulge    "
     type (gal_comp)               :: burst              ! galaxy's burst    "
     type (quasar)                 :: QSO                ! galaxy's QSO
     type (accretion)              :: acc                ! accretion (of cold gas) props of a galaxy.
     type (accretion)              :: wind               ! ejected gas which falls back onto a galaxy
     type (gal_bookkeep)           :: gbk                ! bookkeeping of quantities integrated over timestep
     real(kind=8),dimension(21)    :: uv_1200_gal        ! UV spectrum from 1105 to 1305A (in 10^11Msun erg/s/µm)
     real(kind=8),dimension(21)    :: uv_1200_dblg       ! UV spectrum from 1105 to 1305A (in 10^11Msun erg/s/µm)
     real(kind=8),dimension(21)    :: uv_1200_burst      ! UV spectrum from 1405 to 1605A (in 10^11Msun erg/s/µm)
     real(kind=8),dimension(21)    :: uv_1500_burst      ! UV spectrum from 1405 to 1605A (in 10^11Msun erg/s/µm)
     type (accretion)              :: burst_outflow      ! mass ejected by burst
     type (shell)                  :: shell
     real(kind=8)                  :: E_int_sn           ! integrated mechanical energy injected by SNe
     real(kind=8)                  :: wind1
     real(kind=8)                  :: wind2
     real(kind=8)                  :: wind10
     real(kind=8)                  :: wind20
     real(kind=8)                  :: wind50
     real(kind=8)                  :: wind100
     real(kind=8),dimension(200)   :: galej
     real(kind=8),dimension(200)   :: burst_ej
     integer(kind=4)               :: nsubsteps
     real(kind=8)                  :: age_lum_weighted   ! age weighted by L_1200

  end type galaxy

  type baryon            ! baryonic content of a halo
     real(kind=8)                  :: mgaz               ! total mass of gas (hot + cold) in the halo
     real(kind=8)                  :: mhotgaz            ! mass of hot gas (includes metals)
     real(kind=8)                  :: mcoldgaz           ! mass of cold gas (includes metals)
     real(kind=8)                  :: mhotz              ! mass of metals in hot phase
     real(kind=8)                  :: mcoldz             ! mass of metals in cold phase
     real(kind=8)                  :: mgazout            ! mass of gas ejected from the halo (cumulative in time) -> igm 
     real(kind=8)                  :: metalsout          ! mass of metas ejected from the halo -> igm 
     real(kind=8)                  :: rvir               ! virial radius (DM)
     real(kind=8)                  :: mvir               ! virial mass (DM)
     real(kind=8)                  :: tvir               ! virial temperature (DM)
     real(kind=8)                  :: cvel               ! circular velocity
     integer(kind=4)               :: nbgal              ! total number of galaxies sitting in halo
     type (galaxy),pointer         :: liste_galaxies(:)  ! list of these galaxies
  end type baryon

  type son               ! son of a halo
     integer(kind=4)               :: nb_sons            ! total number of sons       
     integer(kind=4),pointer       :: list_sons(:)       ! list of indices
  end type son

  type hprofile          ! parameters defining halo density profile (DM)
     real(kind=8)                  :: rho_0              ! central density (virial density for TSIS)
     real(kind=8)                  :: r_c                ! core radius (virial radius for TSIS)
  end type hprofile

  type halo              ! DM halo
     integer(kind=4)               :: my_number          ! number of halo in list
     integer(kind=4)               :: my_timestep        ! number of timestep at which halo is present
     real(kind=8)                  :: my_form_aexp       ! exp factor when halo got half its mass
     integer(kind=4)               :: ncont              ! contamination parameter
     type (tree_info)              :: tree               ! links to progenitors & sons via simple tree
     type (baryon)                 :: datas              ! baryonic content (including galaxies)
     type (father)                 :: my_fathers         ! list of progenitors
     type (son)                    :: my_sons            ! list of sons
     type (shape)                  :: sh                 ! triaxiality of halo information
     type (vector)                 :: p                  ! position vector of center of mass in box
     type (vector)                 :: v                  ! proper velocity vector of center of mass
     type (vector)                 :: L                  ! angular momentum vector
     type (hprofile)               :: halo_profile       ! density profile parameters of halo
     real(kind=8)                  :: mfof               ! fof mass of halo (DM only)
     real(kind=8)                  :: macc               ! mass accreted by a halo from diffuse particles
     real(kind=8)                  :: rfof               ! distance of most distant part. to center of halo
     real(kind=8)                  :: spin               ! spin parameter of halo
     real(kind=8)                  :: ek                 ! kinetic energy
     real(kind=8)                  :: ep                 ! potential energy
     real(kind=8)                  :: et                 ! total energy
     type(igm)                     :: igm                ! igm content in the neighbourhood of the halo
#ifdef DEFINE_IDS 
     integer(kind=4)               :: BushID             ! ID of the independant tree the halo is in (read from treemaker)
     integer(kind=8)               :: HaloID             ! Unique ID of a halo in the whole simulation (defined here)
#endif
  end type halo

  type halo_ptr          ! points on a single halo
     type (halo),pointer           :: p 
  end type halo_ptr

  ! it is necessary to define an 'average' column density for the models.  In practice, eg for an expon disc, the
  ! average density is zero (there is a tail out to infinity) --> compute luminosity-weighted average density.   
  ! obsc_rad is the ratio of scale-length to radius that should be used for computing the average column density.
  type comp_info    ! info about galaxy components (disc, bulge, burst) 
     real(kind=8)                  :: s_to_m             ! ratio of scale length to median mass length
     real(kind=8)                  :: esc_para           ! ratio of escape to characteristic velocity
     real(kind=8)                  :: obsc_rad           ! radius for mass-weighted average column density
     real(kind=8)                  :: typ_vdisp          ! typical velocity dispersion (turb. pres.) for cold gas
  end type comp_info

  type tshalo            ! 'timestep' structure: list of halos, & imp variables that would need separate allocation 
     type (halo),pointer           :: liste_halos(:)     ! list of all halos at present timestep  
     real(kind=8)                  :: mass_bary          ! baryonic mass in the simulation box
     real(kind=8)                  :: mass_halo          ! halo      "    "  "     "         "
     real(kind=8)                  :: mass_cold          ! cold gas  "    "  "     "         "   (includes stars, in fact)
     real(kind=8)                  :: mass_star          ! stellar   "    "  "     "         "
#ifdef DUAL_IMF
     real(kind=8)                  :: mass_star2         ! ... with second IMF
#endif
     real(kind=8)                  :: mass_metals        ! metals    "    "  "     "         "
     real(kind=8)                  :: mass_gal_gas       ! mass of gas in galaxies summed over whole simulation volume
     real(kind=8)                  :: mass_hot_gas       ! total mass of hot gas in all halos
     real(kind=8)                  :: mass_hot_mets      ! total mass of metals in hot gas of halos
     real(kind=8)                  :: average_tau        ! average merging time between satellite galaxies   
     integer(kind=4)               :: n_total            ! total number of galaxies
     integer(kind=4)               :: n_ss               ! total number of satellite-satellite galaxy mergers
     integer(kind=4)               :: n_dynfric          ! total number of dynamical friction mergers
     integer(kind=4)               :: n_form             ! total number of newly formed galaxies
     integer(kind=4)               :: n_lost             ! total number of galaxies lost bc host halo has no descendent
     integer(kind=4)               :: n_biglost          ! number of galaxies lost in halos with mass > 1e12 M_sun
     integer(kind=4)               :: n_unst             ! number of galaxies which disk was unstable during timestep
     integer(kind=4)               :: tau_count          ! number of satellite-satellite mergers in a halo
     integer(kind=4)               :: nb_of_halos        ! total number of halos 
     integer(kind=4)               :: nb_maj_mergs       ! total number of major mergers (mass ratio > xi)
     integer(kind=4)               :: nb_min_mergs       ! total number of minor mergers (mass ratio < xi)
     integer(kind=4)               :: nb_halo_mergs      ! total number of halo mergers
     real(kind=8)                  :: aexp               ! expansion factors 
     real(kind=8)                  :: age_univ           ! age of the universe
     real(kind=8)                  :: age_univm2         ! age of the universe minus 2Myr
     real(kind=8)                  :: age_univm10        ! age of the universe minus 10Myr
     real(kind=8)                  :: age_univm20        ! age of the universe minus 20Myr
     real(kind=8)                  :: age_univm50        ! age of the universe minus 50Myr
     real(kind=8)                  :: age_univm100       ! age of the universe minus 100Myr     
     real(kind=8)                  :: omega_t            ! omega_matter(t)
     real(kind=8)                  :: global_SF          ! global star formation rate in the simulation box
#ifdef DUAL_IMF
     real(kind=8)                  :: global_SF2         ! global SFR with second IMF 
#endif
     real(kind=8)                  :: t_igm              ! temperature of the IGM 
#ifdef RENUMBERING 
     integer(kind=4), allocatable  :: list_halos_number(:) ! data duplication useful when calling renumbering
#endif
  end type tshalo
     
  !-------------------------------------------------------------------------------------------------------------  
  ! I/O global variables
  !-------------------------------------------------------------------------------------------------------------  
#ifndef MAXPATHSIZE
  integer(kind=4):: MAXPATHSIZE
  parameter(MAXPATHSIZE = 512)
#endif
  character(MAXPATHSIZE)       :: data_dir                        ! directory for param files and data storage
  character(MAXPATHSIZE)       :: treefile                        ! directory for tree file 
  character(MAXPATHSIZE)       :: inputpath                       ! directory for spectra and metal ejection 
  character(MAXPATHSIZE)       :: filterpath                      ! directory for filters (U,B, ....) 
  integer(kind=4)              :: n_ascii_want                    ! number of timesteps for ascii (rf_mags) outputs
  integer(kind=4),allocatable  :: ascii_list(:)                   ! list of these timesteps
  integer(kind=4)              :: n_sfh_gal                       ! number of galaxies for which we want to output the sfh_tabs
  integer(kind=4),allocatable  :: sfh_gal_list(:,:)               ! IDs of these galaxies
  integer(kind=4)              :: ind_sfh_gal                     ! index of the last galaxy read in the list
  character(MAXPATHSIZE)       :: sfh_dir                         ! directory where shf files wil be dumped
#ifdef MOMAF_INPUTS
  character(MAXPATHSIZE)       :: momaf_snap_dir                  ! where to output momaf input files
#endif

  !-------------------------------------------------------------------------------------------------------------  
  ! parameters (either explicitly defined here or read in the user input file) 
  !-------------------------------------------------------------------------------------------------------------  
  integer(kind=4)              :: nsteps                           ! total number of steps in tree
  integer(kind=4)              :: nsteps_do                        ! number of steps to be done
  integer(kind=4)              :: iseed                            ! random number generator initialization seed
  integer(kind=4)              :: seed_gal_pos                     ! random seed for galaxy positions 
  integer(kind=4),parameter    :: ncomp     = 3                    ! number of galaxy components (disc,bulge,burst)
  real(kind=8),parameter       :: pi        = 3.1415927d0                   
  real(kind=8),parameter       :: gravconst = 430.1d0              ! in Mpc km^2 (10^11 M_sun)^-1 s^-2 
  real(kind=8),parameter       :: root2     = 1.4142136d0
  real(kind=8)                 :: alphapar                         ! star formation efficiency
  real(kind=8)                 :: epsilon                          ! SN feedback efficiency
  real(kind=8)                 :: eta_SN                           ! SN fraction for 1 M_sun of stars (IMF dependent)
  real(kind=8)                 :: omega_b                          ! cosmic baryons in units of critical density
  real(kind=8)                 :: omega_0                          ! cosmic matter content   "   "    "       "  
  real(kind=8)                 :: hubble                           ! reduced hubble parameter: H_0/100
  real(kind=8)                 :: Lbox_phys                        ! size of the box @ z=0, in physical Mpc
  real(kind=8)                 :: nu                               ! fraction of cold gas already in a new born halo
  real(kind=8)                 :: psi                              ! satellite-satellite merger normalization
  real(kind=8)                 :: chi                              ! galaxy merger power law (mass transfer recipe)
  real(kind=8)                 :: delay_fountain                   ! starts the fountain delta_fountain times the halo t_dyn after
  real(kind=8)                 :: duration_fountain                ! extends the fountain delta_fountain times the halo t_dyn 
  real(kind=8)                 :: meta_stable_fudge                ! parameter controling meta-stability of discs
  real(kind=8)                 :: disc_instability_threshhold      ! parameter controling stability of discs (formerly hard-ccoded  ~0.3)
#ifdef RENUMBERING
  real(kind=8)                 :: global_t_igm                    ! temperature of IGM just after reionisation (K) 
#endif

  ! hot/cold accretion parameters 
  real(kind=8)                 :: logMmin ! units (10^11 Msun)
  real(kind=8)                 :: logMmax 
  real(kind=8)                 :: ColdAccretionZ  ! metallicity (in units of mcoldz/mcold) of cold-phase accretion


#ifdef DUAL_IMF
  character(20)                :: imfname2                         ! name of second Initial Mass Function  
  real(kind=8)                 :: eta_sn2                          ! as eta_SN for the second IMF 
#if (DUAL_IMF == 1) 
  !! variables relative to implementation 1... 
  real(kind=8)                 :: z_threshold                      ! if global_redshift > z_threshold t
  real(kind=8)                 :: tdyn_threshold                   ! and tgyn < tfyn_threshold (units in the input file: Myr)
#elif (DUAL_IMF == 2)
  !! variables relative to implementation 2... 
  real(kind=8)                 :: z_threshold                      ! if global_redshift > z_threshold t
  real(kind=8)                 :: m_threshold                      ! if m_gal > m_threshold then use the top heavy IMF
#endif
#endif

#ifdef CLUMPY_SF
  real(kind=8)                 :: frac_acc_to_burst                ! fraction of stream accretion going directly to the burst component.
#endif
#ifdef TH4
  character(3)                 :: imfname                          ! name of Initial Mass Function   ! modif_imf : TH4 is 3 characters...
#else
  character(4)                 :: imfname                          ! name of Initial Mass Function 
#endif

  real(kind=8)                 :: z_reionisation                   ! redshift of reionization
  real(kind=8)                 :: alpha_reheat                     ! efficiency of IGM reheating by supernovae
  real(kind=8),parameter       :: age_for_transp = 0.1d0           ! transfer stars from burst to bulge older than this only
  integer(kind=4)              :: age_index_transp                 ! corresponding age index
  real(kind=8),parameter       :: already_star_fac = 0.0d0         ! fraction of stars in a new born galaxy
  real(kind=8),parameter       :: mu    = 0.1d0                    ! parameter common to all QSOs
  real(kind=8),parameter       :: gamma = 0.00d0                   ! idem
  real(kind=8),parameter       :: t_acc = 0.03d0                   ! idem    
  character(40),parameter      :: profile = 'TSIS'                 ! type of density profile for DM halos
  real(kind=8),parameter       :: ais  = 21.38d0                   ! parameters for non-singular isothermal 
  real(kind=8),parameter       :: a2is = 9.08d0                    ! sphere density profile (see appendix A 
  real(kind=8),parameter       :: bis  = 19.81d0                   ! of Shapiro, Iliev and Raga 1999, MNRAS:
  real(kind=8),parameter       :: b2is = 14.62d0                   ! ais = A; a2is = a^2; bis = B; b2is = b^2)
  real(kind=8),parameter       :: log_L_sun = 22.583d0             ! log10(3.83E22) = L_bol_sun [in ergs/s] ... 
                                                                   ! ... /10^11 bc our masses are [in 10^11 M_sun]
  real(kind=8),parameter       :: Z_sun = 0.0167d0                 ! solar metallicity value
  type (comp_info),parameter   :: disc_params  = comp_info(1.68d0,1.30d0,2.83d0,10.0d0)                      
  ! contains parameters common to all discs
  type (comp_info),parameter   :: bulge_params = comp_info(1.0d0+root2,root2,2.47d0,100.0d0)                   
  ! contains parameters common to all bulges
  real(kind=8),parameter       :: min_size = 1e-4                  ! minimum size for a galaxy (disc,bulge,burst,halo core) 
  ! corresponds to 100 pc --> large HII region
  real(kind=8),parameter       :: delta_m_sn   = 10.0d0            ! average mass of SN ejecta in M_sun units
  real(kind=8),parameter       :: vel_disp_fid = 20.0d0            ! fiducial velocity dispersion for feedback in km/s.
  real(kind=8),parameter       :: rel_prec = 1.e-6                 ! relative precision (single precision optimized O2)
  ! obtained using F90 precision intrisic function
  real(kind=8),parameter       :: time_res = 1e-3                  ! maximum time resolution for SSP --> 1 Myr
  integer(kind=4),parameter    :: iluv_min = 142                   ! index 142 of the alamb array corresponding to 1105 A (defined to compute UV intensity around lya)
  integer(kind=4),parameter    :: iluv_max = 162                   ! index 162 of the alamb array corresponding to 1305 A (defined to compute UV intensity around lya)
  integer(kind=4),parameter    :: iburstuv_min = 172               ! index 172 of the alamb array corresponding to 1405 A (defined to compute UV intensity around 1500A)
  integer(kind=4),parameter    :: iburstuv_max = 192               ! index 192 of the alamb array corresponding to 1605 A (defined to compute UV intensity around 1500A)
  integer(kind=4),parameter    :: ibeta_min = 196               ! index 196 of the alamb array corresponding to 1645 A (defined to compute beta slope of SEDs)
  integer(kind=4),parameter    :: ibeta_max = 262               ! index 262 of the alamb array corresponding to 2305 A (defined to compute beta slope of SEDs)
  real(kind=8)                 :: maxerror

  !-------------------------------------------------------------------------------------------------------------
  ! global variable containing all time step galaxies and halo information
  !-------------------------------------------------------------------------------------------------------------
  real(kind=8)                 :: global_redshift                  ! redshift of current timestep
  type (tshalo),allocatable    :: tsno(:)                          ! list of all halos at all timesteps 

  !-------------------------------------------------------------------------------------------------------------
  ! variables pertaining to metallicity + cooling files
  !-------------------------------------------------------------------------------------------------------------

  integer(kind=4)              :: ncool                            ! number of temperatures for cooling curves
  integer(kind=4)              :: nrej                             ! number of timesteps for metals and SSP spec evol
  integer(kind=4)              :: nmets                            ! number of metallicities for cooling curves
  integer(kind=4)              :: nelements                        ! number of heavy elements to follow (1 for simply Z)
  integer(kind=4)              :: nfile                            ! number of metallicites for metals and SSP spec
  real(kind=8),allocatable     :: templ(:)                         ! values of temperature for cooling curves
  real(kind=8),allocatable     :: coolcurlmet(:,:)                 ! cooling curves as a function of metallicity
  real(kind=8),allocatable     :: coolcurl(:)                      ! cooling curve
  real(kind=8),allocatable     :: tabmetcool(:)                    ! metallicity array for cooling curves
  real(kind=8),allocatable     :: ztab(:,:)                        ! metal ejecta for an instantaneous burst of stars
  real(kind=8),allocatable     :: gaztab(:,:)                      ! gas ejecta for an instantaneous burst of stars 
#ifdef DUAL_IMF
  real(kind=8),allocatable     :: ztab2(:,:)                      ! metal ejecta for a top heavy IMF burst
  real(kind=8),allocatable     :: gaztab2(:,:)                    ! gas ejecta for a top heavy IMF burst 
#endif
  real(kind=8),allocatable     :: timetab(:)                       ! time at each timestep for ejecta and SSP spec evol
  real(kind=8),allocatable     :: delta_timetab(:)                 ! interval between timesteps for ejecta & spec evol
  real(kind=8),allocatable     :: tabmetspec(:)                    ! metallicity table for Simple Stellar Population 

  !-------------------------------------------------------------------------------------------------------------
  ! global variables for spectra and dust
  !-------------------------------------------------------------------------------------------------------------

  ! stellar population synthesis (optical) spectra for an instantaneous burst of star formation (1 M_sun)
  integer(kind=4)              :: lmax                             ! number of stellar (UV --> near-IR) wavelengths read in.  
  real(kind=8),allocatable     :: alamb(:),dalamb(:)               ! wavelength values  
  real(kind=8),allocatable     :: sburst(:,:,:)                    ! sburst is 3d - (lambda,time,metallicity)
#ifdef DUAL_IMF
  real(kind=8),allocatable     :: sburst2(:,:,:)                  ! sburst2 is also 3d  (for second IMF)
#endif
  ! filters used to compute observed magnitudes
  integer(kind=4)              :: nftot                            ! Total number of filters
  type filter                                                      ! The filter object has several attributes:
     integer(kind=4)           :: lftot                            ! number of wavelengths in the filter
     real(kind=8)              :: cfil,aire                        ! the width + normalization 
     real(kind=8),pointer      :: wave(:),trans(:)                 ! list of wavelengths + transmission values.  
     character(20)             :: name                             ! name of filter 
  end type filter
  type (filter),allocatable    :: filt(:)

  ! dust absorption and emission properties  
  integer(kind=4)              :: find,fins,nbre                   ! nb wavelengths for ext,emis & nb spectra
  real(kind=8),allocatable     :: alambda(:),exti(:,:)             ! dust extinction curve (alamba is wavelength)
  real(kind=8),allocatable     :: lum(:),alam_ir(:),spec_ir(:,:)   ! dust emission spectra (alam_ir "       "   )

  ! optical depth
  real(kind=8),allocatable     :: tau_hires(:,:),albedo(:)         ! names are explicit
  integer(kind=4),parameter    :: nseries = 10, ntau = 120         ! nb of terms (ser exp) & nb opt depth 
  real(kind=8)                 :: series(nseries)                  ! series expansion
  real(kind=8),dimension(ntau) :: tau_tab,ext_tab,slab_tab,slab_int_tab ! arrays opt depths for 2 dust distrib 

  ! total spectra - IR(dust) + optical (stars): names are quite explicit
  real(kind=8),allocatable     :: interp_arr(:),wavelength(:),restframe_wave(:),tot_wave(:)
  real(kind=8),allocatable     :: tot_wave2(:),cosmic_extinction2(:),cosmic_extinction(:)

  ! magnitudes: absolute, rest frame, with cosmic extinction and average mag redshift gradient across the sim box
  real(kind=8),allocatable     :: amabs_rf(:),amabs_rf_noext(:),amabs_diff(:)
#ifdef MOMAF_INPUTS
  real(kind=8),allocatable     :: amabs_ce(:)
#endif

  logical(kind=4)              :: rf_flag


! useless stuff to clean up at some point 
  logical(kind=4)              :: output_fits_spectra              ! defined in baryons.dat. If true, output spectra in fits files.
  character(MAXPATHSIZE)       :: spec_dir                         ! directory where to output the spectra

! tibo - Desired output
#ifdef DESIRED_OUTPUTS
  integer(kind=4),allocatable     :: desired_ts(:)
  integer(kind=4)                 :: ndes_ts
#endif

end module GLOB_DEFS


  
