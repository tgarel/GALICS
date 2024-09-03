#!/usr/bin/perl
use strict;
use warnings;

# Tree files parameters
# --- WMAP1 256 option ---
# my $treeDir      = "/simulations/Horizon/256-100h-1Mpc/TreeMaker/FOF/";
# my $treeFileMin  = 1;
# my $treeFileMax  = 24;
# my $treeFileBase = "tree_file_073";
# my $nsteps       = 73; # 80;  # number of timesteps 
# --- WMAP1 512 option --- 
# my $treeDir      = "/data/blaizot/Horizon/512-100h-1Mpc/TreeMakerAcc/";
# my $treeFileMin  = 1;
# my $treeFileMax  = 15;
# my $treeFileBase = "tree_file_080";
# my $nsteps       = 80;  # number of timesteps 
# --- WMAP3 256-20Mpc option ---
# my $treeDir      = "/data/blaizot/Horizon/256-20h-1Mpc-W3/TreeMakerAcc/";
# my $treeFileMin  = 1;
# my $treeFileMax  = 14;
# my $treeFileBase = "tree_file_095";
# my $nsteps       = 95;  # number of timesteps 
# --- WMAP3 512 option --- 
# my $treeDir      = "/data/blaizot/Horizon/512-100h-1Mpc-W3/TreeMakerAcc/FOF/bigTrees/";
# my $treeFileMin  = 1;
# my $treeFileMax  = 15;
# my $treeFileBase = "tree_file_090";
# my $nsteps       = 90;  # number of timesteps 
# --- WMAP3 1024 option --- 
# my $treeDir      = "/data/blaizot/Horizon/1024-100h-1Mpc-W3/TreeMaker/alternate-take/";
# my $treeFileMin  = 1;
# my $treeFileMax  = 111;
# my $treeFileBase = "tree_file_095";
# my $nsteps       = 95;  # number of timesteps 
# --- WMAP5 1024 option --- 
# my $treeDir      = "/simu2/garel/thibault_W5_100h-1Mpc-1024/TreeMaker/";
# my $treeFileMin  = 1;
# my $treeFileMax  = 100;
# my $treeFileBase = "tree_file_021";
# my $nsteps       = 21;  # number of timesteps

# --- newSimu => z=0  WMAP5 1024 option ---                                                                                                                    
my $treeDir      = "/cral2/garel/simu3_garel/W5-100h-1Mpc_thibault_z0/DM-LSS-1024/TreeMaker200/";
my $treeFileMin  = 1;
my $treeFileMax  = 177;
my $treeFileBase = "tree_file_094";
my $nsteps       = 94;  # number of timesteps 


# GalaxyMaker run path 
# my $gmDir        = "/data/blaizot/TEST/newKennicutt/fullRun1/";
# my $gmDir        = "/simu1/data2/garel/GM-RUNS/1024-100h-1Mpc-W3/kenn_sfr/analytic_gas_evol/sf100_fd045_clumpy_facc1/";
#my $gmDir        = "/simu3/garel/GM-RUNS/1024-100h-1Mpc-W5/kenn_sfr/analytic_gas_evol/sf5_fd170_clumpy_faccBG2bst1_tdyn_lmax21_lmin10_NOTcorr_durf05_fastBGacc_2/";
# my $gmDir        = "/simu3/garel/GM-RUNS/1024-100h-1Mpc-W5/kenn_sfr/analytic_gas_evol/sf1_fd09_nsub/";
my $gmDir        = "/cral2/garel/simu3_garel/GM-RUNS/1024-100h-1Mpc-W5/NewSimu_z0/sf5_zmigr_newFilters2_beta_SFthresh21/";

my $logFile      = "log.out"; # where stderr is re-directed...

# comment/uncomment the following lines as you need:
clear_dirs();
prepare_runs();
launch_runs(1,177); # 100 # 177
check_runs();


# CLEAR DIRS
# ----------
# clear directories if they already exist.
sub clear_dirs {
    my $cmd = "rm -rf $gmDir ";
    system($cmd);
    $cmd = "mkdir $gmDir";
    system($cmd);
}


# PREPARE RUNS 
# ------------
# create directories to run GalaxyMaker and write in each of them the
# parameter files that are needed (as defined in settings.pl). Also make a
# copy of the executable GalaxyMaker in each dir.
sub prepare_runs {
    for (my $irun = $treeFileMin; $irun <= $treeFileMax; $irun++) 
    {
	my $runDir = "$gmDir"."/$irun";
	my $cmd = "mkdir $runDir $runDir"."/Snapshots $runDir"."/spectra";
	system($cmd);
	$cmd = "cp ../main/GalaxyMaker $runDir"."/GM$irun"; #give specific name to see in top
	system($cmd);
	$cmd = "cp -r ../main $runDir"."/."; # also copy code and makefile to check
	system($cmd);
# tibo
	$cmd = "cp ../main/desired_outputs.dat $runDir"."/main/."; # also copy desired outputs file
	system($cmd);

	my $num = "$irun";
	if ($irun < 100) {$num = "0"."$num";}
	if ($irun < 10) {$num = "0"."$num";}
	$cmd = "ln -s $treeDir"."/$treeFileBase".".$num"." $runDir"."/tree_file";
	system($cmd);
	
	write_baryons_dat($runDir,$nsteps);
	write_filters_dat($runDir);
	write_ascii_ts_list_dat($runDir,$nsteps);
    }
}

# LAUNCH RUNS 
# -----------
# start runs from i1 to i2 
sub launch_runs {
    my $i1 = $_[0];
    my $i2 = $_[1];

    my $fichier = "cmd_file.sh";
    open(CFILE,">$fichier");

    for (my $irun = $i1; $irun <= $i2; $irun++) 
    {
#	my $cmd = "cd $gmDir"."/$irun  ; ./GM$irun ./ > $gmDir"."/$irun"."/$logFile & ";
#	system($cmd);
	my $cmd = "cd $gmDir"."/$irun  &&  ./GM$irun ./ > $gmDir"."/$irun"."/$logFile ";
	print CFILE "$cmd \n";
    }
    close(CFILE);
}

# CHECK JOBS
# ----------
sub check_runs {
    my $nsuccess = 0;
    my $nfailed  = 0;
    my $nmissing = 0;
    my $file = "$gmDir"."/successRun.list";
    open(SFILE,">$file");
    $file = "$gmDir"."/failedRun.list";
    open(FFILE,">$file");
    $file = "$gmDir"."/missingsRun.list";
    open(MFILE,">$file");
    for (my $irun = 1; $irun <= $treeFileMax; $irun++) 
    {
	$file = "$gmDir"."/$irun"."/$logFile";
	if (-e $file) # log file exists
	{
	    my $cmd = "grep \"END OF GalaxyMaker\" $file";
	    if (!system($cmd)) {
		$nsuccess++;
		print SFILE "$irun \n";
	    }else{
		$nfailed++;
		print FFILE "$irun \n";
	    }
	} 
	else  # log file does not exist 
	{
	    $nmissing++;
	    print MFILE "$irun \n";
	}
    }
    close(SFILE);
    close(FFILE);
    close(MFILE);
    print "nb of jobs completed              : $nsuccess \n";
    print "nb of jobs failed (or incomplete) : $nfailed \n";
    print "nb of jobs missing                : $nmissing \n";
}


#------------------- SUBS ---------------------# 

# GalaxyMaker parameters (baryons.dat file)
sub write_baryons_dat {
    my $filename = "$_[0]"."/baryons.dat";
    my $nsteps   = $_[1];
    open(FILE,">$filename");
# --- WMAP1 simulations ----
#    print FILE "omega_0        = 0.300         ! amount of dark matter \n";
#    print FILE "omega_b        = 0.044         ! amount of baryons  \n";
#    print FILE "hub            = 0.700         ! reduced hubble parameter  \n";
#    print FILE "Lbox           = 142.85714     ! size of the box in comoving Mpc \n";
# --- WMAP3-100Mpc simulations ----
#     print FILE "omega_0        = 0.240         ! amount of dark matter \n";
#     print FILE "omega_b        = 0.044         ! amount of baryons  \n";
#     print FILE "hub            = 0.730         ! reduced hubble parameter  \n";
#     print FILE "Lbox           = 136.9863      ! size of the box in comoving Mpc \n";

# --- WMAP3-20Mpc simulations ----
#     print FILE "omega_0        = 0.240         ! amount of dark matter \n";
#     print FILE "omega_b        = 0.044         ! amount of baryons  \n";
#     print FILE "hub            = 0.730         ! reduced hubble parameter  \n";
#     print FILE "Lbox           = 27.39726      ! size of the box in comoving Mpc \n";
# --- WMAP5-100Mpc simulations (Komatsu+09) ----
    print FILE "omega_0        = 0.277         ! amount of dark matter \n";
    print FILE "omega_b        = 0.0459        ! amount of baryons  \n";
    print FILE "hub            = 0.702         ! reduced hubble parameter  \n";
    print FILE "Lbox           = 142.45014     ! size of the box in comoving Mpc \n";

    print FILE "alphapar       = 5.0         ! star formation efficiency \n"; # 25.00
    print FILE "epsilon        = 0.20         ! feedback efficiency  \n"; #0.045 # 0.10
    print FILE "psi            = 0.017         ! SS merging normalization \n";
    print FILE "chi            = 3.333         ! power-law for merging galaxies mass distribution \n";
    print FILE "IMF            = kenn           ! initial mass function : kennicutt, salpeter or scalo \n";
    print FILE "z_reionisation = 10.0          ! redshift at which reionisation is completed \n";
    print FILE "alpha_reheat   = 0.001         ! efficiency of SN heating on the IGM \n";
    print FILE "meta_stable    = 1.2           ! disc stability parameter \n";
    print FILE "disc_instability_threshhold  = 0.4  ! stability limit for discs \n";
    print FILE "logMmin = 1.  \n"; # 0.5
    print FILE "logMmax = 2.1   \n"; # 2.5
    print FILE "ColdAccretionZ = 0.001         ! metallicity of cold flows \n";
    print FILE "frac_acc_to_burst = 1.0        ! fraction of gas accreted directly onto the burst component \n"; # 0.2
    print FILE "delay_fountain    = 0.5           ! starts fountain delay_fountain times the halo t_dyn later \n ";  # 0.5
    print FILE "duration_fountain = 0.5           ! extends the fountain delay_fountain times the halo t_dyn \n "; # 2.5 # 1.


    print FILE "treefile       = ./tree_file \n";
    print FILE "inputpath      = /cral2/garel/simu3_garel/Inputs_GM/ \n";
#    print FILE "filterpath     = /data/blaizot/Filters/ \n";
    print FILE "filterpath     = /cral2/garel/simu3_garel/Filters_tibo/ \n";
    print FILE "spec_dir       = ./spectra/ \n";
    print FILE "momaf_snap_dir = ./Snapshots/ \n";
#    print FILE "momaf_snap_dir = /data/garel/GM-RUNS/512-100h-1Mpc-W3/kenn_sfr/sf20_feedB002_best3/momaf_cone/ \n";
    print FILE "nsteps_do      = 24   ! number of timesteps to compute \n"; #25
    #print FILE "nsteps_do      = $nsteps   ! number of timesteps to compute \n";
   close(FILE);
}

# GalaxyMaker filter list (filters.dat file)
sub write_filters_dat {
    my $filename = "$_[0]"."/filters.dat";
    open(FILE,">$filename");
    print FILE "28 \n";   # 32
    print FILE "STEIDEL_Un \n";
    print FILE "STEIDEL_G \n";
    print FILE "STEIDEL_R \n";
    print FILE "STEIDEL_I \n";
#    print FILE "NUV_galex \n";
#    print FILE "FUV_galex \n"; 
#    print FILE "SDSS_U \n";
#    print FILE "SDSS_g \n";
#    print FILE "SDSS_r \n";
#    print FILE "SDSS_i \n";
#    print FILE "SDSS_z \n";
#    print FILE "IRAS_12mic \n";
#    print FILE "IRAS_25mic \n";
#    print FILE "IRAS_60mic \n";
#    print FILE "IRAS_100mic \n";
    print FILE "SCUBA_850mic \n";
#    print FILE "MIPS_160mic \n";
#    print FILE "MIPS_24mic \n";
#    print FILE "MIPS_70mic \n";
    print FILE "UV_1500A_Gab04 \n";
#    print FILE "UV_2800A_Gab04 \n";
#    print FILE "JOHNSON_KAB \n";
#    print FILE "i775 \n";
#    print FILE "z850 \n";
#    print FILE "iwata07_1570A \n";  
#    print FILE "UV_1600Ang \n";
    print FILE "z_hudf \n";
    print FILE "Y_hudf \n";
    print FILE "V_hudf \n";
    print FILE "J_hudf \n";
    print FILE "i_hudf \n";
    print FILE "H_hudf \n";
    print FILE "B_hudf \n";
    print FILE "CFHTLS_U \n";
    print FILE "CFHTLS_G \n";
    print FILE "CFHTLS_R \n";
    print FILE "CFHTLS_I \n";
    print FILE "CFHTLS_Z \n";
    print FILE "B_sub_ouchi \n";
    print FILE "R_sub_ouchi \n";
    print FILE "V_sub_ouchi \n";
    print FILE "iprime_sub_ouchi \n";
    print FILE "zprime_sub_ouchi \n";
#    print FILE "nb387haruka \n";
#    print FILE "nb503ouchi \n";
#    print FILE "nb570ouchi \n";
#    print FILE "nb711ouchi \n";
#    print FILE "nb816ouchi \n";
#    print FILE "nb921ouchi \n";
    print FILE "F105_takuya_rebin \n";
    print FILE "F125_takuya_rebin \n";
    print FILE "F140_takuya_rebin \n";
    print FILE "F160_takuya_rebin \n";
    print FILE "F850_takuya_rebin \n";

   close(FILE);
}

# GalaxyMaker list of rest-frame mag outputs
sub write_ascii_ts_list_dat {
    my $filename = "$_[0]"."/ascii_ts_list.dat";
    my $nsteps   = $_[1];    
    open(FILE,">$filename");
    print FILE "! This file contains the list of timesteps for which rf mags will be outputed.\n";
#    print FILE "4 ! bla bla \n";
#    print FILE "8 \n";
#    print FILE "11 \n";
#    print FILE "13 \n";
#    print FILE "16 \n";
    print FILE "$nsteps  ! nb of timesteps to output ... \n"; # compute mags at last output (z=0)
   for (my $i = 1; $i <= $nsteps; $i++) {
	print FILE "$i \n";
   }
    close(FILE);
}

