#
# Makefile for GalaxyMaker
#

# Possible options : 
# ------------------
# - BIG_RUN      : for a normal (as opposed to resimulation) simulation
# - RENUMBERING  : ...
# - OLD_TREE_FMT : for reading tree.dat file in ancient format ... (default is new format)
# - DALE_HELOU   : Use DALE&HELOU's updated (F)IR templates instead of Julien's aging phd ... :) 
# - DUAL_IMF=x   : Use two IMFs (typically a normal one + a top-heavy one)
#		 : <x> value has to be 1 (tdyn criterion) or 2 (mass threshold)  (e.g. -DDUAL_IMF=1)
# - MOMAF_INPUTS : output MOMAF input files (which now replace the obsolete ce_mags_... files)
#                : 
# - DEFINE_IDS   : should be defaut -> define unique (and indexed) IDs for halos and galaxies.
#                : NB: changes IOs as it reads BushID in the tree files. 
#
# - NO_SN_FEEDBACK  : sets SN feedback to zero.
#
# - NO_DISC_INSTAB  : All the discs are considered stable (-> mergers are the only thing that fuels buxxxs).
#
# - NO_MERGERS      : No mergers ! (and actually no dynamical friction -> sat galaxies never reach the center of a halo ...) 
#
# - NO_REIONIZATION : All haloes can trap gas (no adiab-factor etc invoked).
# 
# - READ_HALO_ACCRETION : if set, reads accretion rates outputed by TreeMaker (ran with COMPUTE_ACCRETION option).
#
# - RECORD_SFR : if set, store and output sfr_tab (equivalent of sfh_tab without recycling ... suited for CB08 models)
#
# - DEBUG_OUTPUTS : ouput loads of info... 
#
# - CLUMPY_SF : in fact, this is about directing part of the accretion to the burst component in order to mimick clumpy SF in discs at high z ... 
#
# - ALL2BURST : put all the accretion (true_acc + fountain) into the burst
#
# - TIBO : if TreeMaker has not been run until z=0, set this option to calculate correctly the redshift
#
# - MIN_IO : sets minimum outputs 
#
# - TH4 : Use TH4 (or TH2, 8) IMF
#
# - DESIRED_OUTPUTS : write-read file containing timesteps to perform outputs
#
# - Z_MIGRATION : extension of CLUMPY_SF to model faster clump migration to burst at higher z => frac_acc2bst increases like ((1+z)/3)^1.5: a la Bournaud+11 (ApJ letters)
#

#OPTIONS = -cpp	
#OPTIONS = -DBIG_RUN -DDALE_HELOU -DDEFINE_IDS 
OPTIONS = -DBIG_RUN -DRENUMBERING -DDALE_HELOU -DDEFINE_IDS -DREAD_HALO_ACCRETION -DCLUMPY_SF -DTIBO -DNO_DISC_INSTAB -DALL2BURST -DZ_MIGRATION -DDEBUG_OUTPUTS #-DDESIRED_OUTPUTS #-DMOMAF_INPUTS #-DDEBUG_OUTPUTS # -DMIN_IO

#MACHINE    = "meso"
#F90C       = pgf90
#CC         = gcc 
#DEBUG      = -Minform=warn -g -C -Mnofprelaxed -Mchkptr
#OPTIM      = -O2 -Mnofprelaxed
#FFLAGS     = -tp k8-64 -mcmodel=medium -Mlarge_arrays -Mpreprocess -Mextend  ${OPTIM}

#MACHINE    = "democrite"   # also good for octopus
#F90C       = ifort
#CC         = gcc 
#DEBUG      = -g -CB -traceback -pg -fpe0 # -warn
#OPTIM      = -O2 -fp-model precise
#FFLAGS     = -cpp ${OPTIM}

# MACHINE    = "rene"   # also good for octopus
# F90C       = ifort
# CC         = gcc 
# DEBUG      = -g -CB -traceback -warn -pg -fpe0
# OPTIM      = -O2 
# FFLAGS     = -cpp ${OPTIM}
# OPTIONS = -DBIG_RUN -DRENUMBERING -DDALE_HELOU -DDEBUG_OUTPUTS -DNO_DISC_INSTAB  -DASCII_TREE -DMOMAF_INPUTS # -DDEFINE_IDS  -DREAD_HALO_ACCRETION  

# MACHINE    = "other"
# F90C       = gfortran
# CC         = gcc 
# DEBUG      = -g -Wall -Wextra -fbounds-check -fimplicit-none -pedantic  -Wsurprising 
# OPTIM      = -O2 
# FFLAGS     = ${OPTIM}

# MACHINE    = "lyocral"
F90C         = ifort
CC           = gcc
DEBUG        = -g -traceback -fpp -check all -debug -warn all -ftrapuv
#OPTIM        = -O2 -fp-model precise
OPTIM        = -O3 -ftz -cpp -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS       = ${OPTIM} #-free -cpp
FFLAGS      += $(OPTIONS)

#F90C    = ifort
#OPTIM   = -O2 -fp-model precise #-xavx
#DEBUG   = -g -CB -fpe0 -traceback -ftrapuv #-warn unused #-pg
#FFLAGS  = ${DEBUG}
#CC      = cc #-m32
#FFLAGS  = -cpp $(OPTIM)



BARYSUBS   =  momaf_ios.o gal_bookkeep_type.o accretion_type.o GLOB_DEFS.o UTILS_BARYONS.o INIT_BARYONS.o CHECK_BARYONS.o BOOKKEEP_BARYONS.o IO_BARYONS.o HAL_N_GAL_BARYONS.o HALOS_BARYONS.o SPECTRA_BARYONS.o GALS_BARYONS.o IDS.o GalaxyMaker.o

.SUFFIXES: .o .or .f90 


ifeq ("$(F90C)", "gfortran")
FFLAGS_COMP = $(FFLAGS) -x f95-cpp-input
else
FFLAGS_COMP = $(FFLAGS)
endif


.f90.o:	
	$(F90C) $(FFLAGS_COMP) $(OPTIONS) -c $<
.c.o:
	$(CC) -c $<

#all:	GalaxyMaker 	
all:	clean2 GalaxyMaker clean


default: all 

clean: 	
	rm -f *.o
	rm -f *.mod
	rm -f *~

clean2: 	
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -f GalaxyMaker

GalaxyMaker:   ${BARYSUBS} 
	${F90C} ${FFLAGS} ${OPTIONS} ${BARYSUBS} -o GalaxyMaker

