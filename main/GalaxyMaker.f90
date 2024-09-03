program GalaxyMaker

  use GALS_BARYONS
  use IDS

  implicit none

  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
  write(errunit,*)
  write(errunit,*) '      GalaxyMaker                                                      '
  write(errunit,*) '      -----------                                                      ' 
  write(errunit,*)
  write(errunit,*)
  write(errunit,*) ' first version             : S. Ninin                  (1999)          '
  write(errunit,*) ' modif.                    : J. Devriendt              (1999-2002)     '
  write(errunit,*) ' modif.                    : B. Lanzoni & S. Hatton    (2000-2001)     '
  write(errunit,*) ' ultimate version          : J. Blaizot & J. Devriendt (2002)          '
  write(errunit,*) ' pre-post-ultimate version : J. Blaizot & J. Forero    (2008)          '
  write(errunit,*) ' post-ultimate version     : J. Blaizot & J. Devriendt (2008)          '
  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'  
  write(errunit,*) 

!!$  print*,'precision =',precision(1.0),precision(1.0d0)
!!$  print*,'rel_prec = ',rel_prec, 10.**(-precision(1.0))

  if (iargc().ge.1) then
     call getarg(1,data_dir)       ! directory where to store the results
  else
     data_dir = "./"
  endif

  call init_put_baryons            ! reads some input data files, including baryons.dat 

#ifdef ASCII_TREE
  call read_tree_ascii             ! reads in the halo tree.dat file
#else
  call read_tree                   ! reads in the halo tree.dat file
#endif


  call build_halotree              ! build halo simpler tree out of complete tree from tree.dat

  call write_snaps_redshifts       ! output ascii file with the snapshot redshift list

#ifdef DEFINE_IDS
  ! jeje 
!!$  if (nsteps == nsteps_do) then 
  call define_haloids 
!!$  else
!!$     write(errunit,*) '> WARNING : I will not define IDs because nsteps_do differs from nsteps'
!!$  end if
  ! ejej
#endif

  ! jeje 
 ! call output_halo_accinfo

! tibo adds:
  call output_halo_props
  
  call evolve_galaxies             ! main loop. does all the physical evolution of baryons

#ifdef DEFINE_IDS
  ! jeje
!!$  if (nsteps == nsteps_do) then 


  call define_galids

!!$  else
!!$     write(errunit,*) '> WARNING : I will not define IDs because nsteps_do differs from nsteps'
!!$  end if
  ! ejej
#endif

  write(errunit,*)
  write(errunit,*)
  write(errunit,*) '> baryon implementation finished '
  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
  write(errunit,*)
  write(errunit,*) '      END OF GalaxyMaker                                              '
  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'  
  write(errunit,*) 

  stop
  
end program GalaxyMaker


