module GALS_BARYONS

  use HALOS_BARYONS
  use SPECTRA_BARYONS

contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!*****************************************************************************************************************
  subroutine evolve_galaxies

  ! main routine for calculating all the halo/galaxy properties.  

    implicit none
  
    type (halo),pointer      :: h  
    integer(kind=4)          :: i,st
    integer(kind=4)          :: tdiffs,tdiffm
    real(kind=8)             :: tavant,tapres,tapres2
    character(MAXPATHSIZE)   :: fname,file
#ifdef DESIRED_OUTPUTS
    integer(kind=4)          :: ierr,i_st
#endif
    integer(kind=4)          :: do_output

    ! initialise IGM temperature calculation 
    st = 1
    call compute_IGM_props(st)

!tibo
#ifdef DESIRED_OUTPUTS
       write(file,'(a,a)') trim(data_dir),'/main/desired_outputs.dat'
       call test_stop(file)                                 ! test for files existence.                                                                       
       open(unit=desired_unit,status='old',file=file)
       read(desired_unit,*) ndes_ts
       allocate(desired_ts(ndes_ts),stat=ierr)
       if (ierr /= 0) stop '> problem allocating Desired outputs'
       do i_st = 1,ndes_ts
          read (desired_unit,'(i3)') desired_ts(i_st)
          write(errunit,'(i3)') desired_ts(i_st)
       end do
       close (desired_unit)
       !    deallocate(desired_ts)                                                             
#endif



    ! loop on all timesteps defined by user (nsteps_do is set in file baryons.dat)
    do st = 1,nsteps_do         
       write(errunit,*)          '> '
       write(errunit,'(a12,i3)') '> timestep ',st
       write(errunit,*)          '> ------------'
       call cpu_time(tavant)     ! it looks so cool when CPU time prints out on the screen ...

       ! tibo: do not evolve galaxies below z =2.8:                                                                  
!if (st <20) then

       
       ! define some timestep-related global variables
#ifdef TIBO
       global_redshift = 1.d0 / tsno(st)%aexp - 1.0d0
#else
       global_redshift = (tsno(nsteps)%aexp/tsno(st)%aexp) - 1.0d0
#endif

       seed_gal_pos    = -st -100  ! seed for random position generator (within halos)

       
       ! evolve each halo from ts-1 to ts (or during a dynamical time for new haloes)
       do i = 1,tsno(st)%nb_of_halos                ! loop on all halos of current time-step
          !print*,'NH',tsno(st)%nb_of_halos
          h => tsno(st)%liste_halos(i)              ! h is an alias for current halo
          if (h%tree%frag /= 0) cycle               ! skip halo if screwed up.
          call compute_accretion(h)                 ! compute galaxy tree (includes all galaxy properties)
          call stats_st_update(h,st)                ! update timestep statistics on halos and galaxies
       end do

       call cpu_time(tapres)
       tdiffs = nint(tapres-tavant) ;  tdiffm = tdiffs/60 ; tdiffs = mod(tdiffs,60)
       write(errunit,'(a33,i3,a3,i3,a3)') ' > computing star formation took ',tdiffm,' m ',tdiffs,' s '
       

       ! output timestep results
       call output_galaxy_info(st)                  ! output gal_results etc. for st
       call output_halo_info(st)                    ! output halo-gas_results etc for st

#ifndef BIG_RUN
       call contamination_out(st)                   ! output contamination of gals and halos
#endif
       call compute_IGM_props(st)                   ! compute temperature of IGM @ st       
       call timestep_stats(st)                      ! update ts bookkeeping

       call output_ts_info(st)                      ! write stuff about SFR, gas content, mergers etc for each ts  

#ifdef DEBUG_OUTPUTS
       if (st > 1) then
          call debug_outputs(st)
          call output_gbk(st)
       endif
#endif
      
#ifdef RECORD_SFR 
       call output_sfr_tabs(st)
#endif

       ! compute/output spectra of all gals at st
       do_output = 0
       
#ifdef DESIRED_OUTPUTS
       do i_st = 1,ndes_ts
          if (st .ne. desired_ts(i_st)) then
             do_output = 1
          else
             do_output = 0
             go to 4
          end if
       end do
#endif
       
4      continue
       
       if (do_output == 0) then ! global output, i.e output "all" info          
          call compute_spectra(st) 
       endif

#ifdef MOMAF_INPUTS
       call write_momaf_snap(st)
#endif
       if (n_sfh_gal > 0) then
          if (ind_sfh_gal <= n_sfh_gal) then
             call output_sfh_tab(st)   ! output sfh_tabs for all galaxies listed at this timestep
          end if
       end if
       call cpu_time(tapres2)          
       tdiffs = nint(tapres2-tapres) ;  tdiffm = tdiffs/60 ; tdiffs = mod(tdiffs,60)
       write(errunit,'(a33,i3,a3,i3,a3)') ' > computing spectra took        ',tdiffm,' m ',tdiffs,' s '
      
       ! output accretion histories for galaxies of the last timestep
       if (st == nsteps_do) then 
          call output_acc_stories(st)
       end if

       ! write tree and clean up to save memory
       if (st > 1) then 
          call output_gal_tree(st-1)       ! output galaxy tree 
          call remove_tree_branch(st-1)    ! erase roots of the tree
       end if
       
       call cpu_time(tapres)
       tdiffs = nint(tapres-tavant) ;  tdiffm = tdiffs/60 ; tdiffs = mod(tdiffs,60)
       write(errunit,'(a33,i3,a3,i3,a3)') ' > total CPU time for timestep   ',tdiffm,' m ',tdiffs,' s '

!tibo: dont evolve galaxies below z =2.8                                                                              
!    endif
       
    end do
 
    print*,'MaxError = ',maxerror

    ! deallocations
!tibo
#ifdef DESIRED_OUTPUTS
    deallocate(desired_ts)
#endif

    allocate(h) ; deallocate(h)
    deallocate(interp_arr,wavelength,restframe_wave,tot_wave)
    deallocate(cosmic_extinction)
    deallocate(tot_wave2,cosmic_extinction2)
    deallocate(amabs_diff,amabs_rf_noext,amabs_rf)
#ifdef MOMAF_INPUTS
    deallocate(amabs_ce)
#endif
    deallocate(templ,coolcurlmet,coolcurl,tabmetcool,ztab,gaztab,timetab,delta_timetab,tabmetspec)
    deallocate(sburst,alamb,dalamb,alambda,exti)
#ifdef DUAL_IMF 
    deallocate(sburst2,ztab2,gaztab2)
#endif

    close(ce_unit)


    return

  end subroutine evolve_galaxies

!******************************************************************************************************************
  subroutine remove_tree_branch(st)
    
    implicit none

    integer(kind=4)          :: st
    integer(kind=4)          :: j !,new_dim
    !type(tshalo),allocatable :: tshalo_temp(:)

    do j=1,tsno(st)%nb_of_halos
       call clear_halo(tsno(st)%liste_halos(j))
    end do
    deallocate(tsno(st)%liste_halos)

! jeje : i comment what's below because it isn't really useful ... or is it ? :) 
    ! deallocate and reallocate tsno 
!!$    new_dim = nsteps - st                   ! change nsteps to nsteps_do...
!!$    allocate(tshalo_temp(new_dim))
!!$    tshalo_temp = tsno(st+1:nsteps)
!!$    deallocate(tsno)
!!$    allocate(tsno(st+1:nsteps)) 
!!$    tsno(st+1:nsteps) = tshalo_temp(1:new_dim)
!!$    deallocate(tshalo_temp)

    return

  end subroutine remove_tree_branch

!******************************************************************************************************************
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

end module GALS_BARYONS
