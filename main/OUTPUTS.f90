module OUTPUTS

  use GLOB_DEFS

  public

contains

  !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  !*****************************************************************************************************************

  subroutine output_galaxy_binary(g,unit)

    implicit none

    type(galaxy)    :: g
    integer(kind=4) :: unit

    write(unit) g%my_number,g%hno,g%nb_merg,g%inst_bulg,g%tgal,g%tbirth,g%disturb,g%r,g%inclination,g%p%x,g%p%y,g%p%z,&
         & g%v%x,g%v%y,g%
    call output_gal_comp_binary(g%disc,unit)
    call output_gal_comp_binary(g%bulge,unit)
    call output_gal_comp_binary(g%burst,unit)


    return

  end subroutine output_galaxy_binary

  !*****************************************************************************************************************

  subroutine output_gal_comp_binary(c,unit)
    
    implicit none 
    
    type(gal_comp)  :: c
    integer(kind=4) :: unit

#ifdef DUAL_IMF
    write(unit) c%mgal,c%rgal,c%mcold,c%tdyn,c%minstar,c%rstrip,c%sfr1,c%sfr10,c%sfr100,c%mcoldz,c%speed,c%transp, &
         & c%minstar2,c%sfr21,c%sfr210,c%sfr2100,c%transp2
#else
    write(unit) c%mgal,c%rgal,c%mcold,c%tdyn,c%minstar,c%rstrip,c%sfr1,c%sfr10,c%sfr100,c%mcoldz,c%speed,c%transp
#endif

    return
    
  end subroutine output_gal_comp_binary

  !*****************************************************************************************************************
  !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

end module OUTPUTS
