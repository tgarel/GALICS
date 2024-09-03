module gal_bookkeep_type

  ! simple book-keeping type for galaxies. We save quantities integrated over a timestep here. 

  public
  
  type gal_bookkeep
     real(kind=8) :: macc_stream
     real(kind=8) :: macc_fountain
     real(kind=8) :: mejcold
     real(kind=8) :: mejhot
     real(kind=8) :: mejout
     real(kind=8) :: mstarform
     ! TIBO
     real(kind=8) :: mstarform_bulge
     ! OBIT
     real(kind=8) :: mej_tab(300)
  end type gal_bookkeep

contains

  !*****************************************************************************************************************

  subroutine clear_gal_bookkeep(gb)
    
    implicit none
    
    type(gal_bookkeep) :: gb

    gb%macc_stream   = 0.0d0
    gb%macc_fountain = 0.0d0
    gb%mejcold       = 0.0d0
    gb%mejhot        = 0.0d0
    gb%mejout        = 0.0d0
    gb%mstarform     = 0.0d0
    gb%mstarform_bulge = 0.0d0
    return

  end subroutine clear_gal_bookkeep

  !*****************************************************************************************************************

  subroutine merge_gal_bookkeep(gb1,gb2)

    ! add gb2 to gb1 
    
    implicit none
    
    type(gal_bookkeep) :: gb1,gb2
    
    gb1%macc_stream   = gb1%macc_stream   + gb2%macc_stream
    gb1%macc_fountain = gb1%macc_fountain + gb2%macc_fountain
    gb1%mejcold       = gb1%mejcold       + gb2%mejcold
    gb1%mejhot        = gb1%mejhot        + gb2%mejhot
    gb1%mejout        = gb1%mejout        + gb2%mejout
    gb1%mstarform     = gb1%mstarform     + gb2%mstarform
    gb1%mstarform_bulge     = gb1%mstarform_bulge     + gb2%mstarform_bulge

    return

  end subroutine merge_gal_bookkeep

  !*****************************************************************************************************************

  subroutine copy_gal_bookkeep(gb1,gb2)

    ! copy gb2 into gb1 (erase gb1)
    
    implicit none
    
    type(gal_bookkeep) :: gb1,gb2
    
    gb1%macc_stream   = gb2%macc_stream
    gb1%macc_fountain = gb2%macc_fountain
    gb1%mejcold       = gb2%mejcold
    gb1%mejhot        = gb2%mejhot
    gb1%mejout        = gb2%mejout
    gb1%mstarform     = gb2%mstarform
    gb1%mstarform_bulge     = gb2%mstarform_bulge

    return

  end subroutine copy_gal_bookkeep

  !*****************************************************************************************************************

  subroutine inc_gbk_macc_stream(gb,m)
    
    implicit none
    
    type(gal_bookkeep) :: gb
    real(kind=8)       :: m

    gb%macc_stream = gb%macc_stream + m

    return

  end subroutine inc_gbk_macc_stream

  !*****************************************************************************************************************

  subroutine inc_gbk_macc_fountain(gb,m)
    
    implicit none
    
    type(gal_bookkeep) :: gb
    real(kind=8)       :: m

    gb%macc_fountain = gb%macc_fountain + m

    return

  end subroutine inc_gbk_macc_fountain

  !*****************************************************************************************************************

  subroutine inc_gbk_mej(gb,mcold,mhot,mout)
    
    implicit none
    
    type(gal_bookkeep)      :: gb
    real(kind=8),intent(in) :: mcold,mhot,mout
    
    gb%mejcold = gb%mejcold + mcold
    gb%mejhot  = gb%mejhot  + mhot
    gb%mejout  = gb%mejout  + mout
    return

  end subroutine inc_gbk_mej

  !*****************************************************************************************************************

  subroutine inc_gbk_mstarform(gb,m)

    implicit none
    
    type(gal_bookkeep) :: gb
    real(kind=8)       :: m

    gb%mstarform = gb%mstarform + m
    
    return

  end subroutine inc_gbk_mstarform

  !*****************************************************************************************************************

    subroutine inc_gbk_mstarform_bulge(gb,m)

    implicit none
    
    type(gal_bookkeep) :: gb
    real(kind=8)       :: m

    gb%mstarform_bulge = gb%mstarform_bulge + m
    
    return

  end subroutine inc_gbk_mstarform_bulge

  !*****************************************************************************************************************

end module gal_bookkeep_type

