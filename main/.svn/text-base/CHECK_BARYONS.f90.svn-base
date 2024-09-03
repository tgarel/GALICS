module CHECK_BARYONS

  ! contains routines used to check calculations
  
  use GLOB_DEFS
  use UTILS_BARYONS
  
  public
  
contains

!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!******************************************************************************************************************
  subroutine check_tdyn(name,comp)

    implicit none

    type(gal_comp)           :: comp
    character(*), intent(in) :: name

    if (comp%tdyn <= 0.0d0 .and. comp%mgal > 0.0d0 .and. comp%mcold > 0.0d0) then
       write(errunit,*) '> zero tdyn in ',name,':',comp%tdyn,comp%mgal,comp%mcold,comp%speed,comp%rgal
       stop
    end if

    return 

  end subroutine check_tdyn

!*****************************************************************************************************************
  subroutine check_orbit_rfof(h)

  ! Just check here that all h's galaxies' orbits are within rfof (actually set orbit radius to rfof 
  ! if larger than this).

    implicit none

    type(halo)      :: h
    integer(kind=4) :: i

    do i= 1,h%datas%nbgal
       h%datas%liste_galaxies(i)%r = min(h%datas%liste_galaxies(i)%r,h%rfof)
    end do

    return

  end subroutine check_orbit_rfof

!******************************************************************************************************************
  subroutine check_orbit_radius(g)

  ! routine to check if the galaxy's orbital radius is less than the radius of its components. 
  ! if so, set orbital radius to zero.
  
    implicit none

    type(galaxy) :: g    
               
    if (g%r <= disc_params%s_to_m  * g%disc%rgal ) then
       g%r = 0.0d0
       return 
    end if
    if (g%r <= bulge_params%s_to_m * g%bulge%rgal) then 
       g%r = 0.0d0
       return
    end if
    if (g%r <= bulge_params%s_to_m * g%burst%rgal) then
       g%r = 0.0d0
       return
    end if

    return

  end subroutine check_orbit_radius

!******************************************************************************************************************
    subroutine check_radius(name,comp)

      implicit none

      type(gal_comp)           :: comp
      character(*), intent(in) :: name

      if (comp%rgal == 0.0d0.and. comp%mgal /= 0.0d0) then 
         write(errunit,*) '> ',name,' radius zero:',comp%rgal,comp%transp,comp%mgal
         stop
      end if

      return 

    end subroutine check_radius
  
!*****************************************************************************************************************  
    subroutine check_mcold(name,comp)

      implicit none

      type(gal_comp)           :: comp
      character(*), intent(in) :: name

      if (comp%mcold < 0.0d0) then 
         write(errunit,*) '> -ve mcold in ',name,':',comp%mcold
         stop
      end if

      return 

    end subroutine check_mcold

!*****************************************************************************************************************
  subroutine check_eps(a,b)

  ! this sub is necessary cos we have 1-X in the previous sub, where X 
  ! can be very small.  1 - X can therefore be beyond machine precision 
  ! to compute.  we therefore round up or down.  

    implicit none 

    real(kind=8),parameter :: one = 1.0d0 - rel_prec
    real(kind=8)           :: a,b

    if (a < rel_prec) then 
       a = 0.0d0
       b = 1.0d0
    else if (a > one) then 
       a = 1.0d0
       b = 0.0d0
    end if

    return 

  end subroutine check_eps

!*****************************************************************************************************************
  subroutine check_newgal(comp,name)

    implicit none

    character(4),intent(in) :: name
    type(gal_comp)          :: comp
    real(kind=8)            :: x,mstar

    if (comp%mcold < 0.0d0) then 
       write(errunit,*) '> mcold -ve in ',name,':', comp%mcold
       stop
    end if

#ifdef DUAL_IMF 
    mstar = comp%minstar + comp%minstar2 
#else
    mstar = comp%minstar
#endif
    if (comp%mcold + mstar /= comp%mgal) then 
       if (comp%mgal == 0.0d0) then 
          write(errunit,*) '> mgal = 0.0 in ',name,':',comp%mcold,mstar,comp%mgal
          stop
       else
          x = (comp%mcold + mstar - comp%mgal)/comp%mgal
          if (abs(x) > rel_prec) then 
             write(errunit,*)  '> mcold+mstar /= mgal in ',name,':',comp%mcold,mstar,comp%mgal,x
             stop
          end if
       end if
    end if

    if (comp%mcoldz > comp%mcold) then
       if (comp%mcold > rel_prec) then
          x = (comp%mcoldz - comp%mcold)/comp%mcold
          if (x > rel_prec) then 
             write(errunit,*)  '> 1st if mcoldz > mcold in ',name,':',comp%mcold,mstar,comp%mgal,comp%mcoldz
             stop
          end if
       else
          if (comp%mcoldz > rel_prec) then
             write(errunit,*) '> mcoldz > mcold in ',name,':',comp%mcold,mstar,comp%mgal,comp%mcoldz
             stop
           else
              ! numerical error due to a difference in single precision --> reset to 0
              comp%mcold  = 0.0d0  
              comp%mcoldz = 0.0d0
           endif   
       endif   
    end if

    if (comp%mgal < 0.0d0) then 
       write(errunit,*) '> mgal < 0 in ',name,':',comp%mcold,mstar,comp%mgal,comp%mcoldz
       stop
    end if

    return

  end subroutine check_newgal

!*****************************************************************************************************************
!!$  subroutine check_halo_gas(h)
!!$
!!$    implicit none
!!$
!!$    type(halo)             :: h
!!$    integer(kind=4)        :: j
!!$    real(kind=8)           :: x,ncf,ms,mets
!!$    type(galaxy)           :: g
!!$
!!$    ncf = 0.0 ; ms = 0.0 ; mets = 0.0 
!!$    do j = 1,h%datas%nbgal
!!$       ncf  = ncf  + h%datas%liste_galaxies(j)%cooling_frac
!!$       ms   = ms   + total_mass(h%datas%liste_galaxies(j))
!!$       mets = mets + total_metal_mass(h%datas%liste_galaxies(j))
!!$    end do
!!$
!!$    if (ncf > 1.0) then 
!!$       write(errunit,*) '> cooling fraction > 1 in check_halo_gas',ncf
!!$       write(errunit,*) (h%datas%liste_galaxies(j)%cooling_frac,j=1,h%datas%nbgal)
!!$       stop
!!$    end if
!!$
!!$    if (ms > 0.0) then 
!!$       x = abs(ms-h%datas%mcoldgaz)/ms
!!$    else
!!$       if (h%datas%mcoldgaz > 0.0) then
!!$          write(errunit,*) '> fatal error in check_halo_gas: ms = 0.0 but halo has cold gas'
!!$          write(errunit,*) '> nb of galaxies in the halo :',h%datas%nbgal
!!$          ms = -1.0
!!$          ms = sqrt(ms)
!!$          stop
!!$       else
!!$          x = 0.0
!!$       endif   
!!$    endif   
!!$
!!$    ! check numerical precision 
!!$    if (x > ncomp*h%datas%nbgal*rel_prec) then 
!!$       ! assume this is the maximum relative error we can tolerate in a timestep 
!!$       ! (for a halo containing 1000 galaxies this is 3e-3  i.e. 0.3 %) 
!!$       write(errunit,*) '> fatal error in check_halo_gas: cooling rounding error too big for halo: ',h%my_number
!!$       write(errunit,*) '> number of galaxies in this halo :',h%datas%nbgal
!!$       write(errunit,*) '> relative error, mcoldgas_gal and mcoldgas_halo: ',x,ms,h%datas%mcoldgaz
!!$       g = h%datas%liste_galaxies(1)
!!$       ms = -1.0
!!$       ms = sqrt(ms)
!!$       stop
!!$    else
!!$       ! reset the value of metals and cold gas in halo so that errors cannot accumulate over timesteps
!!$       h%datas%mcoldgaz = ms
!!$       h%datas%mcoldz   = mets
!!$       call get_mhalo(h)
!!$    end if
!!$
!!$    return
!!$
!!$  end subroutine check_halo_gas

!*****************************************************************************************************************
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

end module CHECK_BARYONS
