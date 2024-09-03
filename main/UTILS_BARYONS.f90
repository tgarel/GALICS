module UTILS_BARYONS

  use GLOB_DEFS

  public

contains

!*****************************************************************************************************************
  subroutine get_central_galaxy(h,ig)

    implicit none 

    type(halo)      :: h
    integer(kind=4) :: ig,i
    real(kind=8)    :: minrad

    minrad = 10000.0
    ig     = -1
    do i = 1,h%datas%nbgal
       if (h%datas%liste_galaxies(i)%r < minrad) then
          minrad = h%datas%liste_galaxies(i)%r
          ig = i
       end if
    end do

    if (ig < 0) then 
       write(errunit,*) '> No central galaxy when there really should be one ...'
       write(errunit,*) '> Terminating ...'
       stop
    end if

    return
    
  end subroutine get_central_galaxy

  !*****************************************************************************************************************
  function cold_gas_in_halo(h)

    ! cold gas in halo is the sum of stars, ism (comp%coldgas), cold streams (gal%acc), and cold winds (gal%wind)

    implicit none

    real(kind=8)    :: cold_gas_in_halo
    type(halo)      :: h
    real(kind=8)    :: tgal
    integer(kind=4) :: ig
    
    cold_gas_in_halo = 0.0d0
    do ig = 1,h%datas%nbgal
       tgal = h%datas%liste_galaxies(ig)%tgal
       cold_gas_in_halo = cold_gas_in_halo + total_mass(h%datas%liste_galaxies(ig)) &
            & + accretion_remaining_mass(h%datas%liste_galaxies(ig)%acc,tgal) &
            & + accretion_remaining_mass(h%datas%liste_galaxies(ig)%wind,tgal)
    end do

    return
    
  end function cold_gas_in_halo

  !*****************************************************************************************************************

  !*****************************************************************************************************************
  function cold_gas_in_halo_streams(h)

    ! cold gas in halo is the sum of stars, ism (comp%coldgas), cold streams (gal%acc), and cold winds (gal%wind)

    implicit none

    real(kind=8)    :: cold_gas_in_halo_streams
    type(halo)      :: h
    real(kind=8)    :: tgal
    integer(kind=4) :: ig
    
    cold_gas_in_halo_streams = 0.0d0
    do ig = 1,h%datas%nbgal
       tgal = h%datas%liste_galaxies(ig)%tgal
       cold_gas_in_halo_streams = cold_gas_in_halo_streams + &
            & + accretion_remaining_mass(h%datas%liste_galaxies(ig)%acc,tgal) &
            & + accretion_remaining_mass(h%datas%liste_galaxies(ig)%wind,tgal)
    end do

    return
    
  end function cold_gas_in_halo_streams

  !*****************************************************************************************************************

  function accflow_mass(g,tnow)
    
    implicit none

    type(galaxy) :: g
    real(kind=8) :: tnow, tlate ! compute accretion from tnow to z=0
    real(kind=8) :: accflow_mass, acc, accz

    tlate = 1000000.0d0
    call accretion_mass(g%acc,tnow,tlate,acc,accz)
    accflow_mass = acc

    return

  end function accflow_mass

  !*****************************************************************************************************************
  subroutine locate(xx,n,x,j)

    ! subroutine which locates the position of a value x in an array xx chap 3

    implicit none

    integer(kind=4) ::  n,j,jl,ju,jm
    real(kind=8)    ::  xx(n),x

    jl = 0
    ju = n+1

    do while (ju-jl > 1) 
       jm = (ju+jl)/2
       if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    j = jl

    return

  end subroutine locate

  !*****************************************************************************************************************
  subroutine hunt(xx,n,x,jlo)

    ! same as above but with correlated values for x cf chap 3

    implicit none

    integer(kind=4) ::  inc,jhi,n,jlo,jm
    real(kind=8)    ::  xx(n),x
    logical         :: ascnd

    ascnd = (xx(n) >= xx(1))
    if ((jlo <= 0) .or. (jlo > n)) then 
       jlo = 0
       jhi = n+1
       goto 3
    end if
    inc = 1
    if ((x >= xx(jlo)) .eqv. ascnd) then 
1      jhi = jlo +inc
       if (jhi > n) then 
          jhi = n+1
       else if ((x >= xx(jhi)) .eqv. ascnd) then
          jlo = jhi
          inc = inc+inc
          goto 1
       end if
    else
       jhi = jlo
2      jlo = jhi-inc
       if (jlo < 1) then 
          jlo = 0
       else if ( (x < xx(jlo)) .eqv. ascnd) then 
          jhi = jlo
          inc = inc+inc
          goto 2
       end if
    end if
3   if (jhi-jlo == 1) then 
       if (x == xx(n)) jlo = n-1
       if (x== xx(1))  jlo = 1
       return
    end if
    jm = (jhi+jlo)/2
    if ((x>=xx(jm)) .eqv. ascnd) then 
       jlo = jm
    else 
       jhi = jm
    end if
    goto 3

  end subroutine hunt

  !*****************************************************************************************************************
  subroutine expint(n,x,res)

    ! exponential integral function from chap 6.

    implicit none

    integer(kind=4)            :: n,i,ii,nm1
    real(kind=8)               :: x,a,b,c,d,del,fact,h,psi,res
    integer(kind=4), parameter :: maxit = 100
    real(kind=8),parameter     :: eps = 1.e-7, fpmin = 1.e-30, euler = .5772156649d0

    nm1 = n-1
    if(n < 0 .or.x < 0.0d0 .or.(x ==0.0d0 .and. (n == 0 .or. n == 1))) then

       write(*,*) '> bad arguments in expint'
       read *

    else if(n == 0)then

       res = exp(-x)/x

    else if(x == 0.d0)then

       res = 1.d0/nm1

    else if(x > 1.0d0)then

       b = x+n
       c = 1.d0/fpmin
       d = 1.d0/b
       h = d
       do i=1,maxit
          a   = -i*(nm1+i)
          b   = b+2.d0
          d   = 1.d0/(a*d+b)
          c   = b+a/c
          del = c*d
          h   = h*del
          if(abs(del-1.d0) < eps)then
             res = h*exp(-x)
             return
          endif
       end do
       write(*,*) '> continued fraction failed in expint'
       read *

    else

       if(nm1 /= 0)then
          res = 1.d0/nm1
       else
          res = -log(x)-euler
       endif
       fact = 1.d0
       do i=1,maxit
          fact = -fact*x/i
          if(i /= nm1)then
             del = -fact/(i-nm1)
          else
             psi = -euler
             do ii=1,nm1
                psi = psi+1.d0/ii
             end do
             del = fact*(-log(x)+psi)
          endif
          res = res+del
          if(abs(del) < abs(res)*eps) return
       end do
       write(*,*) '> series failed in expint'
       read *
       
    endif

    return

  end subroutine expint

  !*****************************************************************************************************************
  subroutine ran1(idum,res)

    ! minimal random number generator of period 10^8 (chap 7)

    implicit none

    integer(kind=4), parameter     :: k4b = selected_int_kind(9)
    integer(kind=4), intent(inout) :: idum
    real(kind=8)                   :: res
    integer(k4b), parameter        :: ia = 16807, im = 2147483647, iq = 127773, ir = 2836
    real(kind=8), save             :: am
    integer(k4b), save             :: ix = -1, iy = -1, k

    if (idum <= 0 .or. iy < 0) then
       am   = nearest(1.0d0,-1.0d0)/im
       iy   = ior(ieor(888889999,abs(idum)),1)
       ix   = ieor(777755555,abs(idum))
       idum = abs(idum)+1
    end if

    ix   = ieor(ix,ishft(ix,13))
    ix   = ieor(ix,ishft(ix,-17))
    ix   = ieor(ix,ishft(ix,5))
    k    = iy/iq
    iy   = ia*(iy-k*iq)-ir*k
    if (iy < 0) iy = iy+im
    res = am*ior(iand(im,ieor(ix,iy)),1)

    return

  end subroutine ran1

  !*****************************************************************************************************************
  subroutine indexx(n,arr,indx)

    ! returns indices indx of an array arr such that arr(indx) is in ascending order (chap 8) 

    implicit none 

    integer(kind=4)            :: n,indx(n)
    real(kind=8)               :: arr(n)
    integer(kind=4), parameter :: m = 7, nstack = 50
    integer(kind=4)            :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
    real(kind=8)               :: a

    do j=1,n
       indx(j) = j
    enddo

    jstack = 0
    l      = 1
    ir     = n
1   if (ir-l < m) then

       do j=l+1,ir
          indxt = indx(j)
          a     = arr(indxt)
          do i=j-1,1,-1
             if (arr(indx(i)) <= a) goto 2
             indx(i+1) = indx(i)
          enddo
          i         = 0
2         indx(i+1) = indxt
       enddo
       if (jstack == 0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack-2

    else

       k         = (l+ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp
       if (arr(indx(l+1)) > arr(indx(ir))) then
          itemp     = indx(l+1)
          indx(l+1) = indx(ir)
          indx(ir)  = itemp
       endif
       if (arr(indx(l)) > arr(indx(ir))) then
          itemp    = indx(l)
          indx(l)  = indx(ir)
          indx(ir) = itemp
       endif
       if (arr(indx(l+1)) > arr(indx(l))) then
          itemp     = indx(l+1)
          indx(l+1) = indx(l)
          indx(l)   = itemp
       endif
       i     = l+1
       j     = ir
       indxt = indx(l)
       a     = arr(indxt)
3      continue
       i     = i+1
       if (arr(indx(i)) < a) goto 3
4      continue
       j     = j-1
       if (arr(indx(j)) > a) goto 4
       if (j < i) goto 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       goto 3
5      indx(l) = indx(j)
       indx(j) = indxt
       jstack  = jstack+2
       if (jstack > nstack) then
          write(*,*) 'nstack too small in indexx'
          read *
       endif
       if (ir-i+1 >= j-l) then
          istack(jstack)   = ir
          istack(jstack-1) = i
          ir               = j-1
       else
          istack(jstack)   = j-1
          istack(jstack-1) = l
          l                = i
       endif

    endif

    goto 1

  end subroutine indexx

  !*****************************************************************************************************************
  function gasdev(idum)

    ! this function is used to draw perpendicular velocity component at random

    implicit none      

    integer(kind=4) :: idum
    real(kind=8)    :: gasdev
    integer(kind=4) :: iset
    real(kind=8)    :: fac,gset,rsq,v1,v2
    save iset,gset
    data iset/0/

    rsq = 100.0d0
    if (iset == 0) then
       do while (rsq >= 1.d0 .or. rsq == 0.d0)
          call ran1(idum,v1)
          v1  = 2.d0*v1-1.d0
          call ran1(idum,v2)
          v2  = 2.d0*v2-1.d0
          rsq = v1**2+v2**2
       end do
       fac    = sqrt(-2.d0*log(rsq)/rsq)
       gset   = v1*fac
       gasdev = v2*fac
       iset   = 1
    else
       gasdev = gset
       iset   = 0
    endif

    return

  end function gasdev

  !*****************************************************************************************************************
  function xinterp(x1,y1,x2,y2,x)         
    ! simple interpolation routine.

    implicit none

    real(kind=8) :: xinterp,coeff_1,coeff_2,y2,y1,x1,x2,x

    coeff_1 = (y2-y1)/(x2-x1)
    coeff_2 = y1 - coeff_1*x1
    xinterp = coeff_1*x + coeff_2  

    return

  end function xinterp

  !*****************************************************************************************************************
  subroutine get_mgal(gal,flag) 

    ! get mgal from sum of minstar and mcold for each comp

    implicit none

    type(galaxy)           :: gal
    character(*)           :: flag
    real(kind=8)           :: mass_stars
    real(kind=8),parameter :: safety_factor = 10.0d0

    ! check that minstar ~ sum(sfhtab) for each component, and if values are close enough, equate
    if (associated(gal%disc%sfh_tab)) then
       mass_stars = sum(gal%disc%sfh_tab)
    else
       mass_stars = 0.0d0
    end if
    if (mass_stars > 0.0d0 .and. gal%disc%minstar > 0.0d0) then
       if (abs(gal%disc%minstar - mass_stars) / mass_stars > safety_factor * rel_prec) then
          write(errunit,*) '> Pb in get_mgal for disc component, called from ',flag
          write(errunit,*) '> sum(sfh_tab), minstar : ',mass_stars,gal%disc%minstar
          stop
       end if
    else if (gal%disc%minstar > rel_prec .or. mass_stars > rel_prec) then       
       write(errunit,*) '> Pb in get_mgal for disc component, called from ',flag
       write(errunit,*) '> sum(sfh_tab), minstar : ',mass_stars,gal%disc%minstar
       stop
    end if
    gal%disc%minstar = mass_stars

    if (associated(gal%bulge%sfh_tab)) then
       mass_stars = sum(gal%bulge%sfh_tab)
    else
       mass_stars = 0.0d0
    end if
    if (mass_stars > 0.0d0 .and. gal%bulge%minstar > 0.0d0) then
       if (abs(gal%bulge%minstar - mass_stars) / mass_stars > safety_factor * rel_prec) then
          write(errunit,*) '> Pb in get_mgal for bulge component, called from ',flag
          write(errunit,*) '> sum(sfh_tab), minstar : ',mass_stars,gal%bulge%minstar
          stop
       end if
    else if (gal%bulge%minstar > rel_prec .or. mass_stars > rel_prec) then       
       write(errunit,*) '> Pb in get_mgal for bulge component, called from ',flag
       write(errunit,*) '> sum(sfh_tab), minstar : ',mass_stars,gal%bulge%minstar
       stop
    end if
    gal%bulge%minstar = mass_stars

    if (associated(gal%burst%sfh_tab)) then
       mass_stars = sum(gal%burst%sfh_tab)
    else
       mass_stars = 0.0d0
    end if
    if (mass_stars > 0.0d0 .and. gal%burst%minstar > 0.0d0) then
       if (abs(gal%burst%minstar - mass_stars) / mass_stars > safety_factor * rel_prec) then
          write(errunit,*) '> Pb in get_mgal for burst component, called from ',flag
          write(errunit,*) '> sum(sfh_tab), minstar : ',mass_stars,gal%burst%minstar
          stop
       end if
    else if (gal%burst%minstar > rel_prec .or. mass_stars > rel_prec) then       
       write(errunit,*) '> Pb in get_mgal for burst component, called from ',flag
       write(errunit,*) '> sum(sfh_tab), minstar : ',mass_stars,gal%burst%minstar
       stop
    end if
    gal%burst%minstar = mass_stars

#ifdef DUAL_IMF
    if (associated(gal%disc%sfh_tab2)) then
       mass_stars = sum(gal%disc%sfh_tab2)
    else
       mass_stars = 0.0d0
    end if
    if (mass_stars > 0.0d0 .and. gal%disc%minstar2 > 0.0d0) then
       if (abs(gal%disc%minstar2 - mass_stars) / mass_stars > safety_factor * rel_prec) then
          write(errunit,*) '> Pb in get_mgal for disc component, called from ',flag
          write(errunit,*) '> sum(sfh_tab2), minstar2 : ',mass_stars,gal%disc%minstar2
          stop
       end if
    else if (gal%disc%minstar2 > rel_prec .or. mass_stars > rel_prec) then       
       write(errunit,*) '> Pb in get_mgal for disc component, called from ',flag
       write(errunit,*) '> sum(sfh_tab2), minstar2 : ',mass_stars,gal%disc%minstar2
       stop
    end if
    gal%disc%minstar2 = mass_stars

    if (associated(gal%bulge%sfh_tab2)) then
       mass_stars = sum(gal%bulge%sfh_tab2)
    else
       mass_stars = 0.0d0
    end if
    if (mass_stars > 0.0d0 .and. gal%bulge%minstar2 > 0.0d0) then
       if (abs(gal%bulge%minstar2 - mass_stars) / mass_stars > safety_factor * rel_prec) then
          write(errunit,*) '> Pb in get_mgal for bulge component, called from ',flag
          write(errunit,*) '> sum(sfh_tab2), minstar2 : ',mass_stars,gal%bulge%minstar2
          stop
       end if
    else if (gal%bulge%minstar2 > rel_prec .or. mass_stars > rel_prec) then       
       write(errunit,*) '> Pb in get_mgal for bulge component, called from ',flag
       write(errunit,*) '> sum(sfh_tab2), minstar2 : ',mass_stars,gal%bulge%minstar2
       stop
    end if
    gal%bulge%minstar2 = mass_stars

    if (associated(gal%burst%sfh_tab2)) then
       mass_stars = sum(gal%burst%sfh_tab2)
    else
       mass_stars = 0.0d0
    end if
    if (mass_stars > 0.0d0 .and. gal%burst%minstar2 > 0.0d0) then
       if (abs(gal%burst%minstar2 - mass_stars) / mass_stars > safety_factor * rel_prec) then
          write(errunit,*) '> Pb in get_mgal for burst component, called from ',flag
          write(errunit,*) '> sum(sfh_tab2), minstar2 : ',mass_stars,gal%burst%minstar2
          stop
       end if
    else if (gal%burst%minstar2 > rel_prec .or. mass_stars > rel_prec) then       
       write(errunit,*) '> Pb in get_mgal for burst component, called from ',flag
       write(errunit,*) '> sum(sfh_tab2), minstar2 : ',mass_stars,gal%burst%minstar2
       stop
    end if
    gal%burst%minstar2 = mass_stars
#endif

    call get_mcomp(gal%disc)
    call get_mcomp(gal%bulge)
    call get_mcomp(gal%burst)

    return

  end subroutine get_mgal

  !******************************************************************************************************************
  subroutine get_mhalo(h)

    implicit none

    type(halo) :: h    

    h%datas%mgaz    = h%datas%mhotgaz + h%datas%mcoldgaz 

    return

  end subroutine get_mhalo

  !******************************************************************************************************************
  subroutine get_mcomp(comp)

    ! get mgal from sum of minstar and mcold for 1 comp.  

    implicit none

    type(gal_comp) :: comp

#ifdef DUAL_IMF
    comp%mgal = comp%minstar + comp%mcold + comp%minstar2
#else
    comp%mgal = comp%minstar + comp%mcold
#endif

    return

  end subroutine get_mcomp

  !******************************************************************************************************************
  function disc_tdyn(disc)

    ! compute dynamical timescale of galaxies.  Two different formulae 
    ! depending on:  1) disk rotation velocity
    !                2) sqrt(2)*bulge freefall
    ! take the half crossing-time for circ orbits at half mass radius, r = r_d *1.68

    implicit none

    type(gal_comp) :: disc
    real(kind=8)   :: disc_tdyn

    if (disc%speed == 0.0d0) then 
       disc_tdyn    = 1e5
    else       
       disc_tdyn    = (pi * disc_params%s_to_m * disc%rgal/ disc%speed) * 977.9d0
    end if

    return

  end function disc_tdyn

  !******************************************************************************************************************
  function bulge_tdyn(bulge)

    implicit none

    type(gal_comp) :: bulge
    real(kind=8)   :: bulge_tdyn

    if (bulge%speed == 0.0d0) then 
       bulge_tdyn = 1e5           
    else       
       bulge_tdyn = (bulge_params%s_to_m * bulge%rgal / bulge%speed) * 977.9d0
    end if

    return

  end function bulge_tdyn

  !******************************************************************************************************************  
  function total_mass(gal)

    implicit none

    type(galaxy) :: gal
    real(kind=8) :: total_mass

    total_mass = gal%disc%mgal + gal%bulge%mgal + gal%burst%mgal

    return

  end function total_mass

  !******************************************************************************************************************
  function total_stellar_mass(gal)

    implicit none

    type(galaxy) :: gal
    real(kind=8) :: total_stellar_mass

    total_stellar_mass = gal%disc%minstar + gal%bulge%minstar + gal%burst%minstar

    return

  end function total_stellar_mass

  !******************************************************************************************************************
#ifdef DUAL_IMF
  function total_stellar_mass2(gal)

    implicit none

    type(galaxy) :: gal
    real(kind=8) :: total_stellar_mass2

    total_stellar_mass2 = gal%disc%minstar2 + gal%bulge%minstar2 + gal%burst%minstar2

    return

  end function total_stellar_mass2
#endif
  !******************************************************************************************************************
  function total_metal_mass(gal)

    implicit none

    type(galaxy) :: gal
    real(kind=8) :: total_metal_mass

    total_metal_mass = gal%disc%mcoldz + gal%bulge%mcoldz + gal%burst%mcoldz

    return

  end function total_metal_mass

  !******************************************************************************************************************
  function total_gal_gas_mass(gal)

    implicit none

    type(galaxy) :: gal
    real(kind=8) :: total_gal_gas_mass

    total_gal_gas_mass = gal%disc%mcold + gal%bulge%mcold + gal%burst%mcold

    return

  end function total_gal_gas_mass

  !******************************************************************************************************************
#ifdef RENUMBERING
  subroutine renumbering(renum,tsno_ts)

    ! this routine reorders halos so that renum matches the tsno 
    ! list number (needed if a big simulation is split in chunks)

    implicit none 

    integer(kind=4) :: renum,tsno_ts
    integer(kind=4) :: index, n

    if(renum.gt.0)then
       n = tsno(tsno_ts)%nb_of_halos
       call locate_int(tsno(tsno_ts)%list_halos_number,n,renum,index)
       if(.not.(renum.eq.tsno(tsno_ts)%list_halos_number(index)))then
          write(errunit,*) '> Fatal error in renum routine', renum, tsno(tsno_ts)%list_halos_number(index)
          write(errunit,*) tsno(tsno_ts)%list_halos_number(index-1),tsno(tsno_ts)%list_halos_number(index+1)
          write(errunit,*) tsno_ts,tsno(tsno_ts)%nb_of_halos
          stop
       end if
       renum = index

    endif
    return

  end subroutine renumbering

  !*****************************************************************************************************************

  subroutine locate_int(xx,n,x,j)

    ! subroutine which locates the position of a value x in an array xx chap 3

    implicit none

    integer(kind=4) ::  n,j,jl,ju,jm
    integer(kind=4)    ::  xx(n),x

    jl = 1
    ju = n+1

    do while (ju-jl > 1) 
       jm = (ju+jl)/2
       if (x > xx(jm)) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    j = jl

    return

  end subroutine locate_int
#endif
  !*****************************************************************************************************************
  subroutine det_inclination(h)

    implicit none

    integer(kind=4) :: k
    real(kind=8)    :: res
    type(halo)      :: h

    do k = 1, h%datas%nbgal
       call ran1(iseed,res)
       h%datas%liste_galaxies(k)%inclination = res ! cos of inclination angle
    end do

    return

  end subroutine det_inclination

  !*****************************************************************************************************************
  subroutine det_gal_positions(h)

    ! could do something better e.g. assume radial profile for velocities consistent with density profile ...
    ! ... tomorrow ... :) 

    implicit none

    type(halo)      :: h
    integer(kind=4) :: ig
    real(kind=8)    :: dv(3)
    type(vector)    :: dpv
    real(kind=8)    :: fac, l

#ifdef TIBO
    fac  = hubble * 1.d0 / tsno(h%my_timestep)%aexp    ! factor go from physical Mpc to comoving Mpc/h
#else   
    fac  = hubble * tsno(nsteps)%aexp / tsno(h%my_timestep)%aexp    ! factor go from physical Mpc to comoving Mpc/h
#endif

   l    = hubble * lbox_phys / 2.d0                                ! and from [-lbox/2;lbox/2] to [0; lbox]

    do ig = 1, h%datas%nbgal 
       if (h%datas%liste_galaxies(ig)%r < 0.3d0 * min_size) then  ! assume peanuts = zero...

          h%datas%liste_galaxies(ig)%p%x = h%p%x * fac + l
          h%datas%liste_galaxies(ig)%p%y = h%p%y * fac + l
          h%datas%liste_galaxies(ig)%p%z = h%p%z * fac + l
          h%datas%liste_galaxies(ig)%v%x = h%v%x
          h%datas%liste_galaxies(ig)%v%y = h%v%y
          h%datas%liste_galaxies(ig)%v%z = h%v%z          

       else

          call get_dp(dpv,h%datas%liste_galaxies(ig)%r,seed_gal_pos)   ! %r is in phys. Mpc
          h%datas%liste_galaxies(ig)%p%x = (h%p%x + dpv%x) * fac + l
          h%datas%liste_galaxies(ig)%p%y = (h%p%y + dpv%y) * fac + l
          h%datas%liste_galaxies(ig)%p%z = (h%p%z + dpv%z) * fac + l

          call gasdev(dv,seed_gal_pos)
          dv = dv * h%datas%cvel/sqrt(2.0d0) 
          h%datas%liste_galaxies(ig)%v%x = h%v%x + dv(1)
          h%datas%liste_galaxies(ig)%v%y = h%v%y + dv(2)
          h%datas%liste_galaxies(ig)%v%z = h%v%z + dv(3)

       end if
    end do



    return

  contains

    !-----------------------------------------------------------------------------------------------------------------
    subroutine get_dp(dpv,r_orbit,idum)

      implicit none

      type(vector)            :: dpv
      real(kind=8),intent(in) :: r_orbit
      real(kind=8)            :: phi, cos_theta, sin_theta
      integer(kind=4)         :: idum 

      phi       = ran3(idum) * 2.d0 * pi
      cos_theta = (ran3(idum) - 0.5d0) * 2.0d0
      sin_theta = sqrt(1.0d0 - cos_theta**2)
      dpv%z     = cos_theta             * r_orbit
      dpv%x     = sin_theta * cos(phi)  * r_orbit
      dpv%y     = sin_theta * sin(phi)  * r_orbit

      return

    end subroutine get_dp

    !-----------------------------------------------------------------------------------------------------------------
    subroutine gasdev(harvest,idum)

      implicit none

      real(kind=8),dimension(:),intent(out)      :: harvest
      real(kind=8),dimension(size(harvest))      :: rsq,v1,v2
      real(kind=8),allocatable,dimension(:),save :: g
      integer(kind=4)                            :: n,ng,nn,m,i
      integer(kind=4), save                      :: last_allocated=0
      logical, save                              :: gaus_stored=.false.
      logical, dimension(size(harvest))          :: mask
      integer(kind=4)         :: idum

      n=size(harvest)
      if (n /= last_allocated) then
         if (last_allocated /= 0) deallocate(g)
         allocate(g(n))
         last_allocated=n
         gaus_stored=.false.
      end if
      if (gaus_stored) then
         harvest=g
         gaus_stored=.false.
      else
         ng=1
         do
            if (ng > n) exit
            do i = ng,n
               v1(i) = ran3(idum)
            end do
            do i = ng,n 
               v2(i) = ran3(idum)
            end do
            v1(ng:n)=2.0d0*v1(ng:n)-1.0d0
            v2(ng:n)=2.0d0*v2(ng:n)-1.0d0
            rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
            mask(ng:n)=(rsq(ng:n)>0.0d0 .and. rsq(ng:n)<1.0d0)
            call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
            v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
            rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
            ng=ng+nn
         end do
         rsq=sqrt(-2.0d0*log(rsq)/rsq)
         harvest=v1*rsq
         g=v2*rsq
         gaus_stored=.true.
      end if

      return

    end subroutine gasdev

    !-----------------------------------------------------------------------------------------------------------------
    subroutine array_copy(src,dest,n_copied,n_not_copied)

      implicit none

      real(kind=8), dimension(:), intent(in) :: src
      real(kind=8), dimension(:), intent(out) :: dest
      integer(kind=4), intent(out) :: n_copied, n_not_copied

      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)

    end subroutine array_copy

    !-----------------------------------------------------------------------------------------------------------------
  end subroutine det_gal_positions

  !*****************************************************************************************************************
  function ran3(idum)

    implicit none

    integer(kind=4) :: idum
    integer(kind=4) :: mbig,mseed,mz
    real(kind=8)    :: ran3,fac
    parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d0/mbig)
    integer(kind=4) :: i,iff,ii,inext,inextp,k
    integer(kind=4) :: mj,mk,ma(55)
    save iff,inext,inextp,ma
    data iff /0/

    if(idum.lt.0.or.iff.eq.0)then
       iff=1
       mj=mseed-iabs(idum)
       mj=mod(mj,mbig)
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
       end do
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.mz)ma(i)=ma(i)+mbig
          end do
       end do
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56)inext=1
    inextp=inextp+1
    if(inextp.eq.56)inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.mz)mj=mj+mbig
    ma(inext)=mj
    ran3=mj*fac

    return

  end function ran3

  !*****************************************************************************************************************

  subroutine resize_sfhtab(comp,t1,t2,tbirth)

    implicit none 

    type(gal_comp)  :: comp
    real(kind=8)     :: t1,t2,dt_star,tbirth
    integer(kind=4)  :: ind_star,ind_star2
    real(kind=8),allocatable :: star_temp(:,:)

    if (comp%minstar > 0.0d0) then
       dt_star  = t2 - tbirth
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       dt_star   = t2 - t1 + timetab(size(comp%sfh_tab,dim=1))
       call locate(timetab,nrej,dt_star,ind_star2)
       ind_star2 = ind_star2 + 1
       ind_star = min(ind_star,ind_star2)
       allocate(star_temp(ind_star,nfile))
       star_temp = 0.0d0
       star_temp(1:size(comp%sfh_tab,dim=1),1:nfile) = comp%sfh_tab
       deallocate(comp%sfh_tab)
       allocate(comp%sfh_tab(ind_star,nfile))
       comp%sfh_tab = star_temp
       deallocate(star_temp)
    else
       dt_star  = t2 - t1
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       allocate(comp%sfh_tab(ind_star,nfile))
       comp%sfh_tab = 0.0d0
    end if
#ifdef DUAL_IMF
    if (comp%minstar2 > 0.0d0) then
       dt_star  = t2 - tbirth
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       dt_star   = t2 - t1 + timetab(size(comp%sfh_tab2,dim=1))
       call locate(timetab,nrej,dt_star,ind_star2)
       ind_star2 = ind_star2 + 1
       ind_star = min(ind_star,ind_star2)
       allocate(star_temp(ind_star,nfile))
       star_temp = 0.0d0
       star_temp(1:size(comp%sfh_tab2,dim=1),1:nfile) = comp%sfh_tab2
       deallocate(comp%sfh_tab2)
       allocate(comp%sfh_tab2(ind_star,nfile))
       comp%sfh_tab2 = star_temp
       deallocate(star_temp)
    else
       dt_star  = t2 - t1
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       allocate(comp%sfh_tab2(ind_star,nfile))
       comp%sfh_tab2 = 0.0d0
    end if
#endif

    return

  end subroutine resize_sfhtab

  !*****************************************************************************************************************
#ifdef DUAL_IMF
  subroutine add_new_stars_to_sfhtab(comp,mfs,metfs,use_second_imf)
#else
  subroutine add_new_stars_to_sfhtab(comp,mfs,metfs)
#endif

    implicit none
    
    type(gal_comp)  :: comp
    real(kind=8)    :: mfs,metfs
    real(kind=8)    :: metallicity
    integer(kind=4) :: met_ind
#ifdef DUAL_IMF 
    logical(kind=4) :: use_second_imf
#endif
    
#ifdef DUAL_IMF 
    if (use_second_imf) then 
       metallicity  = metfs / mfs
       call locate(tabmetspec,nfile,metallicity,met_ind)    
       if (met_ind == 0) then
          comp%sfh_tab2(1,1)     = comp%sfh_tab2(1,1) +mfs
       else if (met_ind == nfile) then
          comp%sfh_tab2(1,nfile) = comp%sfh_tab2(1,nfile) + mfs
       else
          comp%sfh_tab2(1,met_ind)       = comp%sfh_tab2(1,met_ind) + (metallicity - tabmetspec(met_ind))   / &
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
          comp%sfh_tab2(1,met_ind+1)     = comp%sfh_tab2(1,met_ind+1) + (tabmetspec(met_ind+1) - metallicity) / & 
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
       end if
       comp%minstar2 = comp%minstar2 + mfs
    else
#endif
       metallicity  = metfs / mfs
       call locate(tabmetspec,nfile,metallicity,met_ind)    
       if (met_ind == 0) then
          comp%sfh_tab(1,1)     = comp%sfh_tab(1,1) + mfs
       else if (met_ind == nfile) then
          comp%sfh_tab(1,nfile) = comp%sfh_tab(1,nfile) + mfs
       else
          comp%sfh_tab(1,met_ind)       = comp%sfh_tab(1,met_ind) + (metallicity - tabmetspec(met_ind))   / &
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
          comp%sfh_tab(1,met_ind+1)     = comp%sfh_tab(1,met_ind+1) + (tabmetspec(met_ind+1) - metallicity) / & 
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
       end if
       comp%minstar = comp%minstar + mfs
#ifdef DUAL_IMF
    end if
#endif
    return
    
  end subroutine add_new_stars_to_sfhtab
  
  !*****************************************************************************************************************
#ifdef RECORD_SFR
  subroutine resize_sfrtab(comp,t1,t2,tbirth)

    implicit none 

    type(gal_comp)  :: comp
    real(kind=8)     :: t1,t2,dt_star,tbirth
    integer(kind=4)  :: ind_star,ind_star2
    real(kind=8),allocatable :: star_temp(:,:)

    if (comp%totsfr > 0.0d0) then
       dt_star  = t2 - tbirth
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       dt_star   = t2 - t1 + timetab(size(comp%sfr_tab,dim=1))
       call locate(timetab,nrej,dt_star,ind_star2)
       ind_star2 = ind_star2 + 1
       ind_star = min(ind_star,ind_star2)
       allocate(star_temp(ind_star,nfile))
       star_temp = 0.0d0
       star_temp(1:size(comp%sfr_tab,dim=1),1:nfile) = comp%sfr_tab
       deallocate(comp%sfr_tab)
       allocate(comp%sfr_tab(ind_star,nfile))
       comp%sfr_tab = star_temp
       deallocate(star_temp)
    else
       dt_star  = t2 - t1
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       allocate(comp%sfr_tab(ind_star,nfile))
       comp%sfr_tab = 0.0d0
    end if
#ifdef DUAL_IMF
    if (comp%totsfr2 > 0.0d0) then
       dt_star  = t2 - tbirth
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       dt_star   = t2 - t1 + timetab(size(comp%sfr_tab2,dim=1))
       call locate(timetab,nrej,dt_star,ind_star2)
       ind_star2 = ind_star2 + 1
       ind_star = min(ind_star,ind_star2)
       allocate(star_temp(ind_star,nfile))
       star_temp = 0.0d0
       star_temp(1:size(comp%sfr_tab2,dim=1),1:nfile) = comp%sfr_tab2
       deallocate(comp%sfr_tab2)
       allocate(comp%sfr_tab2(ind_star,nfile))
       comp%sfr_tab2 = star_temp
       deallocate(star_temp)
    else
       dt_star  = t2 - t1
       call locate(timetab,nrej,dt_star,ind_star)
       ind_star = ind_star + 1
       allocate(comp%sfr_tab2(ind_star,nfile))
       comp%sfr_tab2 = 0.0d0
    end if
#endif

    return

  end subroutine resize_sfrtab

  !*****************************************************************************************************************
#ifdef DUAL_IMF
  subroutine add_new_stars_to_sfrtab(comp,mfs,metfs,use_second_imf)
#else
  subroutine add_new_stars_to_sfrtab(comp,mfs,metfs)
#endif

    implicit none
    
    type(gal_comp)  :: comp
    real(kind=8)    :: mfs,metfs
    real(kind=8)    :: metallicity
    integer(kind=4) :: met_ind
#ifdef DUAL_IMF 
    logical(kind=4) :: use_second_imf
#endif
    
#ifdef DUAL_IMF 
    if (use_second_imf) then 
       metallicity  = metfs / mfs
       call locate(tabmetspec,nfile,metallicity,met_ind)    
       if (met_ind == 0) then
          comp%sfr_tab2(1,1)     = comp%sfr_tab2(1,1) +mfs
       else if (met_ind == nfile) then
          comp%sfr_tab2(1,nfile) = comp%sfr_tab2(1,nfile) + mfs
       else
          comp%sfr_tab2(1,met_ind)       = comp%sfr_tab2(1,met_ind) + (metallicity - tabmetspec(met_ind))   / &
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
          comp%sfr_tab2(1,met_ind+1)     = comp%sfr_tab2(1,met_ind+1) + (tabmetspec(met_ind+1) - metallicity) / & 
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
       end if
       comp%totsfr2 = comp%totsfr2 + mfs
    else
#endif
       metallicity  = metfs / mfs
       call locate(tabmetspec,nfile,metallicity,met_ind)    
       if (met_ind == 0) then
          comp%sfr_tab(1,1)     = comp%sfr_tab(1,1) + mfs
       else if (met_ind == nfile) then
          comp%sfr_tab(1,nfile) = comp%sfr_tab(1,nfile) + mfs
       else
          comp%sfr_tab(1,met_ind)       = comp%sfr_tab(1,met_ind) + (metallicity - tabmetspec(met_ind))   / &
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
          comp%sfr_tab(1,met_ind+1)     = comp%sfr_tab(1,met_ind+1) + (tabmetspec(met_ind+1) - metallicity) / & 
               & (tabmetspec(met_ind+1) - tabmetspec(met_ind)) * mfs
       end if
       comp%totsfr = comp%totsfr + mfs
#ifdef DUAL_IMF
    end if
#endif
    return
    
  end subroutine add_new_stars_to_sfrtab

!*****************************************************************************************************************
#ifdef DUAL_IMF
  subroutine transfer_sfr_tabs(gal,sfr_transf,sfr_transf2)
#else 
  subroutine transfer_sfr_tabs(gal,sfr_transf)
#endif

    implicit none 
    
    type(galaxy)             :: gal
    real(kind=8)             :: sfr_transf
#ifdef DUAL_IMF
    real(kind=8)             :: sfr_transf2
#endif
    real(kind=8),allocatable :: transf(:,:)
    integer(kind=4)          :: size_disc,size_bur
    
    size_disc = size(gal%disc%sfr_tab,dim=1)
    allocate(transf(size_disc,nfile))
    ! stuff to transfer
    transf    = gal%disc%sfr_tab * sfr_transf
    ! remove it from disc
    gal%disc%sfr_tab = gal%disc%sfr_tab - transf
    ! add it to burst
    if (gal%burst%totsfr > 0.0d0) then 
       size_bur = size(gal%burst%sfr_tab,dim=1)
       if (size_bur < size_disc) then 
          transf(1:size_bur,:) = transf(1:size_bur,:) + gal%burst%sfr_tab
          deallocate(gal%burst%sfr_tab)
          allocate(gal%burst%sfr_tab(size_disc,nfile))
          gal%burst%sfr_tab = transf
       else
          gal%burst%sfr_tab(1:size_disc,:) = gal%burst%sfr_tab(1:size_disc,:) + transf
       end if
    else
       allocate(gal%burst%sfr_tab(size_disc,nfile))
       gal%burst%sfr_tab = transf
    end if
    gal%burst%totsfr = sum(gal%burst%sfr_tab)
    gal%disc%totsfr  = sum(gal%disc%sfr_tab)
    
    ! repeat in case of dual imf 
#ifdef DUAL_IMF
    size_disc = size(gal%disc%sfr_tab2,dim=1)
    allocate(transf(size_disc,nfile))
    ! stuff to transfer
    transf    = gal%disc%sfr_tab2 * sfr_transf
    ! remove it from disc
    gal%disc%sfr_tab2 = gal%disc%sfr_tab2 - transf
    ! add it to burst
    if (gal%burst%totsfr2 > 0.0d0) then 
       size_bur = size(gal%burst%sfr_tab2,dim=1)
       if (size_bur < size_disc) then 
          transf(1:size_bur,:) = transf(1:size_bur,:) + gal%burst%sfr_tab2
          deallocate(gal%burst%sfr_tab2)
          allocate(gal%burst%sfr_tab2(size_disc,nfile)
          gal%burst%sfr_tab2 = transf
       else
          gal%burst%sfr_tab2(1:size_disc,:) = gal%burst%sfr_tab2(1:size_disc,:) + transf
       end if
    else
       allocate(gal%burst%sfr_tab2(size_disc,nfile))
       gal%burst%sfr_tab2 = transf
    end if
    gal%burst%totsfr2 = sum(gal%burst%sfr_tab2)
    gal%disc%totsfr2  = sum(gal%disc%sfr_tab2)
#endif

    return 

  end subroutine transfer_sfr_tabs

!*****************************************************************************************************************
  
  subroutine add_sfr_tabs(comp1,comp2)
    
    ! perform comp1%sfr_tab = comp1%sfr_tab + comp2%sfr_tab
    ! and comp1%totsfr = comp1%totsfr + comp2%totsfr
    ! with dual_imf case too

    implicit none
    
    type(gal_comp)           :: comp1, comp2
    integer(kind=4)          :: size1,size2
    real(kind=8),allocatable :: star_temp(:,:)

    if (comp2%totsfr > 0.0d0) then
       size2 = size(comp2%sfr_tab,dim=1) 
       if (comp1%totsfr > 0.0d0) then
          size1 = size(comp1%sfr_tab,dim=1) 
       else
          size1 = 0
       end if
       comp1%totsfr = comp1%totsfr  + comp2%totsfr
       if (size2 > size1) then
          if (size1 > 0) then
             allocate(star_temp(size2,nfile))
             star_temp = 0.0d0
             star_temp(1:size1,:) = comp1%sfr_tab
             deallocate(comp1%sfr_tab)
             allocate(comp1%sfr_tab(size2,nfile))
             comp1%sfr_tab      = star_temp
             deallocate(star_temp)
             comp1%sfr_tab = comp1%sfr_tab + comp2%sfr_tab
          else
             allocate(comp1%sfr_tab(size2,nfile))
             comp1%sfr_tab = comp2%sfr_tab
          end if
       else
          comp1%sfr_tab(1:size2,:) = comp1%sfr_tab(1:size2,:) + comp2%sfr_tab
       end if
    end if

#ifdef DUAL_IMF
    if (comp2%totsfr2 > 0.0d0) then
       size2 = size(comp2%sfr_tab2,dim=1) 
       if (comp1%totsfr2 > 0.0d0) then
          size1 = size(comp1%sfr_tab2,dim=1) 
       else
          size1 = 0
       end if
       comp1%totsfr2 = comp1%totsfr2  + comp2%totsfr2
       if (size2 > size1) then
          if (size1 > 0) then
             allocate(star_temp(size2,nfile))
             star_temp = 0.0d0
             star_temp(1:size1,:) = comp1%sfr_tab2
             deallocate(comp1%sfr_tab2)
             allocate(comp1%sfr_tab2(size2,nfile))
             comp1%sfr_tab2       = star_temp
             deallocate(star_temp)
             comp1%sfr_tab2       = comp1%sfr_tab2 + comp2%sfr_tab2
          else
             allocate(comp1%sfr_tab2(size2,nfile))
             comp1%sfr_tab2 = comp2%sfr_tab2
          end if
       else
          comp1%sfr_tab2(1:size2,:) = comp1%sfr_tab2(1:size2,:) + comp2%sfr_tab2
       end if
    end if
#endif

    return

  end subroutine add_sfr_tabs

#endif
!*****************************************************************************************************************

  subroutine add_sfh_tabs(comp1,comp2)

    implicit none
    
    type(gal_comp)           :: comp1, comp2
    integer(kind=4)          :: size1,size2
    real(kind=8),allocatable :: star_temp(:,:)
    
    if (comp2%minstar > 0.0d0) then
       size2 = size(comp2%sfh_tab,dim=1) 
       if (comp1%minstar > 0.0d0) then
          size1 = size(comp1%sfh_tab,dim=1) 
       else
          size1 = 0
       end if
       comp1%minstar = comp1%minstar  + comp2%minstar
       if (size2 > size1) then
          if (size1 > 0) then
             allocate(star_temp(size2,nfile))
             star_temp = 0.0d0
             star_temp(1:size1,:) = comp1%sfh_tab
             deallocate(comp1%sfh_tab)
             allocate(comp1%sfh_tab(size2,nfile))
             comp1%sfh_tab      = star_temp
             deallocate(star_temp)
             comp1%sfh_tab = comp1%sfh_tab + comp2%sfh_tab
          else
             allocate(comp1%sfh_tab(size2,nfile))
             comp1%sfh_tab = comp2%sfh_tab
          end if
       else
          comp1%sfh_tab(1:size2,:) = comp1%sfh_tab(1:size2,:) + comp2%sfh_tab
       end if
    end if

#ifdef DUAL_IMF
    if (comp2%minstar2 > 0.0d0) then
       size2 = size(comp2%sfh_tab2,dim=1) 
       if (comp1%minstar2 > 0.0d0) then
          size1 = size(comp1%sfh_tab2,dim=1) 
       else
          size1 = 0
       end if
       comp1%minstar2 = comp1%minstar2  + comp2%minstar2
       if (size2 > size1) then
          if (size1 > 0) then
             allocate(star_temp(size2,nfile))
             star_temp = 0.0d0
             star_temp(1:size1,:) = comp1%sfh_tab2
             deallocate(comp1%sfh_tab2)
             allocate(comp1%sfh_tab2(size2,nfile))
             comp1%sfh_tab2      = star_temp
             deallocate(star_temp)
             comp1%sfh_tab2 = comp1%sfh_tab2 + comp2%sfh_tab2
          else
             allocate(comp1%sfh_tab2(size2,nfile))
             comp1%sfh_tab2 = comp2%sfh_tab2
          end if
       else
          comp1%sfh_tab2(1:size2,:) = comp1%sfh_tab2(1:size2,:) + comp2%sfh_tab2
       end if
    end if
#endif

    return

  end subroutine add_sfh_tabs

  !*****************************************************************************************************************
  ! Tibo
  subroutine det_halo_positions(h)
    
    implicit none
    
    type(halo)      :: h
    real(kind=8)    :: fac, l
    
#ifdef TIBO
    fac  = hubble * 1.d0 / tsno(h%my_timestep)%aexp                 ! factor go from physical Mpc to comoving Mpc/h
#else   
    fac  = hubble * tsno(nsteps)%aexp / tsno(h%my_timestep)%aexp    ! factor go from physical Mpc to comoving Mpc/h
#endif
    
    l    = hubble * lbox_phys / 2.d0                                ! and from [-lbox/2;lbox/2] to [0; lbox]
    
    h%p%x =  h%p%x * fac + l
    h%p%y =  h%p%y * fac + l
    h%p%z =  h%p%z * fac + l

    return

  end subroutine det_halo_positions
!*****************************************************************************************************************
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
end module UTILS_BARYONS




