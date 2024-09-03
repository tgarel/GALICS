module shell_type

  use GLOB_DEFS
  use accretion_type

  public

  real(kind=8),parameter    :: G = 6.667d0 * 10.0**(-11)              ! SI
  real(kind=8),parameter    :: Msun = 1.9891d0 * 10.0**30             ! kg
  real(kind=8),parameter    :: M_H = 1.674d0 * 10.0**(-27)            ! kg
  real(kind=4),parameter    :: mu_neutral = 1.22d0                            ! mean molecular weight ofr fully neutral medium
  real(kind=8),parameter    :: oneyr_in_sec = 31558149.54d0
  real(kind=8),parameter    :: mass_in_kg = 1.9891*10.0**30 * 1.0d11
  real(kind=8),parameter    :: oneMpc_in_m = 1.0d6 * 3.08568025d0 * 1.0d16
  real(kind=8),parameter    :: V_f = 20000.0d0                        ! m.s^-1 = velocity dispersion of the shell = terminal velocity for the shell
  real(kind=8),parameter    :: kb = 1.3806503d-23                     ! SI
  real(kind=8),parameter    :: E_sn = 10.0d0**43                      ! Joules

contains

!*****************************************************************************************************************

  subroutine clear_shell(s)
    
    implicit none
    
    type(shell) :: s
    
    s%radius       = 0.0d0
    s%speed        = 0.0d0
    s%age          = 0.0d0
    s%zmet         = 0.0d0
    s%nh           = 0.0d0
    s%tau_dust     = 0.0d0
    s%Rfrag        = 0.0d0
    s%Rmax         = 0.0d0
    s%Mshell       = 0.0d0
    s%Mshellz      = 0.0d0
    s%tstart       = 0.0d0
    s%onoff        = .false.
    s%radius_subm1 = 0.0d0
    s%c_factor     = 1.0d0
    s%mcold_gal    = 0.0d0
    s%rgal         = 0.0d0
    s%ii           = 0
    s%oo           = 0
    s%a            = 0
    s%v_esc_halo   = 0.0d0
    s%vv           = 0
    s%c            = 0
    s%Rfinal       = 0

    return
    
  end subroutine clear_shell

!*****************************************************************************************************************

  subroutine copy_shell(s1,s2)
    
    implicit none
    
    type(shell) :: s1,s2

    s1%radius       = s2%radius
    s1%speed        = s2%speed
    s1%age          = s2%age
    s1%zmet         = s2%zmet
    s1%nh           = s2%nh
    s1%tau_dust     = s2%tau_dust
    s1%Rfrag        = s2%Rfrag
    s1%Rmax         = s2%Rmax
    s1%Mshell       = s2%Mshell
    s1%Mshellz      = s2%Mshellz
    s1%tstart       = s2%tstart
    s1%onoff        = s2%onoff
    s1%radius_subm1 = s2%radius_subm1
    s1%c_factor     = s2%c_factor
    s1%mcold_gal    = s2%mcold_gal
    s1%rgal         = s2%rgal
    s1%ii           = s2%ii
    s1%oo           = s2%oo
    s1%a            = s2%a
    s1%v_esc_halo   = s2%v_esc_halo
    s1%vv           = s2%vv
    s1%c            = s2%c
    s1%Rfinal       = s2%Rfinal

    return

  end subroutine copy_shell

!*****************************************************************************************************************

  subroutine launch_shell(t1,mcold,r,m,mstar,mfs,s,E_int_sn)

    implicit none
    
    type(shell)               :: s

    real(kind=8)              :: mcold,mstar,r,mfs,t1,E_pot,E_th,t_ff,m,E_int_sn
    real(kind=8),parameter    :: vgas = 20000.0d0                            

    ! Wind starts as soon as the injected kinetic SN energy = (E_pot + E_thermal)
    !print*,'test',E_int_sn,mfs
    E_int_sn     = E_int_sn + mfs * 1.d11 * eta_sn * E_sn * 0.2                                  ! assume 20% of the SN energy is converted into kinetic energy 
    E_pot          = G * mcold * mass_in_kg * mstar *  mass_in_kg / (r * oneMpc_in_m)        ! potential energy of the gas
    E_th           = 0.5 * mcold * mass_in_kg * vgas**2                                                   ! thermal energy of the gas
!    print*,'E_int_sn',E_int_sn
!    print*,'E_pot',E_pot
!    print*,'E_th',E_th
!    print*,'mfs',mfs*1.d11
 !   print*,mstar*1.d11
  !  print*,mcold*1.d11

    if (E_int_sn .gt. (E_pot + E_th)) then
        s%onoff = .true.
        s%tstart = t1
        !print*,'true'
    endif
    if (E_int_sn .lt. (E_pot + E_th)) then
       !print*,'false'
    endif
    s%ii = 0
  end subroutine launch_shell


    
!*****************************************************************************************************************

  subroutine compute_shell_props(gal,s,t,mcold,mcoldz,r,m,mstar,mfs,h,delta_t)

    implicit none

    type(halo)                :: h
    type(galaxy)              :: gal
    type(shell)               :: s
    real(kind=8)              :: tstart,t,mcold,mcoldz,r,m,rho,mstar,mfs,Erate,accrate,accratez,t1wind,t2wind,E_pot,E_th,t_ff,v_esc_halo,M_igm,M_igmz,delta_t
    integer(kind=4)           :: j,i
    real(kind=8),parameter    :: vgas = 20000.0d0                            
    real(kind=8),parameter    :: dt = 0.001                           

    s%zmet   = mcoldz / mcold                                                                                                  

    E_pot        = G * mcold * mass_in_kg * mstar *  mass_in_kg / (r * oneMpc_in_m)        ! potential energy of the gas
    E_th         = 0.5d0 * mcold * mass_in_kg * vgas**2                                                   ! thermal energy of the gas

    if (s%ii .eq. 0) then
       s%a         = 1
       s%mcold_gal = mcold
       s%rgal      = r
    endif
   
    s%a = s%a + 1
    !print*,'rgal',s%rgal

    rho        = s%mcold_gal * mass_in_kg / (4.0d0 / 3.0d0 * pi * (s%rgal * oneMpc_in_m)**3)                                                                   ! kg/m^3 
    s%age      = (s%a - 1) * delta_t * 10.0d0**9 * oneyr_in_sec !(t-s%tstart) * 10.0d0**9 * oneyr_in_sec       ! s
   
    !print*,'age',s%age/oneyr_in_sec/10.0d0**6
    if (s%ii .eq. 0) then 
       !Erate = (E_pot + E_th) / (10.0d0**6 * oneyr_in_sec)
       !print*,'erate',Erate
    else
       !Erate = eta_sn * E_sn * mfs * 1.d11 / (10.0d0**6 * oneyr_in_sec)           
       !print*,'Erate',Erate                                                                                            
    endif
    Erate = 10.0d0**27 * mfs * 10.0d0**11 
    !print*,'Erate',Erate

    s%ii  = 1
   
    !call shell_postshock_T(s,Tshock)

    if (s%radius .le. r * oneMpc_in_m) then
       !print*,'h1',r* oneMpc_in_m
       s%radius_subm1 = s%radius      
       s%radius   = 0.76d0 * (Erate / rho)**0.2 * s%age**(3.0d0/5.0d0)                                                                                    ! m    
       !s%speed    = 0.528d0 * (Erate / rho)**0.2 * s%age**(-2.0d0/5.0d0)                                                                                 ! m/s  
       s%speed    = 3.0d0 / 5.0d0 * s%radius / s%age                                                                                                      ! m/s  
       !print*,'v',s%speed/1000.0d0
       !print*,'r kpc',s%radius/3.08d19
       !print*,'radius',s%radius
       
       s%Mshell   = s%Mshell + 4.0d0 / 3.0d0 * pi * (s%radius**3 - s%radius_subm1**3) * rho
       s%Mshellz  = s%Mshellz + 4.0d0 / 3.0d0 * pi * (s%radius**3 - s%radius_subm1**3) * rho * s%zmet

       if (s%radius .eq. 0.0d0) then
          s%nh = 88.0d0
          print*,'radius eq 0...'
       else
          s%nh   = s%Mshell / M_H / mu_neutral / (4.0d0 * pi * s%radius * s%radius * 1.0d4) ! cm^-2
       endif

       s%tau_dust = 3.43d0 * (s%zmet / z_sun)**1.35 * (s%nh / (2.1d0*10.0**21))                                        
       
    else
       !print*,'h2'
       if (s%oo .eq. 0) then
          s%v_esc_halo = sqrt(2.0d0 * abs(4 * pi * G * h%halo_profile%rho_0 * mass_in_kg / oneMpc_in_m**3 * (oneMpc_in_m * h%halo_profile%r_c)**2 * log(s%radius / h%halo_profile%r_c) - G * h%datas%mvir * mass_in_kg / h%halo_profile%r_c / oneMpc_in_m))
          !s%Rfrag = 1.0d0 / (1.0d0 / (r * oneMpc_in_m) + (s%speed**2 - V_f**2) / (2.0d0 * G * m * mass_in_kg))            
          s%oo = 1
       endif
       
       !if (s%speed .gt. 3.0d0) then
       if (s%vv .eq. 0 ) then
          !print*,'h3'
          if (s%speed .gt. 20000.0d0) then
             s%c_factor = 1.0d0
          else
             if (s%c .eq. 0) then
                s%Rfrag = s%radius
                s%c      = 1 
             endif
             s%c_factor = (s%Rfrag / s%radius)**2
             !print*,'C_fact=',s%c_factor
          endif
          !print*,'v momentum ',s%speed
          !s%radius_subm1 = s%radius
          do j=1,1000
             if (s%speed - G * (mstar) * mass_in_kg / s%radius**2 * dt * 10.0d0**6 * oneyr_in_sec .lt. 0.0d0) then 
                s%vv = 1
                exit
             endif
             s%speed  = s%speed - G * (mstar) * mass_in_kg / s%radius**2 * dt * 10.0d0**6 * oneyr_in_sec!+ h%datas%mvir * s%radius / (h%datas%rvir * oneMpc_in_m)
             s%radius = s%radius + s%speed * dt * 10.0d0**6 * oneyr_in_sec
             !print*,'vfor',s%speed/1000.0d0
          end do
          !print*,'v momentum',s%speed/1000.
          !print*,'r momentum',s%radius/3.08d19
       else
          !print*,'vh',s%v_esc_halo/1000.0d0,s%speed/1000.0d0
          if (s%speed .lt. s%v_esc_halo) then
             call shell_freefall(t_ff,mstar,s%radius)
             accrate = s%Mshell / mass_in_kg / (t_ff / oneyr_in_sec / 10.0d0**9) 
             accratez = s%Mshellz / mass_in_kg / (t_ff / oneyr_in_sec / 10.0d0**9)
             t1wind = (s%age + t_ff) / oneyr_in_sec / 10.0d0**9
             t2wind = (s%age + 2.0d0 * t_ff) / oneyr_in_sec / 10.0d0**9
             !print*,'ACCRATE,z,t1,t2,tf',accrate,accratez,t1wind,t2wind,t_ff
             call add_accretion_event(gal%wind,accrate,accratez,t1wind,t2wind)
             s%onoff = .false.
             !print*,'KILLED'
             call clear_shell(s)
             i = 0
          else
             s%onoff = .false.
             print*,'KILLED esc'
             call clear_shell(s)
             h%igm%mcold  = h%igm%mcold + s%Mshell
             h%igm%mcoldz = h%igm%mcoldz + s%Mshellz
             i = 0
          endif
       endif
    endif
    

  end subroutine compute_shell_props


!*****************************************************************************************************************

  subroutine shell_when_merger(s1,s2)

    implicit none
    
    type(shell)        :: s1,s2 

    s1%Mshell          = s1%Mshell + s2%Mshell

    call clear_shell(s1)
    call clear_shell(s2)

  end subroutine shell_when_merger

!*****************************************************************************************************************

!!$  subroutine shell_postshock_T(s,Tshock)
!!$
!!$    implicit none
!!$
!!$    ! the radiative phase starts when the post-shock temperature becomes less than 6.*10^5K (Castor et al., 1975) => radiative cooling becomes efficient
!!$    
!!$    type(shell)  :: s
!!$
!!$    real(kind=8) :: Tshock
!!$
!!$    Tshock = 3.0d0 * 0.62 * M_H * s%speed**2 / 16.0d0 / kb
!!$
!!$  end subroutine shell_postshock_T

!*****************************************************************************************************************

  subroutine shell_freefall(t_ff,mstar,R)

    implicit none
     
    type(shell)    :: s
    real(kind=8)   :: t_ff,mstar,R

    t_ff = sqrt(pi**2 * R**3 / (8.0d0 * G * mstar * mass_in_kg))
    
  end subroutine shell_freefall

end module shell_type
