program build_db_tables
  
  use utils

  implicit none

  ! local general variables
  integer(kind=4)        :: ts,isubdir,EndLine
  character(3)           :: numdir
  character(MAXPATHSIZE) :: filename

  ! halo variables 
  integer(kind=8) :: BushID,TreeID,HaloID,DescendantID,FirstProgenitorID,NextProgenitorID,LastProgenitorID
  integer(kind=4) :: hno,nbgal,ih,nh,ihm
  real(kind=4)    :: mfof,mvir,mcoldgas,mtemp,mhotgas,mhotz,mcolz,rcool,rfof,rvir,tvir,fhtemp,rc,rho0,tcmt,mstar,metsinstars
  real(kind=4)    ::  x,y,z,vx,vy,vz,lx,ly,lz,ek,ep,et,spin

  EndLine = 1
  call get_info  ! define global variables from input parameter file
  
  do ts = 1,nsteps
     ! open output files (for DB only)
     if (build_halo_tables)   call open_halo_table_files(ts)
     if (build_halo_momaf) then 
        call init_halo_momaf(ts)
        ihm = 1
     end if
     if (build_galaxy_tables) call open_galaxy_table_files(ts)
     
     do isubdir = 1,nsubdirs
        numdir = irun2string(isubdir)
        
        if (build_halo_tables .or. build_halo_momaf) then 
           if (build_halo_tables) then 
              write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numdir),'/halo-gas_results.',ts
              open(unit=hgrunit,file=filename,status='old',form='formatted')
              read(hgrunit,*) ! skip header
           end if
           write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numdir),'/halo_dm_results.',ts
           open(unit=hdmunit,file=filename,status='old',form='formatted')
           read(hdmunit,*) ! skip header
           write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numdir),'/HaloIDs.',ts
           open(unit=hidunit,file=filename,status='old',form='formatted')
           read(hidunit,*) nh
           do ih = 1,nh
              !read(hidunit,'(7(i16,1x))') BushID,TreeID,HaloID,DescendantID,FirstProgenitorID,NextProgenitorID,LastProgenitorID
              read(hidunit,*) BushID,TreeID,HaloID,DescendantID,FirstProgenitorID,NextProgenitorID,LastProgenitorID
              if (build_halo_tables) then 
                 read(hgrunit,'(i8,1x,11(e14.6,1x),i6,1x,5(e14.6,1x))') &
                      hno,mfof,mvir,mcoldgas,mtemp,mhotgas,mhotz,mcolz,rfof,rvir,tvir,fhtemp,nbgal,rc,rho0,mstar,metsinstars
              end if
              read(hdmunit,'(13(e14.6,1x))') x,y,z,vx,vy,vz,lx,ly,lz,ek,ep,et,spin
              if (build_halo_tables) then 
                 call write_int16(halodbunit,HaloID)
                 call write_int8(halodbunit,ts)
                 call write_int16(halodbunit,TreeID)
                 call write_int16(halodbunit,BushID)
                 call write_int16(halodbunit,LastProgenitorID)
                 call write_int16(halodbunit,DescendantID)
                 call write_int16(halodbunit,FirstProgenitorID)
                 call write_int16(halodbunit,NextProgenitorID)
                 call write_real4(halodbunit,x)
                 call write_real4(halodbunit,y)
                 call write_real4(halodbunit,z)
                 call write_real4(halodbunit,vx)
                 call write_real4(halodbunit,vy)
                 call write_real4(halodbunit,vz)
                 call write_real4(halodbunit,lx)
                 call write_real4(halodbunit,ly)
                 call write_real4(halodbunit,lz)
                 call write_real4(halodbunit,spin)
                 call write_real4(halodbunit,mfof)
                 call write_real4(halodbunit,rfof)
                 call write_real4(halodbunit,mvir)
                 call write_real4(halodbunit,rvir)
                 call write_real4(halodbunit,tvir)
                 call write_real4(halodbunit,mcoldgas)
                 call write_real4(halodbunit,mhotgas)
                 call write_real4(halodbunit,mhotz)
                 call write_real4(halodbunit,mcolz,EndLine)
              end if ! build_halo_tables
              
              if (build_halo_momaf) then 
                 h_id(ihm) = HaloID
                 h_pos(:,ihm) = (/x,y,z/)
                 h_vel(:,ihm) = (/vx,vy,vz/)
                 ihm = ihm + 1
              end if
           end do
           close(hgrunit) ; close(hdmunit) ; close(hidunit)
        end if

     end do
     
     if (build_halo_tables) close(halodbunit)
     if (build_halo_momaf)  call dump_halo_momaf_file(ts)
     if (build_galaxy_tables) then 
        close(galtabunit)
        close(disctabunit)
        close(bulgetabunit)
        close(bursttabunit)
     end if
  end do

  stop

end program build_db_tables
