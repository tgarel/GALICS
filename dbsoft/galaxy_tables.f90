program galaxy_tables

  use utils

  implicit none 

  character(200) :: filename,dbdir,numrun
  integer(kind=4) :: nruns,irun
  integer(kind=4) :: nsteps,ts,ng,ig
  integer(kind=8) :: BushID,TreeID,HaloID,DescendantID,FirstProgenitorID,NextProgenitorID,LastProgenitorID,GalaxyID
  real(kind=4)    :: discmgal,discmcold,discsfr,discminstar,discmcoldz,discrgal,disctdyn,discspeed,disctransp,discinstsfr
  real(kind=4)    :: bulgemgal,bulgemcold,bulgesfr,bulgeminstar,bulgemcoldz,bulgergal,bulgetdyn,bulgespeed,bulgetransp,bulgeinstsfr
  real(kind=4)    :: burstmgal,burstmcold,burstsfr,burstminstar,burstmcoldz,burstrgal,bursttdyn,burstspeed,bursttransp,burstinstsfr
  real(kind=4)    :: disturb,rorbit,cvel,inclination,tbirth,x,y,z,bolum,irbolum
  integer(kind=4) :: halonum,nbmerg
  real(kind=4)    :: totmgal,totmcold,totsfr,totminstar,totmcoldz,totrgal,totinstsfr

  nsteps = 80
  nruns  = 8
  write(dbdir,'(a)') "/data9/blaizot/512-100h-1Mpc-W3/DBTables/"
  write(gmdir,'(a)') "/data9/blaizot/512-100h-1Mpc-W3/DBGM/"
!!$  nsteps = 73
!!$  nruns  = 4
!!$  write(dbdir,'(a)') "/data9/blaizot/256-100h-1Mpc/DBTables/"
!!$  write(gmdir,'(a)') "/data9/blaizot/256-100h-1Mpc/DBGM/"

  ! write votable
!!$  write(filename,'(a,a,a)') trim(dbdir),'galaxy_table.xml'
!!$  open(unit=metaunit,file=filename,status='unknown',form='formatted')
!!$  call write_xml_header(metaunit)

  ! read / write data
  do ts = 1,nsteps 
     write(errunit,*) ts

     ! open galaxy table file (and add to script)
     write(filename,'(a,a,i3.3,a)') trim(dbdir),'galaxy_table_',ts,'.db'
     open(unit=galtabunit,file=filename,status='unknown',form='formatted')
     write(filename,'(a,a,i3.3,a)') trim(dbdir),'disc_table_',ts,'.db'
     open(unit=disctabunit,file=filename,status='unknown',form='formatted')
     write(filename,'(a,a,i3.3,a)') trim(dbdir),'bulge_table_',ts,'.db'
     open(unit=bulgetabunit,file=filename,status='unknown',form='formatted')
     write(filename,'(a,a,i3.3,a)') trim(dbdir),'burst_table_',ts,'.db'
     open(unit=bursttabunit,file=filename,status='unknown',form='formatted')

     ! loop on sub-runs to fill the db file
     do irun = 1,nruns
        numrun = irun2string(irun)
        write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numrun),'/GalaxyIDs.',ts
        open(unit=galidunit,file=filename,status='old',form='formatted')
        read(galidunit,*) ng
        write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numrun),'/gal_results.',ts
        open(unit=galresunit,file=filename,status='old',form='formatted')
        read(galresunit,*) ! skip header
        write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numrun),'/gal_sfr.',ts
        open(unit=instsfrunit,file=filename,status='old',form='formatted')
        read(instsfrunit,*) ! skip header        
        write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numrun),'/bol_lum.',ts
        open(unit=lumunit,file=filename,status='old',form='formatted')
        read(lumunit,*) ! skip header        
        do ig = 1,ng
           read(galidunit,'(8(i16,1x))') BushID,TreeID,GalaxyID,DescendantID,FirstProgenitorID,NextProgenitorID,LastProgenitorID,HaloID
           ! read stupid gal_results format 
           bulgemgal=0.0; bulgemcold=0.0; bulgesfr=0.0; bulgeminstar=0.0; bulgemcoldz=0.0; bulgergal=0.0; bulgetdyn=0.0; bulgespeed=0.0; bulgetransp=0.0
           burstmgal=0.0; burstmcold=0.0; burstsfr=0.0; burstminstar=0.0; burstmcoldz=0.0; burstrgal=0.0; bursttdyn=0.0; burstspeed=0.0; bursttransp=0.0
           read(galresunit,fmt=12,ADVANCE='NO', EOR=700) discmgal,discmcold,discsfr,discminstar,discmcoldz,discrgal,disctdyn,discspeed,disctransp
           read(galresunit,fmt=11,ADVANCE='NO', EOR=700) disturb,rorbit,halonum,cvel,inclination,nbmerg,tbirth
           read(galresunit,fmt=12,ADVANCE='NO', EOR=700) bulgemgal,bulgemcold,bulgesfr,bulgeminstar,bulgemcoldz,bulgergal,bulgetdyn,bulgespeed,bulgetransp
           read(galresunit,fmt=12,ADVANCE='YES') burstmgal,burstmcold,burstsfr,burstminstar,burstmcoldz,burstrgal,bursttdyn,burstspeed,bursttransp
700        continue
11         format (2(e14.6,1x),i7,1x,2(e14.6,1x),1x,i6,1x,e14.6) 
12         format (9(e14.6,1x))
           ! read instantaneous SFR 
           read(instsfrunit,'(3(e14.6))') discinstsfr,bulgeinstsfr,burstinstsfr
           ! read positions
           x = 0.0; y = 0.0; z = 0.0
           ! read bolometric luminosities
           read(lumunit,'(2(e14.6,1x))') bolum,irbolum

           ! define global quantities
           totmgal    = discmgal    + bulgemgal    + burstmgal
           totmcold   = discmcold   + bulgemcold   + burstmcold
           totsfr     = discsfr     + bulgesfr     + burstsfr
           totinstsfr = discinstsfr + bulgeinstsfr + burstinstsfr
           totminstar = discminstar + bulgeminstar + burstminstar
           totmcoldz  = discmcoldz  + bulgemcoldz  + burstmcoldz
           if (totminstar > 0.0) then 
              totrgal    = (discminstar*discrgal + bulgeminstar*bulgergal + burstminstar*burstrgal) / totminstar
           else
              totrgal = 0.0
           end if
           
           ! write it all to the table files
           call write_int16(galtabunit,GalaxyID)
           call write_int16(galtabunit,HaloID)
           call write_int8(galtabunit,ts)
           call write_int16(galtabunit,TreeID)
           call write_int16(galtabunit,BushID)
           call write_int16(galtabunit,LastProgenitorID)
           call write_int16(galtabunit,DescendantID)
           call write_int16(galtabunit,FirstProgenitorID)
           call write_int16(galtabunit,NextProgenitorID)
           call write_real4(galtabunit,totmgal)
           call write_real4(galtabunit,totminstar)
           call write_real4(galtabunit,totmcold)
           call write_real4(galtabunit,totmcoldz)
           call write_real4(galtabunit,totsfr)
           call write_real4(galtabunit,totinstsfr)
           call write_real4(galtabunit,totrgal)
           call write_real4(galtabunit,rorbit)
           call write_real4(galtabunit,inclination)
           call write_real4(galtabunit,tbirth)
           call write_real4(galtabunit,disturb)
           call write_real4(galtabunit,x)
           call write_real4(galtabunit,y)
           call write_real4(galtabunit,z)
           call write_real4(galtabunit,bolum)
           call write_real4(galtabunit,irbolum)

           call write_int16(disctabunit,GalaxyID)
           call write_int16(disctabunit,HaloID)
           call write_int8(disctabunit,ts)
           call write_real4(disctabunit,discmgal)
           call write_real4(disctabunit,discminstar)
           call write_real4(disctabunit,discmcold)
           call write_real4(disctabunit,discmcoldz)
           call write_real4(disctabunit,discsfr)
           call write_real4(disctabunit,discinstsfr)
           call write_real4(disctabunit,discrgal)
           call write_real4(disctabunit,disctdyn)
           call write_real4(disctabunit,discspeed)

           call write_int16(bulgetabunit,GalaxyID)
           call write_int16(bulgetabunit,HaloID)
           call write_int8(bulgetabunit,ts)
           call write_real4(bulgetabunit,bulgemgal)
           call write_real4(bulgetabunit,bulgeminstar)
           call write_real4(bulgetabunit,bulgemcold)
           call write_real4(bulgetabunit,bulgemcoldz)
           call write_real4(bulgetabunit,bulgesfr)
           call write_real4(bulgetabunit,bulgeinstsfr)
           call write_real4(bulgetabunit,bulgergal)
           call write_real4(bulgetabunit,bulgetdyn)
           call write_real4(bulgetabunit,bulgespeed)

           call write_int16(bursttabunit,GalaxyID)
           call write_int16(bursttabunit,HaloID)
           call write_int8(bursttabunit,ts)
           call write_real4(bursttabunit,burstmgal)
           call write_real4(bursttabunit,burstminstar)
           call write_real4(bursttabunit,burstmcold)
           call write_real4(bursttabunit,burstmcoldz)
           call write_real4(bursttabunit,burstsfr)
           call write_real4(bursttabunit,burstinstsfr)
           call write_real4(bursttabunit,burstrgal)
           call write_real4(bursttabunit,bursttdyn)
           call write_real4(bursttabunit,burstspeed)

        end do
        close(galidunit); close(galresunit); close(instsfrunit); close(lumunit)
     end do
     close(galtabunit); close(disctabunit); close(bulgetabunit); close(bursttabunit)

!!$     write(filename,'(a,a,a,i3.3,a)') '<LINK content-role="location" content-type="ascii" href="file:',trim(dbdir),'halo_table_',ts,'.db"/>'
!!$     write(metaunit,'(a)') filename
       
  end do
  
!!$  call write_xml_footer(metaunit)
!!$  close(metaunit)

  stop

end program galaxy_tables
