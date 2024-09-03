program momaf_wrapup

  ! piece of code which reads all snapshot files and puts them in one file per type per snapshot.

  use utils

  implicit none

  character(MAXPATHSIZE)  :: filename,snapdir,numrun,progname
  character(3)   :: nsteps_str, nruns_str
  integer(kind=4) :: nruns,irun,ng,get_ng,ngtot,iglast,ifilt
  integer(kind=4) :: ts
  integer(kind=8),allocatable :: galid(:),haloid(:),tgalid(:),thaloid(:)
  real(kind=4),allocatable    :: pos(:,:),vel(:,:),inc(:),tpos(:,:),tvel(:,:),tinc(:),m(:),dm(:),tm(:),tdm(:)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! We read on the command line the number of steps, the number of runs
  !! snapshot directory, the galaxymaker directory.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  if (iargc().le.3) then
!!$     call getarg(0,progname)
!!$     write(*,*)"Usage: ", adjustl(trim(progname)), " <nb steps> <nb runs> <snapdir> <gmdir>"
!!$     write(*,*)" - nb steps: number of timesteps (3 characters max)"
!!$     write(*,*)" - nb runs: number of runs (3 characters max)"
!!$     write(*,*)" - snapdir: the output directory"
!!$     write(*,*)" - gmdir: the input directory, each run must be under <gmdir>/<run number>"
!!$     write(*,*)""
!!$     write(*,*)"Note that:"
!!$     write(*,*)" - the file filters.dat must me present in <gmdir>"
!!$     write(*,*)" - each run must be present under <gmdir>/<run number>"
!!$
!!$     stop
!!$  endif
!!$
!!$  call getarg(1,nsteps_str)  ! number of steps
!!$  read(nsteps_str,'(i3)') nsteps
!!$  call getarg(2,nruns_str)   ! number of runs
!!$  read(nruns_str,'(i3)') nruns
!!$  call getarg(3,snapdir)     ! snapshot directory
!!$  call getarg(4,gmdir)       ! galaxymaker directory

!   nsteps = 73
!   nruns  = 4
!   write(snapdir,'(a)') "/data/blaizot/Simulations/256-100h-1Mpc/GalaxyMaker/Cattaneo/Snapshots/"
!   write(gmdir,'(a)')   "/data/blaizot/Simulations/256-100h-1Mpc/GalaxyMaker/Cattaneo/"
  nsteps = 29
  nruns  = 177

  write(snapdir,'(a)') "/scratch/garel/Galics/GM-RUNS/1024-100h-1Mpc-W5/NewSimu_z0/sf5_zmigr_newFilters2_z2_BlueMuse/Snapshots/"
  write(gmdir,'(a)') "/scratch/garel/Galics/GM-RUNS/1024-100h-1Mpc-W5/NewSimu_z0/sf5_zmigr_newFilters2_z2_BlueMuse/"
  
    ! read filter names
  write(filename,'(a,a)') trim(gmdir),'/filters.dat'
  call read_filter_names(filename)

  ! Tibo: wrapup only from ts >=3 because previous ones are useless for cones...
  do ts = 3,nsteps
     ngtot = 0
     do irun = 1,nruns
        numrun = irun2string(irun)
        write(filename,'(a,a,a,a,i3.3)') trim(gmdir),'/',trim(numrun),'/Snapshots/IDsPosVel.',ts
        filename = trim(filename)//char(0)
        !print*,get_ng(filename)
        ngtot = ngtot + get_ng(filename)
     end do
     print*,ts,ngtot
     if (ngtot > 0) then 
        
        allocate(galid(ngtot),haloid(ngtot),inc(ngtot),pos(3,ngtot),vel(3,ngtot))
        iglast = 0
        do irun = 1,nruns
           numrun = irun2string(irun)
           write(filename,'(a,a,a,a,i3.3)') trim(gmdir),'/',trim(numrun),'/Snapshots/IDsPosVel.',ts
           filename = trim(filename)//char(0)
           ng       = get_ng(filename)
           if (ng > 0) then 
              allocate(tgalid(ng),thaloid(ng),tinc(ng),tpos(3,ng),tvel(3,ng))
              call read_momaf_propfile(filename,tGalID,tHaloID,tInc,tPos,tVel)
              galid(iglast+1:iglast+ng)  = tgalid(:)
              haloid(iglast+1:iglast+ng) = thaloid(:)
              inc(iglast+1:iglast+ng)    = tinc(:)
              pos(:,iglast+1:iglast+ng)  = tpos(:,:)
              vel(:,iglast+1:iglast+ng)  = tvel(:,:)
              iglast = iglast + ng
              deallocate(tgalid,thaloid,tinc,tpos,tvel)
           end if
        end do

        write(filename,'(a,a,i3.3)') trim(snapdir),'/IDsPosVel.',ts
        filename = trim(filename)//char(0)
        call write_momaf_propfile(filename,ngtot,galid,haloid,inc,pos,vel)
        deallocate(galid,haloid,inc,pos,vel)

        allocate(m(ngtot),dm(ngtot))
        ! join disc mag files 
        do ifilt=1,nfilters
           iglast = 0
           do irun = 1,nruns
              numrun = irun2string(irun)
              write(filename,'(a,a,a,a,a,a,i3.3)') trim(gmdir),'/',trim(numrun),'/Snapshots/disc_',trim(filters(ifilt)%name),'.',ts
              filename = trim(filename)//char(0)
              ng       = get_ng(filename)
              if (ng > 0) then 
                 allocate(tm(ng),tdm(ng))
                 call read_momaf_magfile(filename,tm,tdm)
                 m(iglast+1:iglast+ng)  = tm(:)
                 dm(iglast+1:iglast+ng) = tdm(:)
                 iglast = iglast + ng
                 deallocate(tm,tdm)
              end if
           end do
           write(filename,'(a,a,a,a,i3.3)') trim(snapdir),'/disc_',trim(filters(ifilt)%name),'.',ts
           filename = trim(filename)//char(0)
           call write_momaf_magfile(filename,m,dm,ngtot)
        end do
        ! join bulge mag files 
        do ifilt=1,nfilters
           iglast = 0
           do irun = 1,nruns
              numrun = irun2string(irun)
              write(filename,'(a,a,a,a,a,a,i3.3)') trim(gmdir),'/',trim(numrun),'/Snapshots/bulge_',trim(filters(ifilt)%name),'.',ts
              filename = trim(filename)//char(0)
              ng       = get_ng(filename)
              if (ng > 0) then 
                 allocate(tm(ng),tdm(ng))
                 call read_momaf_magfile(filename,tm,tdm)
                 m(iglast+1:iglast+ng)  = tm(:)
                 dm(iglast+1:iglast+ng) = tdm(:)
                 iglast = iglast + ng
                 deallocate(tm,tdm)
              end if
           end do
           write(filename,'(a,a,a,a,i3.3)') trim(snapdir),'/bulge_',trim(filters(ifilt)%name),'.',ts
           filename = trim(filename)//char(0)
           call write_momaf_magfile(filename,m,dm,ngtot)
        end do
        ! join burst mag files 
        do ifilt=1,nfilters
           iglast = 0
           do irun = 1,nruns
              numrun = irun2string(irun)
              write(filename,'(a,a,a,a,a,a,i3.3)') trim(gmdir),'/',trim(numrun),'/Snapshots/burst_',trim(filters(ifilt)%name),'.',ts
              filename = trim(filename)//char(0)
              ng       = get_ng(filename)
              if (ng > 0) then 
                 allocate(tm(ng),tdm(ng))
                 call read_momaf_magfile(filename,tm,tdm)
                 m(iglast+1:iglast+ng)  = tm(:)
                 dm(iglast+1:iglast+ng) = tdm(:)
                 iglast = iglast + ng
                 deallocate(tm,tdm)
              end if
           end do
           write(filename,'(a,a,a,a,i3.3)') trim(snapdir),'/burst_',trim(filters(ifilt)%name),'.',ts
           filename = trim(filename)//char(0)
           call write_momaf_magfile(filename,m,dm,ngtot)
        end do
        deallocate(m,dm)
        
     else

        write(filename,'(a,a,i3.3)') trim(snapdir),'/IDsPosVel.',ts
        filename = trim(filename)//char(0)
        call write_momaf_empty_propfile(filename)
        do ifilt = 1, nfilters
           write(filename,'(a,a,a,a,i3.3)') trim(snapdir),'/disc_',trim(filters(ifilt)%name),'.',ts
           filename = trim(filename)//char(0)
           call write_momaf_empty_magfile(filename)
           write(filename,'(a,a,a,a,i3.3)') trim(snapdir),'/bulge_',trim(filters(ifilt)%name),'.',ts
           filename = trim(filename)//char(0)
           call write_momaf_empty_magfile(filename)
           write(filename,'(a,a,a,a,i3.3)') trim(snapdir),'/burst_',trim(filters(ifilt)%name),'.',ts
           filename = trim(filename)//char(0)
           call write_momaf_empty_magfile(filename)
        end do

     end if


  end do
  
  
  stop

end program momaf_wrapup
