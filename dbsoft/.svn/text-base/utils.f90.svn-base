module utils

  public

#ifndef MAXPATHSIZE
  integer(kind=4),parameter :: MAXPATHSIZE = 512
#endif
  
  ! Definition of the run 
  character(MAXPATHSIZE) :: gmdir
  integer(kind=4)        :: nsteps, nsubdirs
  logical(kind=4)        :: build_halo_tables, build_galaxy_tables
  logical(kind=4)        :: build_halo_momaf, build_galaxy_momaf

  ! filters used to compute observed magnitudes
  integer(kind=4)             :: nfilters                         ! Total number of filters
  type filter                                                    ! The filter object has several attributes:
     character(20)            :: name                             ! name of filter 
  end type filter
  type (filter),allocatable   :: filters(:)

  ! halo stuff for momaf outputs : 
  integer(kind=4)             :: nhtot
  real(kind=4),allocatable    :: h_pos(:,:),h_vel(:,:)
  integer(kind=8),allocatable :: h_id(:)


  ! IO units
  integer(kind=4),parameter :: errunit       = 6
  ! output units
  integer(kind=4),parameter :: halodbunit    = 15
  integer(kind=4),parameter :: galtabunit    = 16
  integer(kind=4),parameter :: disctabunit   = 17
  integer(kind=4),parameter :: bulgetabunit  = 18
  integer(kind=4),parameter :: bursttabunit  = 19
  integer(kind=4),parameter :: halomomafunit = 20
  integer(kind=4),parameter :: galmomafunit  = 21
  ! input units
  integer(kind=4),parameter :: filters_unit  = 30
  integer(kind=4),parameter :: tmpunit       = 31  
  integer(kind=4),parameter :: hgrunit       = 32
  integer(kind=4),parameter :: hdmunit       = 33
  integer(kind=4),parameter :: hidunit       = 34
  integer(kind=4),parameter :: galidunit     = 35
  integer(kind=4),parameter :: galresunit    = 36
  integer(kind=4),parameter :: instsfrunit   = 37
  integer(kind=4),parameter :: lumunit       = 38
  
contains

!*****************************************************************************************************************

  subroutine get_info

    implicit none 
    
    character(MAXPATHSIZE)  :: file,line,name,value  
    integer(kind=4)         :: i
    
    write(file,'(a)') 'db_params.dat'
    call test_stop(file)
    open(unit=tmpunit,file=file,status='old',form='formatted')
    do
       read(tmpunit,'(a)',end=2) line
       i = scan(line,'=')
       if (i == 0 .or. line(1:1) == '#') cycle
       name  = trim(adjustl(line(:i-1)))
       value = trim(adjustl(line(i+1:)))
       i     = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))
       select case (trim(name))
       case ('gmdir','GMDir')
          gmdir = trim(value)
       case ('halotables','HaloTables')
          read(value,*) build_halo_tables
       case ('galaxytables','GalaxyTables')
          read(value,*) build_galaxy_tables
       case ('halomomaf','HaloMoMaF')
          read(value,*) build_halo_momaf
       case ('galaxymomaf','GalaxyMoMaF')
          read(value,*) build_galaxy_momaf
       case ('nsteps')
          read(value,*) nsteps
       case ('nsubdirs','NSubDirs')
          read(value,*) nsubdirs
       case default 
          write(errunit,*) 'Do not recognize parameter ',trim(name)
       end select
    end do
2   close(tmpunit)
    
    if (.not. build_halo_tables .and. .not. build_galaxy_tables &
         & .and. .not. build_halo_momaf .and. .not. build_galaxy_momaf) stop

    return

  end subroutine get_info
    
!*****************************************************************************************************************
  
  subroutine open_halo_table_files(ts)
    
    implicit none 
    
    integer(kind=4)        :: ts
    character(MAXPATHSIZE) :: file
    
    write(file,'(a,a,i3.3,a)') trim(gmdir),'/halo_table_',ts,'.db'
    open(unit=halodbunit,file=file,status='unknown',form='formatted')
    
    return
    
  end subroutine open_halo_table_files

!*****************************************************************************************************************

  subroutine open_galaxy_table_files(ts)
    
    implicit none 
    
    integer(kind=4)        :: ts
    character(MAXPATHSIZE) :: file
    
    write(file,'(a,a,i3.3,a)') trim(gmdir),'galaxy_table_',ts,'.db'
    open(unit=galtabunit,file=file,status='unknown',form='formatted')
    write(file,'(a,a,i3.3,a)') trim(gmdir),'disc_table_',ts,'.db'
    open(unit=disctabunit,file=file,status='unknown',form='formatted')
    write(file,'(a,a,i3.3,a)') trim(gmdir),'bulge_table_',ts,'.db'
    open(unit=bulgetabunit,file=file,status='unknown',form='formatted')
    write(file,'(a,a,i3.3,a)') trim(gmdir),'burst_table_',ts,'.db'
    open(unit=bursttabunit,file=file,status='unknown',form='formatted')
    
    return

  end subroutine open_galaxy_table_files

!*****************************************************************************************************************

  subroutine  read_filter_names(filename)
    
    ! read the filter names file.
    
    implicit none                       
    
    character(200)  :: filename
    integer(kind=4) :: i,ierr
    
    call test_stop(filename)                                 ! test for files existence.  
    open(unit=filters_unit,status='old',file=filename)
    
    read(filters_unit,*) nfilters
    allocate(filters(nfilters),stat=ierr)
    if (ierr /= 0) stop '> problem allocating filters'
    do i = 1,nfilters
       read (filters_unit,'(a)') filters(i)%name 
    end do
    close (filters_unit)
    
    return
    
  end subroutine read_filter_names
  
!*****************************************************************************************************************

  subroutine test_stop(file)
    
    ! test for file existence
    
    implicit none
    
    character(*)    :: file
    logical(kind=4) :: exist
    
    inquire(file=file,exist=exist)
    if (.not. exist) then
       write(errunit,*) '> fatal error: ',file,'does not exist.'  
       stop
    endif
    
    return
    
  end subroutine test_stop
  
!*****************************************************************************************************************
  
  function irun2string(irun)
    
    implicit none
    
    integer(kind=4) :: irun
    character(3)    :: irun2string

    if (irun < 10) then 
       write(irun2string,'(i1.1)') irun
    elseif (irun < 100) then 
       write(irun2string,'(i2.2)') irun
    elseif (irun < 1000) then 
       write(irun2string,'(i3.3)') irun
    else
       write(errunit,*) 'number of runs too large '
       stop
    end if
    
    return

  end function irun2string

!*****************************************************************************************************************
  
  subroutine write_real4(unit,x,EndLine)

    implicit none
    
    integer(kind=4)          :: unit
    real(kind=4)             :: x
    integer(kind=4),optional :: EndLine
    character(30)            :: value

    write(value,'(e14.6)') x
    if (present(EndLine)) then 
       write(unit,'(a)') trim(adjustl(value))
    else
       write(unit,'(a,",")',advance='no') trim(adjustl(value))
    end if
    
    return

  end subroutine write_real4

!*****************************************************************************************************************
  
  subroutine write_int16(unit,x,EndLine)

    implicit none
    
    integer(kind=4)          :: unit
    integer(kind=8)          :: x
    integer(kind=4),optional :: EndLine
    character(30)            :: value

    write(value,'(i16)') x
    if (present(EndLine)) then 
       write(unit,'(a)') trim(adjustl(value))
    else
       write(unit,'(a,",")',advance='no') trim(adjustl(value))
    end if
    
    return

  end subroutine write_int16

!*****************************************************************************************************************
 
  subroutine write_int8(unit,x,EndLine)

    implicit none
    
    integer(kind=4)          :: unit
    integer(kind=4)          :: x
    integer(kind=4),optional :: EndLine
    character(30)            :: value

    write(value,'(i8)') x
    if (present(EndLine)) then 
       write(unit,'(a)') trim(adjustl(value))
    else
       write(unit,'(a,",")',advance='no') trim(adjustl(value))
    end if
    
    return

  end subroutine write_int8

!*****************************************************************************************************************

  subroutine init_halo_momaf(ts)

    implicit none
    
    integer(kind=4)        :: ts, isubdir, nh
    character(3)           :: numdir
    character(MAXPATHSIZE) :: filename

    ! count total nb of halos and allocate big array.
    nhtot = 0
    do isubdir = 1,nsubdirs
       numdir = irun2string(isubdir) 
       write(filename,'(a,a,a,i3.3)') trim(gmdir),trim(numdir),'/HaloIDs.',ts
       open(unit=hidunit,file=filename,status='old',form='formatted')
       read(hidunit,*) nh
       close(hidunit)
       nhtot = nhtot + nh
    end do
    if (nhtot > 0) then 
       if (allocated(h_id)) then 
          deallocate(h_id,h_pos,h_vel)
       end if
       allocate(h_id(nhtot),h_pos(3,nhtot),h_vel(3,nhtot))
       h_id  = -1
       h_pos = 0.0
       h_vel = 0.0
    end if

    return

  end subroutine init_halo_momaf

!*****************************************************************************************************************
  
  subroutine dump_halo_momaf_file(ts)
    
    implicit none
    
    integer(kind=4) :: ts
    character(MAXPATHSIZE) :: file
    
    write(file,'(a,a,i3.3,a)') trim(gmdir),'/Halo_IDsPosVel.',ts
    file = trim(file)//char(0)  ! ensures compatibility with C
    if (nhtot > 0) then 
       call write_halo_momaf_file(file,nhtot,h_id,h_pos,h_vel)
    else
       call write_empty_halo_momaf_file(file)
    end if
    
    return

  end subroutine dump_halo_momaf_file

!*****************************************************************************************************************

end module utils
