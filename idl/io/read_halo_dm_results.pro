pro read_halo_dm_results,outputdir,timestep,hdmres

  filename = outputdir + '/halo_dm_results.' + string(timestep,'(i3.3)')
  openr,11,filename
  readf,11,ts,nh,dum
  if nh gt 0 then begin 
     hdmres = define_structure('halo_dm_results')
     hdmres = replicate(hdmres,nh)
     readf,11,hdmres
  endif else begin 
     print,'No halo found ...'
  endelse
  close,11

end
