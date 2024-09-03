pro read_halo_gas_results,outputdir,timestep,hgres,nh,renumbering=renumbering

  filename = outputdir + '/halo-gas_results.' + string(timestep,'(i3.3)')
  openr,11,filename
  readf,11,ts,nh,dum
  if nh gt 0 then begin 
     hgres = define_structure('halo_gas_results',renumbering=renumbering)
     hgres = replicate(hgres,nh)
     readf,11,hgres
  endif else begin 
     print,'No halo found ... '
  endelse
  close,11
end

