pro read_halo_reservoir,outputdir,timestep,hreserv

  filename = outputdir + '/halo-reservoir.'+string(timestep,format='(i3.3)')
  openr,11,filename
  readf,11,ts,nh,dum
  if nh gt 0 then begin 
     hreserv = define_structure('halo_reservoir')
     hreserv = replicate(hreserv,nh)
     readf,11,hreserv
  endif else begin 
     print,'No halo found ... '
  endelse
  close,11

end
