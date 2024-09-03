pro read_gal_sfr,outputdir,timestep,sfr

  filename = outputdir + '/gal_sfr.' + string(timestep,'(i3.3)')
  openr,11,filename
  readf,11,ts,ng,dum,dum2
  if ng gt 0 then begin 
     sfr = define_structure('gal_sfr')
     sfr = replicate(sfr,ng)
     readf,11,sfr
  endif else begin 
     print,'No galaxy found ...'
     sfr = -1
  endelse
  close,11

end

