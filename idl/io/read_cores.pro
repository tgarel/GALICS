pro read_cores,outputdir,timestep,cores
  
  filename = outputdir + '/cores.' + string(timestep,format='(i3.3)')
  openr,11,filename
  readf,11,ts,ng,aexp,age_univ
  if ng gt 0 then begin 
     cores = define_structure('cores')
     cores = replicate(cores,ng)
     readf,11,cores
  endif else begin 
     print,'No galaxy found ...'
  endelse
  close,11

end
