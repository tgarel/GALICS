pro read_galaxyids,outputdir,timestep,gtree,ng

  filename = outputdir + '/GalaxyIDs.'+string(timestep,format='(i3.3)')
  openr,11,filename
  ng = 0L
  readf,11,ng
  if ng gt 0 then begin 
     treestruct = define_structure('galaxyids') 
     gtree = replicate(treestruct,ng)
     readf,11,gtree
  endif else begin 
     print,'No galaxy found ... '
     gtree = -1
  endelse
  close,11

end
