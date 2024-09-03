pro read_haloids,outputdir,timestep,htree

  filename = outputdir + '/HaloIDs.'+string(timestep,format='(i3.3)')
  openr,11,filename
  readf,11,nh
  if nh gt 0 then begin 
     treestruct = define_structure('haloids') 
     htree = replicate(treestruct,nh)
     readf,11,htree
  endif else begin 
     print,'No galaxy found ... '
  endelse
  close,11

end
