pro read_bol_lum,outputdir,timestep,lum

  filename = outputdir + '/bol_lum.' + string(timestep,format='(i3.3)')
  openr,11,filename
  readf,11,dum,ng
  if ng gt 0 then begin 
     bol = fltarr(2,ng,/nozero)
     readf,11,bol
  endif else begin 
     print,'No galaxy found ...'
  endelse 
  close,11

  if ng gt 0 then begin
     lum = create_struct('bol', 0.0, 'ir', 0.0)
     lum = replicate(lum,ng)
     lum.bol = reform(bol(0,*))
     lum.ir = reform(bol(1,*))
  endif
end


