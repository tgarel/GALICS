pro read_gal_uv1200, outputdir, timestep, gal,galinfo

  log_L_sun = 22.583d0 ; solar luminosity in units of 10^11 Msun

  file = outputdir + 'uv1200_spec.'+string(timestep,format='(i3.3)')
  openr,11,file

; read header 
  st = (n_total = 0L) & aexp = ( age_univ = 0.0) 
  readf,11,format='(i3,1x,i6,1x,e14.6,1x,e14.6,1x,e14.6)',st,n_total,aexp,age_univ
  
  galinfo = {st:st,n_total:n_total,aexp:aexp,age_univ:age_univ}
  if n_total gt 0 then begin

     uv1200      = dblarr(21)
     mean        = 0.0d0
     gal         = {mean:mean,uv1200:uv1200}
     gal         = replicate(gal,n_total)
     
     uv1200(*) = 0.0d0
     
     for igal = 0L, n_total - 1L do begin
        readf,11,format='(22(e14.6,1x))',mean,uv1200
        
        gal(igal).mean        = mean * 10.0d0^log_L_sun * 10.0d0^11 / 10.0d0^4 ;in erg/s/A 
        gal(igal).uv1200      = uv1200  * 10.0d0^log_L_sun * 10.0d0^11 / 10.0d0^4 ;in erg/s/A 
                                ;gal(igal).mean        = mean * 10.0d0^log_L_sun * 10.0d0^11 / 10.0d0^4  ;in erg/s/A 

     endfor
     
  endif else begin
     gal = -1
     print,'no galaxy'
  endelse

  close,11

end
