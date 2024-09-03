pro read_gal_ejecta, outputdir, timestep, galej,galeji

  file = outputdir + 'gal_ejecta.'+string(timestep,format='(i3.3)')
  openr,11,file
  
  
  readf,11,format='(i3,1x,i6,1x,e14.6,1x,e14.6)',st,n_total,aexp,age_univ
  galeji = {st:st,ng:n_total,aexp:aexp,age_univ:age_univ}
  print,n_total

  if n_total gt 0 then begin

     nsub  = 0L
     mej   = dblarr(200)
     mej(*)= -99.0d0
     galej = {mej:mej}
     galej = replicate(galej,n_total)

     for igal = 0L, n_total - 1L do begin
        readf,11,format='(i3)',nsub
        ej = dblarr(nsub)
        readf,11,format='(nsub(e14.6))',ej
        galej(igal).mej = ej
     endfor
  endif else begin
     gal = -1
     ;print,'no galaxy'
  endelse
  
  close,11

end
