pro sdss_hod,halo,gals,mags

; code to reproduce some plots from Yang, Mo and Van den Bosch (2008).
; NB: these authors define the "central" galaxy as the brightest (or
; most massive). I reproduce this definition here.

; halo mass sampling (in log)
lmhmin  = 11
lmhmax  = 15
lmhbins = 40

; magnitude limits
maglim = [-21.5,-21.0,-20.5,-20.0,-19.5,-19.0,-18.5,-18.0]
nmags  = n_elements(maglim)

; 1. identify central/satellite galaxies in each halo (and build
; halo->galaxy link).
; --------------------------------------------------------------
galinh = replicate({firstgal:0L,lastgal:0L},n_elements(halo))
for ih = 0L, n_elements(halo)-1L do begin 
   
endfor


; count galaxies in halos
nh        = n_elements(halo)
hod       = fltarr(nh,nmags)
;brightest = fltarr(nh)
;ibright   = lonarr(nh)
sdssr     = mags.sdss_r
for imag = 0,nmags-1 do begin 
   selection = where(sdssr le maglim(imag),ng)
   for ig = 0L,ng-1L do begin 
      ih = gals(selection(ig)).hno - 1L
      hod(ih,imag) = hod(ih,imag) + 1.0
;;       if brightest(ih) gt sdssr(selection(ig)) then begin  ; flag brightest (i.e. central) gal
;;          brightest(ih) = sdssr(selection(ig))
;;          ibright(ih)   = ig
;;       endif
   endfor
endfor

stop

; build HOD for each magnitude limit (i.e. <Ngal> as a function of Mh)
hcnt = ( gcnt = fltarr(lmhbins))
for imag = 0,n_elements(maglim)-1 do begin 
   dx = (lmhmax - lmhmin) / float(lmhbins)
   x  = findgen(lmhbins) * dx + lmhmin + 0.5 * dx
   imh = fix((alog10(halo.mvir)+11.0 - lmhmin)/dx)
   ii = where(imh ge 0 and imh lt lmhbins,ni)
   if ni ne 0 then begin 
      for ih = 0L, ni -1L do begin 
         hcnt(imh(ii(ih))) ++ 
         gcnt(imh(ii(ih))) = gcnt(imh(ii(ih))) + hod(ii(ih),imag)
      endfor
   endif
   ii = where(hcnt gt 0,ni)
   if ni ne 0 then gcnt(ii) = gcnt(ii) / hcnt(ii)
   plot_io,x,gcnt,psym=10,xr=[0.1,1000.],/xs,xtitle='M!dvir!n',ytitle='<Ng>'
endfor

end
