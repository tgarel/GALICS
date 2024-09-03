pro sdss_sizelum,galres,mags,cbd,no_obs_data=no_obs_data

; compares galics data with results from Shen et al. (2003)
; NB: (1) Shen et al use a concentration parameter (R90/R50) to split their
; population into early and late types (at value 2.86). With this cut they
; find 32% early types (ellipticals). I reproduce this selection here.
; (2) Shen et al give a single radius for the whole galaxy, which is R50
; (Petrosian). I do not take into account possible biases (under-estimate of
; bulges?) due to Petrosian technique. I compute a R50 radius integrating over
; the composite profile of each galaxy.

; KEYWORDS: 
; no_obs_data : set true not to overplot observations

; save concentration and bulge to total into cbd
ii = where(mags.sdss_r gt -24 and mags.sdss_r le -16.,ni)
cbd = {c:dblarr(ni),b2t:dblarr(ni)}

; define min/max radii and number of bins for size histograms
histormin  = 0.1
histormax  = 50.
historbins = 100

; loop on magnitude bins
for magbin = 0,15 do begin 
   magmin = -24. + magbin * 0.5
   magmax = magmin + 0.5
   mag = 0.5*(magmax + magmin)
   magsel = where(mags.sdss_r lt magmax and mags.sdss_r gt magmin,ni)
   if ni ne 0 then begin 
      selr50 = ( selr90 = dblarr(ni))
      ; retrieve scale-lengths and luminosities of each component of each
      ; selected galaxy
      rd  = galres(magsel).disc_rgal  * 1000.d0 ; Mpc->kpc
      rbg = galres(magsel).bulge_rgal * 1000.d0
      rbs = galres(magsel).burst_rgal * 1000.d0
      ld  = mags(magsel).disc_sdss_r 
      ld  = 10.d0^(-0.4*ld)
      lbg = mags(magsel).bulge_sdss_r
      lbg = 10.d0^(-0.4*lbg)
      lbs = mags(magsel).burst_sdss_r
      lbs = 10.d0^(-0.4*lbs)
      for igal = 0L,ni-1L do begin 
         ; sample profiles from r=0 to at least 4*rd or 16*rb so as to reach R90
         ; in all cases.
         rmin  = 0.0d0
         rmax  = max([4.d0*rd(igal),16.0d0*rbg(igal),16.d0*rbs(igal)])
         rbins = 1000
         if rd(igal) gt 0  then dprof = expdisc(rd(igal),rmin,rmax,rbins)    $
         else dprof = {r:(dindgen(rbins)+0.5)/float(rbins)*(rmax-rmin)+rmin,s:dblarr(rbins)}
         if rbg(igal) gt 0 then gprof = hernquist(rbg(igal),rmin,rmax,rbins) $
         else gprof = {r:(dindgen(rbins)+0.5)/float(rbins)*(rmax-rmin)+rmin,s:dblarr(rbins)}
         if rbs(igal) gt 0 then sprof = hernquist(rbs(igal),rmin,rmax,rbins) $
         else sprof = {r:(dindgen(rbins)+0.5)/float(rbins)*(rmax-rmin)+rmin,s:dblarr(rbins)}
         ; build composite profile (luminosity weighted)
         cprof = ld(igal)*dprof.s + lbg(igal)*gprof.s + lbs(igal)*sprof.s
         ; normalise to total luminosity
         cprof = cprof / (ld(igal) + lbg(igal) + lbs(igal))
         ; get R50 and R90 
         i50 = (where(cprof gt 0.5))(0)
         i90 = (where(cprof gt 0.9))(0)
         if i50 lt 0 or i90 lt 0 then stop,'R90 not found in sdss_sizelum'
         selr50(igal) = dprof.r(i50)
         selr90(igal) = dprof.r(i90)
      endfor
      ; define concentration index and early/late populations
      concentration = selr90 / selr50
      cbd(magsel).c   = concentration
      cbd(magsel).b2t = (lbg+lbs)/(lbg+lbs+ld)
      early = where(concentration ge 2.86,nearly)
      late  = where(concentration lt 2.86,nlate)
      print,'in mag bin '+strtrim((magmax+magmin)*0.5,2)+' :'
      print,'nearly = '+strtrim(nearly,2)+', and nlate = '+strtrim(nlate,2)
      if not keyword_set(no_obs_data) then begin 
         ; get distributions from Shen et al. (2003)
         shen = shen_sizelum(mag,histormin,histormax,historbins)
         ; compute distributions and plot the stuff (normalised as Shen : area in Log10(r) space is one)
         if nearly gt 0 then begin 
            jeje_histo,alog10(selr50(early)),alog10(histormin),alog10(histormax),historbins,hearly
            plot_oi,shen.x,shen.netg,linestyle=2,thick=2,xr=[0.3,40.],/xs,yr=[0,10],/ys, $
                    xtitle='R50 [kpc]',title='Early, Mr='+strtrim(mag,2)
            oplot,10.0d0^hearly.x,hearly.dn/hearly.dx/hearly.norm,psym=10,thick=3
         endif
         if nlate  gt 0 then begin 
            jeje_histo,alog10(selr50(late)),alog10(histormin),alog10(histormax),historbins,hlate
            plot_oi,shen.x,shen.nltg,linestyle=2,thick=2,xr=[0.3,40.],/xs,yr=[0,10],/ys, $
                    xtitle='R50 [kpc]',title='Late, Mr='+strtrim(mag,2)
            oplot,10.0d0^hlate.x,hlate.dn/hlate.dx/hlate.norm,psym=10,thick=3
         endif
      endif else begin 
         if nearly gt 0 then begin 
            jeje_histo,alog10(selr50(early)),alog10(histormin),alog10(histormax),historbins,hearly
            plot_oi,10.0d0^hearly.x,hearly.dn/hearly.dx/hearly.norm,psym=10,thick=3,xr=[0.3,40.],/xs,$
                    yr=[0,10],/ys,xtitle='R50 [kpc]',title='Early, Mr='+strtrim(mag,2)
         endif 
         if nlate  gt 0 then begin 
            jeje_histo,alog10(selr50(late)),alog10(histormin),alog10(histormax),historbins,hlate
            plot_oi,10.0d0^hlate.x,hlate.dn/hlate.dx/hlate.norm,psym=10,thick=3,$
                    xr=[0.3,40.],/xs,yr=[0,10],/ys,xtitle='R50 [kpc]',title='Late, Mr='+strtrim(mag,2)
         endif
      endelse
   endif                        ; if gals in mag bin
endfor ; end loop on magbins


end
