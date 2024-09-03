pro sdss_cmr,mags,volume,hubble,no_obs_data=no_obs_data,pop1=pop1,pop2=pop2

; INPUTS: 
; - mags: structure array containing at least SDSS_r and SDSS_u tags.
; - volume : volume filled up by the galaxies of mags, in Mpc^3
; - hubble : in units of 100km/s/Mpc
; KEYWORDS : 
; - no_obs_data : if set, do not overplot observational data on top of simulation.
; - pop1, pop2  : list of indices (to mags) of sub-populations to
;   overplot.

cmin  = 0.0
cmax  = 3.5
cbins = 100

for imag = 0,15 do begin 

   mag = -23.25 + imag * 0.5 ; bins from Baldry et al 04
   if not keyword_set(no_obs_data) then begin 
      bcd = baldry_coldist(mag,cmin,cmax,cbins)
      plot,bcd.col,bcd.nred+bcd.nblue,linestyle=2,xtitle='u-r',ytitle='# / Mpc!u3!n / d(u-r)',title='Mr = '+strtrim(mag,2),charsize=2
   endif
   ii  = where(mags.sdss_r lt mag+0.25 and mags.sdss_r ge mag-0.25,ni)
   if ni ne 0 then begin 
      c = mags(ii).sdss_u - mags(ii).sdss_r
      jeje_histo,c,cmin,cmax,cbins,hc
      if not keyword_set(no_obs_data) then begin 
         oplot,hc.x,hc.dn/hc.dx/volume,psym=10,thick=3
      endif else begin 
         plot,hc.x,hc.dn/hc.dx/volume,psym=10,thick=3,xtitle='u-r',ytitle='# / Mpc!u3!n / d(u-r)',title='Mr = '+strtrim(mag,2),charsize=2
      endelse
      
      if keyword_set(pop1) then begin 
         ii = where(mags(pop1).sdss_r lt mag+0.25 and mags(pop1).sdss_r ge mag-0.25,ni)
         if ni ne 0 then begin 
            c  = mags(pop1(ii)).sdss_u - mags(pop1(ii)).sdss_r
            jeje_histo,c,cmin,cmax,cbins,hc
            oplot,hc.x,hc.dn/hc.dx/volume,psym=10,thick=3,linestyle=2,col=150
         endif
      endif

      if keyword_set(pop2) then begin 
         ii = where(mags(pop2).sdss_r lt mag+0.25 and mags(pop2).sdss_r ge mag-0.25,ni)
         if ni ne 0 then begin 
            c  = mags(pop2(ii)).sdss_u - mags(pop2(ii)).sdss_r
            jeje_histo,c,cmin,cmax,cbins,hc
            oplot,hc.x,hc.dn/hc.dx/volume,psym=10,thick=3,linestyle=3,col=220
         endif
      endif

   endif
endfor


end
