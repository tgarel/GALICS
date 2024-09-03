pro sdss_mzr, met, mstar, title=title, no_obs_data=no_obs_data

; INPUTS : 
; met : array of metallicities (solar units)
; mstar : array of corresponding stellar masses
; KEYWORDS: 
; - title : title to the plot (e.g. description of sample)
; - no_obs_data : set to not overplot observational results.


jeje_hist2d, alog10(mstar)+11., alog10(met)+8.69, 7.7, 11.7,8.0,9.5,100,100,map,xaxis,yaxis
; NB: 8.69 is Zsun in units of 12+Log(O/H) (footnote from Tremonti, Sec. 3.1.)
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='12+Log(O/H)',xtitle='Log(Mstar)',nodata=0.0,rgb_nodata=[255,255,255],/noerase,title=title

if not keyword_set(no_obs_data) then begin
   mzr = read_tremonti()
   oplot,mzr.mstar,mzr.median,thick=4,col=120,linestyle=2
   oplot,mzr.mstar,mzr.p16,thick=2,col=120,linestyle=2
   oplot,mzr.mstar,mzr.p84,thick=2,col=120,linestyle=2
endif

end
