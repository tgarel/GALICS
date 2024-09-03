pro analyse,hgr,run=run,psfile=psfile,noread=noread,datadir=datadir

loadct,39
if keyword_set(psfile) then begin 
   !p.multi=[0,2,2]
   set_plot,'ps'
   device,filename=psfile,/color,bits_per_pixel=8,xsize=20,ysize=20
endif else begin 
   device,decomposed=0
endelse

;; if not keyword_set(run) then begin 
;;    run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
;;    run.outputdir = '/Users/blaizot/Documents/ASTRO/data/GM2/'
;;    run.snapzdir  = '/Users/blaizot/Documents/ASTRO/data/GM2/1/Snapshots/'
;;    run.dirmin    = 1
;;    run.dirmax    = 3
;;    run.name      = '512!u3!n std'
;;    run.volumempc = 142.85714^3 /60.*3.
;;    run.hubble    = 0.7
;; endif

if not keyword_set(run) then begin 
   run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
   run.outputdir='/data/blaizot/GM-RUNS/256-100h-1Mpc/WSFFB4/'
   run.snapzdir='/data/blaizot/GM-RUNS/256-100h-1Mpc/WSFFB4/1/'
   run.dirmin=1
   run.dirmax=24
   run.name='WSFFB'
   run.volumempc = 142.85714^3.
   run.hubble=0.7
endif

; ------------------------------------------
; baryonic content of halos at z=0
; ------------------------------------------
z = 0.
timestep = z2ts(z,run.snapzdir)
if not keyword_set(noread) then begin
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + string(i,format='(i1.1)')+'/'
   if i eq run.dirmin then begin
      read_halo_gas_results,outputdir,timestep,hgr,/renumbering
   endif else begin 
      read_halo_gas_results,outputdir,timestep,hgrt,/renumbering
      hgr  = [hgr,hgrt]
   endelse
endfor
endif

;; plot_oo,hgr.mfof*1.d11,(hgr.mcoldgaz+hgr.mhotgaz)*1.d11,psym=3,xtitle='M!dFoF!n (M!d!9n!3!n)',$
;;       ytitle='M!dbar!n (M!d!9n!3!n)',yr=[1.e7,1.e15],/ys,charsize=1.5
;; oplot,hgr.mfof*1.d11,hgr.mcoldgaz*1.d11,psym=3,col=200
;; oplot,hgr.mfof*1.d11,hgr.mhotgaz*1.d11,psym=3,col=240
;; oplot,hgr.mfof*1.d11,hgr.mfof*1.d11*0.044/0.3 ; baryon fraction 

;; plot_oi,hgr.mvir*1.d11,hgr.mcoldgaz/(hgr.mcoldgaz+hgr.mhotgaz),psym=3,yr=[-0.1,1.1],/ys
;; ii = where(hgr.mcoldgaz eq 0.0,ni)
;; if ni ne 0 then begin 
;;    cfrac = hgr(ii).mcoldgaz - 0.05
;;    oplot,hgr(ii).mvir*1.d11,cfrac,psym=3
;; endif
;; oplot,hgr.mvir*1.d11,hgr.mhotgaz/(hgr.mcoldgaz+hgr.mhotgaz),psym=3,col=245

; fraction of baryons in the hot phase
x = alog10(hgr.mvir)+11.0
y = hgr.mhotgaz/(hgr.mcoldgaz+hgr.mhotgaz)
jeje_hist2d,x,y,10,15,0.0,1.01,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='fraction of baryons in hot gas',$
          xtitle='Log!d10!n(M!dVir!n) [M!d!9n!n!3]',nodata=0.0,rgb_nodata=[255,255,255],/noerase

; metallicity of the hot phase
x = alog10(hgr.mvir)+11.0
y = hgr.mhotz/hgr.mhotgaz/0.016
jeje_hist2d,x,y,10,15,0.0,0.6,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='Z/Zsun (hot gas)',$
          xtitle='Log!d10!n(M!dVir!n) [M!d!9n!n!3]',nodata=0.0,rgb_nodata=[255,255,255],/noerase


; ------------------------------------------
;         z = 0.1 : SDSS comparison
; ------------------------------------------
z = 0.1
if not keyword_set(noread) then begin 
timestep = z2ts(z,run.snapzdir)
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir +$
             string(i,format='(i1.1)')+'/'
   if i eq run.dirmin then begin
      read_gal_sfr, outputdir, timestep, sfr
      read_gal_results,outputdir,timestep,gr
      read_halo_gas_results,outputdir,timestep,hgr,/renumbering
      read_gal_mags,outputdir,timestep,bumask,filters,mags,/total
   endif else begin 
      read_gal_sfr, outputdir, timestep, sfrt
      read_gal_results,outputdir,timestep,grt
      read_halo_gas_results,outputdir,timestep,hgrt,/renumbering
      read_gal_mags,outputdir,timestep,bumask,filters,magst,/total
      sfr  = [sfr,sfrt]
      gr   = [gr,grt]
      hgr  = [hgr,hgrt]
      mags = [mags,magst]
   endelse
endfor
sfrtot = 100.*(sfr.disc_sfr+sfr.bulge_sfr+sfr.burst_sfr)
endif


; stellar mass function 
; ---------------------
mstar = gr.disc_minstar + gr.bulge_minstar + gr.burst_minstar
lmstar = alog10(mstar) + 11.0 ; in Msun 
jeje_histo, lmstar, 7., 14., 100, hmstar
plot_io,hmstar.x,hmstar.dn/hmstar.dx/run.volumempc,psym=10,xtitle='Log!d10!n M!dstar!n [M!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nM!dstar!n / Mpc!u3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,0.1],/ys
plot_obs_data,'BDP_DR3MF',run.hubble


; disc specific star formation rate vs. disc size.
;plot_oo,gr.disc_rgal*1000.,sfr.disc_sfr/gr.disc_mgal,psym=3,xr=[1.e-4,1.e-1]*1000.,/xs,xtitle='disc radius (kpc)',charsize=1.5,ytitle='SSFR (Gyr!u-1!n)'


; SDSS luminosity functions
; -------------------------
sdss_lfs,mags,run.volumempc,run.hubble


; Gas/star fractions 
; ------------------
;; plot_io,lmstar,mstar/(gr.disc_mgal+gr.bulge_mgal+gr.burst_mgal),psym=3,xtitle='Log!d10!n(M!dstar!n) [M!d!9n!n!3]',$
;;         ytitle='f!dstar!n',yr=[0.01,1.1],/ys,charsize=1.5,xr=[7,11.5],/xs
sel = where(mstar gt 1.e-4); and (gr.burst_mstar + gr.bulge_mstar)/mstar lt 0.4)
;lfrac = lmstar(sel) - alog10(gr(sel).disc_mgal+gr(sel).bulge_mgal+gr(sel).burst_mgal) - 11.0
lfrac = mstar(sel)/(gr(sel).disc_mgal+gr(sel).bulge_mgal+gr(sel).burst_mgal)
jeje_hist2d,lmstar(sel),1.-lfrac,8.,11.5,0.0,1.01,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='f!dgas!n',xtitle='Log!d10!n(M!dstar!n) [M!d!9n!n!3]',rgb_nodata=[255,255,255],nodata=0.0,/noerase


; Gas metallicity (vs. mstar)
; ---------------------------
sel = where(mstar gt 1.e-4); and (gr.burst_mstar + gr.bulge_mstar)/mstar lt 0.4)
met = (gr(sel).disc_mcoldz + gr(sel).bulge_mcoldz + gr(sel).burst_mcoldz)/(gr(sel).disc_mcold + gr(sel).bulge_mcold + gr(sel).burst_mcold)/0.016
jeje_hist2d,lmstar(sel),alog10(met),8.,11.5,-1.5,0.5,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='Z/Zsun',xtitle='Log!d10!n(M!dstar!n) [M!d!9n!n!3]',rgb_nodata=[255,255,255],nodata=0.0,/noerase

; color magnitude relation 
; ------------------------
x = mags.sdss_r
y = mags.sdss_u - mags.sdss_r
jeje_hist2d,x,y,-23,-14,1,3,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='u-r',xtitle='r',rgb_nodata=[255,255,255],nodata=0.0,/noerase




if keyword_set(psfile) then begin 
   device,/close
   set_plot,'x'
endif

end

