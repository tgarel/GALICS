pro calibration_plots, run, psfile,no_obs_data=no_obs_data,high_z=high_z

; INPUT : 
; - run : a structure describing the GalaxyMaker run to make plots with.
;  -> contains at least the following elements : 
;     - outputdir : where the sub-folders are 
;     - snapzdir  : where the file snaps_redshifts.dat is
;     - dirmin    : first subdir to use
;     - dirmax    : last subdir to use (NB: all dirs between dirmin and
;        dirmax will be read)
;     - name      : a string describing the run (sometimes used for plot
;       titles)
;     - volumempc : the comoving volume to analyse, in Mpc^3 (NB: can be only
;       a fraction of the whole box depending on dirmin and dirmax)
;     - hubble    : the reduced hubble constant (in units of 100 km/s/Mpc)
; - psfile : the name of the postscript file with the plots.
; KEYWORDS : 
; - no_obs_data : if set, do not overplot observational data on top of simulation.

set_plot,'ps'
device,filename=psfile,/color,bits_per_pixel=8,xsize=18,ysize=24,yoffset=3,xoffset=1
!p.multi=[0,3,4]
loadct,39

; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; I.   Low-redshift Universe : SDSS 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; define z=0 snapshot
z        = 0.0
timestep = z2ts(z,run.snapzdir)
; read useful information 
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + strtrim(i,2) + '/'
   if i eq run.dirmin then begin
      read_gal_mags,outputdir,timestep,filters,mags,/total
      read_gal_results,outputdir,timestep,galres,grinfo
      read_gal_sfr,outputdir,timestep,sfr
      read_bol_lum,outputdir,timestep,bol
   endif else begin 
      read_gal_mags,outputdir,timestep,filters,magst,/total
      mags = [mags,magst]
      read_gal_results,outputdir,timestep,galrest,grinfo
      galres = [galres,galrest]
      read_gal_sfr,outputdir,timestep,sfrt
      sfr = [sfr,sfrt]
      read_bol_lum,outputdir,timestep,bolt
      bol = [bol,bolt]
   endelse
endfor
; read halo properties too
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + strtrim(i,2)+'/'
   if i eq run.dirmin then begin
      read_halo_gas_results,outputdir,timestep,hgr,/renumbering
   endif else begin 
      read_halo_gas_results,outputdir,timestep,hgrt,/renumbering
      hgr  = [hgr,hgrt]
   endelse
endfor


; I.1/ Luminosity functions at z = 0 (compare to Blanton et al., 2005)
; --------------------------------------------------------------------
sdss_lfs,mags,run.volumempc,run.hubble,no_obs_data=no_obs_data

; I.2/ Stellar Mass Function (compare to Panter et al., 2007)
; -----------------------------------------------------------
mstar  = galres.disc_minstar + galres.bulge_minstar + galres.burst_minstar
lmstar = alog10(mstar+1.e-5) + 11.0 ; in Msun 
jeje_histo, lmstar, 7., 14., 100, hmstar
plot_io,hmstar.x,hmstar.dn/hmstar.dx/run.volumempc,psym=10,xtitle='Log!d10!n M!dstar!n [M!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nM!dstar!n / Mpc!u3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,0.1],/ys
if not keyword_set(no_obs_data) then plot_obs_data,'BDP_DR3MF',run.hubble

; I.3/ Color-magnitude relation (compare to Baldry et al., 2004)
; --------------------------------------------------------------
sats = where(galres.r_orbit gt 0.0,ns)
cent = where(galres.r_orbit eq 0.0,nc)
sdss_cmr,mags,run.volumempc,run.hubble,no_obs_data=no_obs_data,pop1=sats,pop2=cent

; I.4/ Stellar-Mass / Gas-Metallicity relation (compare to Tremonti et
; al.,2004)
; --------------------------------------------------------------------
; total metallicity (average over three components, no SFR selection)
gasmet = galres.disc_mcoldz + galres.bulge_mcoldz + galres.burst_mcoldz
gas    = galres.disc_mcold  + galres.bulge_mcold  + galres.burst_mcold
met    = gasmet / gas / 0.016 ; metallicity in solar units
sdss_mzr,met, mstar,title='total (disc+bulge+burst) - no SFR selection',no_obs_data=no_obs_data
; disc metallicity
sel = where(sfr.disc_sfr10 gt 0.,nsel)
print,nsel
if nsel gt 0 then begin 
   gasmet = galres(sel).disc_mcoldz
   gas    = galres(sel).disc_mcold
   met    = gasmet / gas / 0.016 ; metallicity in solar units
   sdss_mzr,met, galres(sel).disc_minstar,title='disc - SFR10 > 0',no_obs_data=no_obs_data
endif
gasmet = galres.disc_mcoldz
gas    = galres.disc_mcold
met    = gasmet / gas / 0.016   ; metallicity in solar units
sdss_mzr,met, galres.disc_minstar,title='disc - all',no_obs_data=no_obs_data

; bulge metallicity
sel = where(sfr.bulge_sfr10 gt 0.,nsel)
print,nsel
if nsel gt 0 then begin 
   gasmet = galres(sel).bulge_mcoldz
   gas    = galres(sel).bulge_mcold 
   met    = gasmet / gas / 0.016 ; metallicity in solar units
   sdss_mzr,met, galres(sel).bulge_minstar,title='bulge - SFR10 > 0',no_obs_data=no_obs_data
endif
gasmet = galres.bulge_mcoldz
gas    = galres.bulge_mcold 
met    = gasmet / gas / 0.016 ; metallicity in solar units
sdss_mzr,met, galres.bulge_minstar,title='bulge' ,no_obs_data=no_obs_data

; burst metallicity
sel = where(sfr.burst_sfr10 gt 0.,nsel)
print,nsel
if nsel gt 0 then begin 
   gasmet = galres(sel).burst_mcoldz
   gas    = galres(sel).burst_mcold
   met    = gasmet / gas / 0.016 ; metallicity in solar units
   sdss_mzr,met, galres(sel).burst_minstar,title='burst - SFR10 > 0',no_obs_data=no_obs_data
endif
gasmet = galres.burst_mcoldz
gas    = galres.burst_mcold
met    = gasmet / gas / 0.016 ; metallicity in solar units
sdss_mzr,met, galres.burst_minstar,title='burst',no_obs_data=no_obs_data

; I.5/ Size-luminosity distribution (compare to Shen et al., 2003)
; ----------------------------------------------------------------
sdss_sizelum,galres,mags,cbd,no_obs_data=no_obs_data
; plot concentration versus B/T 
jeje_hist2d,cbd.c,cbd.b2t,2.,10.,0.,1.,100,100,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,xtitle='C (R90/R50)',ytitle='B/T (r-band)',nodata=0,rgb_nodata=[255,255,255],/noerase
oplot,[2.86,2.86],[0,1],linestyle=2

; I.6/ Gas content of galaxies (compare to Zwaan et al., 2003)
; ------------------------------------------------------------
; HI mass function 
b = 3.1/8.1                     ; M_HII/M_HI, cosmic measurements
a = 0.76                        ; fraction of hydrogen
disc_HI  =  galres.disc_mcold - galres.disc_mcoldz ;
disc_HI  = a*disc_HI                               ;
disc_HI  = disc_HI/(1.0 + b)                       ;
disc_HI  = alog10(disc_HI + 2.e-10) + 11.0          ; (2.e-10 to avoid zero...)
jeje_histo,disc_hi,7.,11.,21,h
plot,h.x,alog10(h.dn/h.dx/run.volumempc),psym=10,thick=3,xtitle='Log!d10!n M!dHI!n (M!d!9n!n!3)',$
     ytitle="Log!d10!n [dN / dLog M / Mpc!u3!n]"
plot_obs_data,'HI_MF',run.hubble

; gas fraction in each component 
y = galres.disc_mcold / (galres.disc_mcold + galres.disc_minstar)
x = alog10(galres.disc_minstar) + 11.0
jeje_hist2d,x,y,6,12,0.0,1.01,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='f!dgas!n', title='disc',$
          xtitle='Log!d10!n(M!dstar!n) [M!d!9n!n!3]',nodata=0.0,rgb_nodata=[255,255,255],/noerase
y = galres.bulge_mcold / (galres.bulge_mcold + galres.bulge_minstar)
x = alog10(galres.bulge_minstar) + 11.0
jeje_hist2d,x,y,6,12,0.0,1.01,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='f!dgas!n', title='bulge',$
          xtitle='Log!d10!n(M!dstar!n) [M!d!9n!n!3]',nodata=0.0,rgb_nodata=[255,255,255],/noerase
y = galres.burst_mcold / (galres.burst_mcold + galres.burst_minstar)
x = alog10(galres.burst_minstar) + 11.0
jeje_hist2d,x,y,6,12,0.0,1.01,50,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='f!dgas!n', title='burst',$
          xtitle='Log!d10!n(M!dstar!n) [M!d!9n!n!3]',nodata=0.0,rgb_nodata=[255,255,255],/noerase

; I.7/ Halo Occupation Distribution (compare to Yang, Mo & Van den
; Vosch 2008)
; ----------------------------------------------------------------
sdss_hod,hgr,galres,mags


; I.10/ Halo gas content 
; ----------------------
; fraction of baryons in the hot phase
x = alog10(hgr.mvir)+11.0
y = hgr.mhotgaz/(hgr.mcoldgaz+hgr.mhotgaz)
jeje_hist2d,x,y,8,15,0.0,1.01,70,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='fraction of baryons in hot gas',$
          xtitle='Log!d10!n(M!dVir!n) [M!d!9n!n!3]',nodata=0.0,rgb_nodata=[255,255,255],/noerase

; metallicity of the hot phase
x = alog10(hgr.mvir)+11.0
y = hgr.mhotz/hgr.mhotgaz/0.016
jeje_hist2d,x,y,8,15,0.0,0.6,70,50,map,xaxis,yaxis
tvim_true,map,xrange=xaxis,yrange=yaxis,ytitle='Z/Zsun (hot gas)',$
          xtitle='Log!d10!n(M!dVir!n) [M!d!9n!n!3]',nodata=0.0,rgb_nodata=[255,255,255],/noerase


; I.11/ Properties of low-z disc galaxies
; ----------------------------------------
; select discs
discs = where(galres.disc_minstar gt 2.*(galres.bulge_minstar + galres.burst_minstar),ndiscs)
if ndiscs ne 0 then begin 
; plot gas vs. mass
gasfrac = galres(disc).disc_mcold / galres(disc).disc_mgal
plot_oo,galres(disc).disc_minstar*1.d11,gasfrac,xr=[1.d8,1.d12],yr=[0.01,1],/xs,/ys

endif ; ndiss ne 0

; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; II.   High-redshift Universe 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; define z=2.5 snapshot
z        = 2.5
timestep = z2ts(z,run.snapzdir)
; read useful information 
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + strtrim(i,2) + '/'
   if i eq run.dirmin then begin
      read_gal_sfr,outputdir,timestep,sfr
   endif else begin 
      read_gal_sfr,outputdir,timestep,sfrt
      sfr = [sfr,sfrt]
   endelse
endfor
; distribution of SFRs 
sfrtot = (sfr.disc_sfr100 + sfr.bulge_sfr100 + sfr.burst_sfr100) * 100. ; in Msun/yr
jeje_histo, sfrtot, 0, 1000, 1000, h
plot,h.x,h.dn,psym=10,xtitle='SFR (Msun/yr)',ytitle='N',charsize=1.5



; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; high-z constraints
if keyword_set(high_z) then begin 
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   



; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
endif ; (keyword high_z)
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


device,/close
set_plot,'x'


end

