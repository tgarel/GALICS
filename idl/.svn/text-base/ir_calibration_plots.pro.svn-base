pro ir_calibration_plots, run, psfile,no_obs_data=no_obs_data

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

;++++++++++++++++++++++++++++++++++++++++++++++++++++++
;           Infrared galaxies Local Universe
;                       z = 0 
;++++++++++++++++++++++++++++++++++++++++++++++++++++++
; define z=0 snapshot
z        = 0.0
timestep = z2ts(z,run.snapzdir)
; read useful information 
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + strtrim(i,2) + '/'
   if i eq run.dirmin then begin
      read_bol_lum,outputdir,timestep,bol
      read_gal_mags,outputdir,timestep,filters,mags,/total
 endif else begin 
      read_bol_lum,outputdir,timestep,bolt
      bol = [bol,bolt]
      read_gal_mags,outputdir,timestep,filters,magst,/total
      mags = [mags,magst]
   endelse
endfor


; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; I.1   IR bolometric luminosities 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

irbol = bol.ir*(run.hubble^2)
lirbol = alog10(irbol)
jeje_histo,  lirbol, 7.0, 14.0, 100, hirbol
plot_io,hirbol.x,hirbol.dn/hirbol.dx/(run.volumempc*run.hubble^3.0),psym=10,xtitle='Log!d10!n L!dir!n [h!u-2!n L!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nL!dir!n / Mpc!u3!nh!u-3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys

; overplot the analytical fit to the Takeuchi 2003 data, correctic the 60mic lum
; to bol lum 
   lum_star = 8.85e8
   phi_star = 2.34e-2
   alpha = 1.23
   sigma = 0.724   
   lum_x = hirbol.x + alog10(2.5)
   lum_func = phi_star*(10.0^lum_x/lum_star)^(1.0-alpha)*exp(-(1.0/(2.0*sigma^2))*alog10(1.0+(10.0^lum_x/lum_star))^2.0)
   oplot,lum_x,lum_func,psym=7,symsize=0.3    



;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; I.2 Set of infrared filters
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ir_lfs,mags,run.volumempc*(run.hubble^3),run.hubble,z, no_obs_data=no_obs_data





;++++++++++++++++++++++++++++++++++++++++++++++++++++++
;           Infrared galaxies High Redshift universe
;                       z = 1.0 
;++++++++++++++++++++++++++++++++++++++++++++++++++++++
; define z=1.00 snapshot
z        = 1.0
timestep = z2ts(z,run.snapzdir)
; read useful information 
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + strtrim(i,2) + '/'
   if i eq run.dirmin then begin
      read_bol_lum,outputdir,timestep,bol
      read_gal_mags,outputdir,timestep,filters,mags,/total
 endif else begin 
      read_bol_lum,outputdir,timestep,bolt
      bol = [bol,bolt]
      read_gal_mags,outputdir,timestep,filters,magst,/total
      mags = [mags,magst]
   endelse
endfor


; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; II.1   IR bolometric luminosities 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

irbol = bol.ir*(run.hubble^2)
lirbol = alog10(irbol)
jeje_histo,  lirbol, 7.0, 14.0, 100, hirbol
plot_io,hirbol.x,hirbol.dn/hirbol.dx/(run.volumempc*run.hubble^3.0),psym=10,xtitle='Log!d10!n L!dir!n [h!u-2!n L!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nL!dir!n / Mpc!u3!nh!u-3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys
plot_obs_data,'Caputi_IRbol_z1',run.hubble
plot_obs_data,'LeFloch_IRbol_z1',run.hubble

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; II.2 Set of infrared filters
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ir_lfs,mags,run.volumempc*(run.hubble^3),run.hubble,z, no_obs_data=no_obs_data



;++++++++++++++++++++++++++++++++++++++++++++++++++++++
;           Infrared galaxies  High redshift univer
;                       z = 2.0
;++++++++++++++++++++++++++++++++++++++++++++++++++++++
; define z=2.00 snapshot
z        = 2.0
timestep = z2ts(z,run.snapzdir)
; read useful information 
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + strtrim(i,2) + '/'
   if i eq run.dirmin then begin
      read_bol_lum,outputdir,timestep,bol
      read_gal_mags,outputdir,timestep,filters,mags,/total
 endif else begin 
      read_bol_lum,outputdir,timestep,bolt
      bol = [bol,bolt]
      read_gal_mags,outputdir,timestep,filters,magst,/total
      mags = [mags,magst]
   endelse
endfor


; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; III.1   IR bolometric luminosities 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

irbol = bol.ir*(run.hubble^2)
lirbol = alog10(irbol)
jeje_histo,  lirbol, 7.0, 14.0, 100, hirbol
plot_io,hirbol.x,hirbol.dn/hirbol.dx/(run.volumempc*run.hubble^3.0),psym=10,xtitle='Log!d10!n L!dir!n [h!u-2!n L!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nL!dir!n / Mpc!u3!nh!u-3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys
plot_obs_data,'Caputi_IRbol_z2',run.hubble

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; III.2 Set of infrared filters
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ir_lfs,mags,run.volumempc*(run.hubble^3),run.hubble,z, no_obs_data=no_obs_data





;++++++++++++++++++++++++++++++++++++++++++++++++++++++
;           Infrared galaxies  High redshift univer
;                       z = 3.0
;++++++++++++++++++++++++++++++++++++++++++++++++++++++
; define z=3.00 snapshot
z        = 3.0
timestep = z2ts(z,run.snapzdir)
; read useful information 
for i = run.dirmin,run.dirmax do begin 
   outputdir= run.outputdir + strtrim(i,2) + '/'
   if i eq run.dirmin then begin
      read_bol_lum,outputdir,timestep,bol
      read_gal_mags,outputdir,timestep,filters,mags,/total
 endif else begin 
      read_bol_lum,outputdir,timestep,bolt
      bol = [bol,bolt]
      read_gal_mags,outputdir,timestep,filters,magst,/total
      mags = [mags,magst]
   endelse
endfor


; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; IV.1   IR bolometric luminosities 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
irbol = bol.ir*(run.hubble^2)
lirbol = alog10(irbol)
jeje_histo,  lirbol, 7.0, 14.0, 100, hirbol
plot_io,hirbol.x,hirbol.dn/hirbol.dx/(run.volumempc*run.hubble^3.0),psym=10,xtitle='Log!d10!n L!dir!n [h!u-2!n L!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nL!dir!n / Mpc!u3!nh!u-3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys


;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; IV.2 Set of infrared filters
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ir_lfs,mags,run.volumempc*(run.hubble^3),run.hubble,z, no_obs_data=no_obs_data



device,/close
set_plot,'x'
end

