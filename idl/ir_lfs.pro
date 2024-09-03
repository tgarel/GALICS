pro ir_lfs,mags,volumempch,hubble,z,no_obs_data=no_obs_data

; INPUTS : 
; - mags      : a structure array containing some IR filters
; - volumempc : the volume occupyied by the galaxies in log10(Lum), in (Mpc/h)^3
; - hubble    : the hubble constant (in units of 100km/s/Mpc)
; - z         : redshift
; KEYWORDS : 
; - no_obs_data: if set, dont overplot observations

mini = 8
maxi = 12
nbins = 100

M_bol_sun  = 3.8
log_L_sun = 33.59
if z eq 0.0 then begin
; IRAS 60mic
   to_frequency = 3.e8/(60.0e-6)
   log_lum = -0.4*((mags.iras_60mic) - 51.63) + alog10(to_frequency)
   log_lum(*) = log_lum(*) - log_L_sun
   log_lum(*) = log_lum(*) + alog10(hubble^2.0)
   jeje_histo,log_lum,mini,maxi,nbins,hu
   plot_io,hu.x,hu.dn/hu.dx/volumempch,psym=10,xtitle='iras_60mic [L/h!u-2!n L!d!9n!n!3]',$
           ytitle='!4U!3 [Mpc!u-3!nh!u3!n log L !u-1!n]',$
           charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys

; overplot the analytical fit to the Takeuchi 2003 data
   lum_star = 8.85e8
   phi_star = 2.34e-2
   alpha = 1.23
   sigma = 0.724   
   lum_func = phi_star*(10.0^hu.x/lum_star)^(1.0-alpha)*exp(-(1.0/(2.0*sigma^2))*alog10(1.0+(10.0^hu.x/lum_star))^2.0)
   oplot,hu.x,lum_func,psym=7,symsize=0.3    


; IRAS 100mic
;   to_frequency = 3.e8/(100.0e-6)
;   log_lum = -0.4*((mags.iras_100mic) - 51.63) + alog10(to_frequency)
;   log_lum(*) = log_lum(*) - log_L_sun
;   log_lum(*) = log_lum(*) + alog10(hubble^2.0)
;   jeje_histo,log_lum,mini,maxi,nbins,hu
;   plot_io,hu.x,hu.dn/hu.dx/volumempch,psym=10,xtitle='iras_100mic [L/h!u-2!n L!d!9n!n!3]',$
;           ytitle='!4U!3 [Mpc!u-3!nh!u3!n log L !u-1!n]',$
;           charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys

   
; MIPS 24mic
   to_frequency = 3.e8/(24.0e-6)
   log_lum = -0.4*((mags.mips_24mic) - 51.63) + alog10(to_frequency)
   log_lum(*) = log_lum(*) - log_L_sun
   log_lum(*) = log_lum(*) + alog10(hubble^2.0)
   jeje_histo,log_lum,mini,maxi,nbins,hu
   plot_io,hu.x,hu.dn/hu.dx/volumempch,psym=10,xtitle='mips_24mic [L/h!u-2!n L!d!9n!n!3]',$
           ytitle='!4U!3 [Mpc!u-3!nh!u3!n log L !u-1!n]',$
           charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys   
   plot_obs_data,'Babbedge_MIPS_24mic_z0',hubble
endif

if z eq 1.0 then begin
; MIPS 24mic
   to_frequency = 3.e8/(24.0e-6)
   log_lum = -0.4*((mags.mips_24mic) - 51.63) + alog10(to_frequency)
   log_lum(*) = log_lum(*) - log_L_sun
   log_lum(*) = log_lum(*) + alog10(hubble^2.0)
   jeje_histo,log_lum,mini,maxi,nbins,hu
   plot_io,hu.x,hu.dn/hu.dx/volumempch,psym=10,xtitle='mips_24mic  [L/h!u-2!n L!d!9n!n!3]',$ 
           ytitle='!4U!3 [Mpc!u-3!nh!u3!n log L !u-1!n]',$
           charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys   
   plot_obs_data,'Babbedge_MIPS_24mic_z1',hubble

; ISOCAM 15mic
;   to_frequency = 3.e8/(15.0e-6)
;   log_lum = -0.4*((mags.isocam_15mic) - 51.63) + alog10(to_frequency)
;   log_lum(*) = log_lum(*) - log_L_sun
;   log_lum(*) = log_lum(*) + alog10(hubble^2.0)
;   jeje_histo,log_lum,mini,maxi,nbins,hu
;   plot_io,hu.x,hu.dn/hu.dx/volumempch,psym=10,xtitle='isocam_15mic  [L/h!u-2!n L!d!9n!n!3]',$
;           ytitle='!4U!3 [Mpc!u-3!nh!u3!n log L !u-1!n]',$
;           charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys   
endif

if z eq 2.0 then begin
; MIPS 24mic
;   to_frequency = 3.e8/(24.0e-6)
;   log_lum = -0.4*((mags.mips_24mic) - 51.63) + alog10(to_frequency)
;   log_lum(*) = log_lum(*) - log_L_sun
;   log_lum(*) = log_lum(*) + alog10(hubble^2.0)
;   jeje_histo,log_lum,mini,maxi,nbins,hu
;   plot_io,hu.x,hu.dn/hu.dx/volumempch,psym=10,xtitle='mips_24mic  [L/h!u-2!n L!d!9n!n!3]',$ 
;           ytitle='!4U!3 [Mpc!u-3!nh!u3!n log L !u-1!n]',$
;           charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys   
;   plot_obs_data,'Babbedge_MIPS_24mic_z2',hubble

endif


if z eq 3.0 then begin
; MIPS 24mic
   to_frequency = 3.e8/(24.0e-6)
   log_lum = -0.4*((mags.mips_24mic) - 51.63) + alog10(to_frequency)
   log_lum(*) = log_lum(*) - log_L_sun
   log_lum(*) = log_lum(*) + alog10(hubble^2.0)
   jeje_histo,log_lum,mini,maxi,nbins,hu
   plot_io,hu.x,hu.dn/hu.dx/volumempch,psym=10,xtitle='mips_24mic  [L/h!u-2!n L!d!9n!n!3]',$ 
           ytitle='!4U!3 [Mpc!u-3!nh!u3!n log L !u-1!n]',$
           charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys   
   plot_obs_data,'Babbedge_MIPS_24mic_z3',hubble

endif


end
