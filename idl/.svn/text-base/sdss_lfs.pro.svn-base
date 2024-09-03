pro sdss_lfs,mags,volumempc,hubble,no_obs_data=no_obs_data

; INPUTS : 
; - mags      : a structure array containing at least the SDSS filters 
;   (mag(*).sdss_x) 
; - volumempc : the volume occupyied by the galaxies in mag, in Mpc^3
; - hubble    : the hubble constant (in units of 100km/s/Mpc)
; KEYWORDS : 
; - no_obs_data: if set, dont overplot observations

mini = -26
maxi = -10
nbins = 100

jeje_histo,mags.sdss_u,mini,maxi,nbins,hu
plot_io,hu.x,hu.dn/hu.dx/volumempc,psym=10,xtitle='u',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
;plot_io,hu.x-5.*alog10(hubble),hu.dn/hu.dx/volumempc/hubble^3,psym=10,xtitle='u - 5 Log!d10!n h',ytitle='!4U!3 [h!u3!nMpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
if not keyword_set(no_obs_data) then plot_obs_data,'SDSS_u_LF',hubble

jeje_histo,mags.sdss_g,mini,maxi,nbins,hg
plot_io,hu.x,hg.dn/hg.dx/volumempc,psym=10,xtitle='g',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
;plot_io,hu.x-5.*alog10(hubble),hg.dn/hg.dx/volumempc/hubble^3,psym=10,xtitle='g - 5 Log!d10!n h',ytitle='!4U!3 [h!u3!nMpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
if not keyword_set(no_obs_data) then plot_obs_data,'SDSS_g_LF',hubble,datadir=datadir

jeje_histo,mags.sdss_r,mini,maxi,nbins,hr
plot_io,hu.x,hr.dn/hr.dx/volumempc,psym=10,xtitle='r',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
;plot_io,hu.x-5.*alog10(hubble),hr.dn/hr.dx/volumempc/hubble^3,psym=10,xtitle='r - 5 Log!d10!n h',ytitle='!4U!3 [h!u3!nMpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
if not keyword_set(no_obs_data) then plot_obs_data,'SDSS_r_LF',hubble,datadir=datadir

jeje_histo,mags.sdss_i,mini,maxi,nbins,hi
plot_io,hu.x,hi.dn/hi.dx/volumempc,psym=10,xtitle='i',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
;plot_io,hu.x-5.*alog10(hubble),hi.dn/hi.dx/volumempc/hubble^3,psym=10,xtitle='i - 5 Log!d10!n h',ytitle='!4U!3 [h!u3!nMpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
if not keyword_set(no_obs_data) then plot_obs_data,'SDSS_i_LF',hubble,datadir=datadir

jeje_histo,mags.sdss_z,mini,maxi,nbins,hz
plot_io,hu.x,hz.dn/hz.dx/volumempc,psym=10,xtitle='z',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
;plot_io,hu.x-5.*alog10(hubble),hz.dn/hz.dx/volumempc/hubble^3,psym=10,xtitle='z - 5 Log!d10!n h',ytitle='!4U!3 [h!u3!nMpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-12,-25],/xs,yr=[1.e-6,0.2],/ys
if not keyword_set(no_obs_data) then plot_obs_data,'SDSS_z_LF',hubble,datadir=datadir

end

