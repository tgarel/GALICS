pro plot_high_z_lfs;,run, nosubdir=nosubdir

; z de gab04: z1,  z2,  z2_3,  z4,  z5_6                                      

;;   window,retain=2
;;   device,decomposed=0
  loadct,39

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  ;run.outputdir='/data/garel/GM-RUNS/512-100h-1Mpc-W3/kenn_sfr/sf20_feedB002_best3/'
  ;run.snapzdir='/data/garel/GM-RUNS/512-100h-1Mpc-W3/kenn_sfr/sf20_feedB002_best3/1/'
  run.outputdir='/data/garel/GM-RUNS/512-100h-1Mpc-W3/kenn_sfr/analyticFB_io_ej/'
  run.snapzdir='/data/garel/GM-RUNS/512-100h-1Mpc-W3/kenn_sfr/analyticFB_io_ej/1/'
  ;run.outputdir='/data/garel/GM-RUNS/512-100h-1Mpc/kenn_sfr/sf1_feedB020/'
  ;run.snapzdir='/data/garel/GM-RUNS/512-100h-1Mpc/kenn_sfr/sf1_feedB020/1/'
  ;run.outputdir='/data/garel/GM-RUNS/1024-100h-1Mpc-W3/kenn_sfr/sf5_feedB002/'
  ;run.snapzdir='/data/garel/GM-RUNS/1024-100h-1Mpc-W3/kenn_sfr/sf5_feedB002/1/'

  run.dirmin=1
  run.dirmax= 15 ;         111
  run.name='a'

  run.volumempc = 136.9863^3.     ; W3
  run.hubble=0.730               ; W3
  ;run.volumempc = 142.85714^3.    ; W1
  ;run.hubble=0.700                ; W1


;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; z = 3
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  set_plot, 'ps'
  device, filename='plotps_lfs/512-100h-1Mpc-W3/kenn_sfr/sf20_feedB002_analytFB_io_ej_noext.ps',/color,bits_per_pixel=8

TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]

z        = 3. 
timestep = z2ts(z,run.snapzdir)
if not keyword_set(nosubdir) then begin 
   for i = run.dirmin,run.dirmax do begin 
      outputdir= run.outputdir + strtrim(i,2) + '/'
      if i eq run.dirmin then begin
         read_gal_mags,outputdir,timestep,filters,mags,/total
         read_gal_mags,outputdir,timestep,filters,magsne,/total,/noext
         read_gal_results,outputdir,timestep,galres,grinfo
      endif else begin 
         read_gal_mags,outputdir,timestep,filters,magst,/total
         mags = [mags,magst]
         read_gal_mags,outputdir,timestep,filters,magstne,/total,/noext
         magsne = [magsne,magstne]
         read_gal_results,outputdir,timestep,galrest,grinfo
         galres = [galres,galrest]
      endelse
   endfor
endif else begin 
   read_gal_mags,run.outputdir,timestep,filters,mags,/total
   read_gal_mags,run.outputdir,timestep,filters,magsne,/total,/noext
   read_gal_results,run.outputdir,timestep,galres,grinfo   
endelse

mini = -26.
maxi = -15.
nbins = 100.
bin = (maxi-mini)/nbins

min_mass = 8.
max_mass = 12.
;-------------------drory07 m_star data ---------------------------------
;; mstar  = galres.disc_minstar + galres.bulge_minstar + galres.burst_minstar
;; lmstar = alog10(mstar+1.e-5) + 11.0 ; in Msun 
;; jeje_histo, lmstar, min_mass, max_mass, nbins, hmstar
;; plot_io,hmstar.x,hmstar.dn/hmstar.dx/run.volumempc,psym=10,xtitle='Log!d10!n M!dstar!n [M!d!9n!n!3]', $  
;;         ytitle='dN / dLog!d10!nM!dstar!n / Mpc!u3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,0.1],/ys, title="Drory3-4"
;; plot_obs_data,'Drory_Mstar_z2',run.hubble, symb=4


;-------------------gab04 1500A filter-------------------------
jeje_histo,mags.UV_1500A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_1500A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1500!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.6,xr=[-18,-24],/xs,yr=[1.e-5,1.e-2],/ys,thick=5,charthick=6,title='UV LF';,title='Gabz3 + reddy_z3'
plot_obs_data,'gab04_1500_z3',run.hubble 
plot_obs_data,'reddy08_1700A_z3',run.hubble,symb=4 
plot_obs_data,'arnouts05_1500A_z3',run.hubble,symb=6
plot_obs_data,'sawicki06_1700A_z3',run.hubble,symb=2
;oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=5
lim = 1./136.9863^3/bin
oplot,[-15.,-26.],[lim,lim],linestyle=1,thick=2

XYOUTS,-22.5,6d-4," z !9A!3 3",charsize=1.9,charthick=6.

legend,['Reddy et Steidel (2008)','Arnouts et al. (2005)','Sawicki et Thompson (2006)','Gabasch et al. (2004)'],textcolors=[135,254,3,2],charsize=1.,charthick=5,psym=[4,6,2,7],colors=[135,254,3,2],/right;,/left,/bottom

;; oplot,[-21.3,-21.3],[0.007,0.007],psym=4,color=135,thick=5
;; XYOUTS,-21.5,0.0065,'Reddy et Steidel (2008)',color=135,charthick=4.
;; oplot,[-21.3,-21.3],[0.005,0.005],psym=6,color=254,thick=5
;; XYOUTS,-21.5,0.0046,'Arnouts et al. (2005)',color=254,charthick=4.
;; oplot,[-21.3,-21.3],[0.0033,0.0033],psym=2,color=3,thick=5
;; XYOUTS,-21.5,0.0031,'Sawicki et Thompson (2006)',color=3,charthick=4.
;; oplot,[-21.3,-21.3],[0.0023,0.0023],psym=7,thick=5,color=2
;; XYOUTS,-21.5,0.00215,'Gabasch et al. (2004)',charthick=4.,color=2
;; RECTANGLE,-21.1,0.0015,-2.9,0.0095

;; oplot,[-16.4,-16.9],[0.0001,0.0001],thick=4
;; XYOUTS,-17.1,0.00009,'Model with extinction',charthick=4.,charsize=1
;; oplot,[-16.4,-16.9],[0.00004,0.00004],thick=4,linestyle=2
;; XYOUTS,-17.1,0.000038,'Model without extinction',charthick=4.,charsize=1

;--------------------gab04 2800A filter-------------------------
jeje_histo,mags.UV_2800A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_2800A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d2800!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='Gabz3'
plot_obs_data,'gab04_2800_z3',run.hubble  
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;; ;-------------------gab04 SDSS_u filter-------------------------
jeje_histo,mags.SDSS_u,mini,maxi,nbins,h
jeje_histo,magsne.SDSS_u,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!du_sdss!n ',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='Gabz3'
plot_obs_data,'gab04_u_z3',run.hubble 
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;; ;-------------------gab04 JOHNSON_B filter-------------------------
jeje_histo,mags.JOHNSON_BAB,mini,maxi,nbins,h
jeje_histo,magsne.JOHNSON_BAB,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!dB!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='Gabz3'
plot_obs_data,'gab04_B_z3',run.hubble    
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;; ;-------------------gab04 SDSS_g filter-------------------------
jeje_histo,mags.SDSS_g,mini,maxi,nbins,h
jeje_histo,magsne.SDSS_g,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!dg_sdss!n ',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='Gabz3'
plot_obs_data,'gab04_g_z3',run.hubble 
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3


;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; z = 4
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
z        = 4
timestep = z2ts(z,run.snapzdir)
if not keyword_set(nosubdir) then begin 
   for i = run.dirmin,run.dirmax do begin 
      outputdir= run.outputdir + strtrim(i,2) + '/'
      if i eq run.dirmin then begin
         read_gal_mags,outputdir,timestep,filters,mags,/total
         read_gal_mags,outputdir,timestep,filters,magsne,/total,/noext
         read_gal_results,outputdir,timestep,galres,grinfo
      endif else begin 
         read_gal_mags,outputdir,timestep,filters,magst,/total
        mags = [mags,magst]
         read_gal_mags,outputdir,timestep,filters,magstne,/total,/noext
         magsne = [magsne,magstne]
         read_gal_results,outputdir,timestep,galrest,grinfo
         galres = [galres,galrest]
      endelse
   endfor
endif else begin 
   read_gal_mags,run.outputdir,timestep,filters,mags,/total
   read_gal_mags,run.outputdir,timestep,filters,magsne,/total,/noext
   read_gal_results,run.outputdir,timestep,galres,grinfo
endelse
mini = -26
maxi = -15
nbins = 100

min_mass = 8.
max_mass = 12.
;-------------------drory07 m_star data ---------------------------------
mstar  = galres.disc_minstar + galres.bulge_minstar + galres.burst_minstar
lmstar = alog10(mstar+1.e-5) + 11.0 ; in Msun 
jeje_histo, lmstar, min_mass, max_mass, nbins, hmstar
plot_io,hmstar.x,hmstar.dn/hmstar.dx/run.volumempc,psym=10,xtitle='Log!d10!n M!dstar!n [M!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nM!dstar!n / Mpc!u3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,0.1],/ys, title="Drory3-4"
plot_obs_data,'Drory_Mstar_z3',run.hubble, symb=4

;-------------------Bouwens07 i_775 filter pour z~4 (1600A) -----------------
jeje_histo,mags.iwata07_1570A,mini,maxi,nbins,h
jeje_histo,magsne.iwata07_1570A,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1600!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='Bouz4'
plot_obs_data,'bouwens07_1600A_z4',run.hubble,sym=5  
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3


;-------------Gab04 2800 filter------------------------------------------------

jeje_histo,mags.UV_2800A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_2800A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d2800!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_2800_z4'
plot_obs_data,'gab04_2800_z4',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------Gab04 1500 filter------------------------------------------------

jeje_histo,mags.UV_2800A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_2800A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1500!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_1500_z4'
plot_obs_data,'gab04_1500_z4',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------Gab04 u filter------------------------------------------------

jeje_histo,mags.SDSS_u,mini,maxi,nbins,h
jeje_histo,magsne.SDSS_u,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='U_sdss',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_u_z4'
plot_obs_data,'gab04_u_z4',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------Gab04 Johnson B filter------------------------------------------------

jeje_histo,mags.JOHNSON_BAB,mini,maxi,nbins,h
jeje_histo,magsne.JOHNSON_BAB,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='Johnson B',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_B_z4'
plot_obs_data,'gab04_B_z4',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------Gab04 g sdss filter----------------------------------------------
jeje_histo,mags.SDSS_g,mini,maxi,nbins,h
jeje_histo,magsne.SDSS_g,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='g_sdss',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_g_z4'
plot_obs_data,'gab04_g_z4',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; z = 5
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
z        = 5
timestep = z2ts(z,run.snapzdir)
if not keyword_set(nosubdir) then begin 
   for i = run.dirmin,run.dirmax do begin 
      outputdir= run.outputdir + strtrim(i,2) + '/'
      if i eq run.dirmin then begin
         read_gal_mags,outputdir,timestep,filters,mags,/total
         read_gal_mags,outputdir,timestep,filters,magsne,/total,/noext
         read_gal_results,outputdir,timestep,galres,grinfo
      endif else begin 
         read_gal_mags,outputdir,timestep,filters,magst,/total
         mags = [mags,magst]
         read_gal_mags,outputdir,timestep,filters,magstne,/total,/noext
         magsne = [magsne,magstne]
         read_gal_results,outputdir,timestep,galrest,grinfo
         galres = [galres,galrest]
      endelse
   endfor
endif else begin 
   read_gal_mags,run.outputdir,timestep,filters,mags,/total
   read_gal_mags,run.outputdir,timestep,filters,magsne,/total,/noext
   read_gal_results,run.outputdir,timestep,galres,grinfo
endelse
mini = -26
maxi = -15
nbins = 100

min_mass = 8.
max_mass = 12.
;-------------------drory07 m_star data ---------------------------------
mstar  = galres.disc_minstar + galres.bulge_minstar + galres.burst_minstar
lmstar = alog10(mstar+1.e-5) + 11.0 ; in Msun 
jeje_histo, lmstar, min_mass, max_mass, nbins, hmstar
plot_io,hmstar.x,hmstar.dn/hmstar.dx/run.volumempc,psym=10,xtitle='Log!d10!n M!dstar!n [M!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nM!dstar!n / Mpc!u3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,0.1],/ys, title="Drory4-5"
plot_obs_data,'Drory_Mstar_z4',run.hubble, symb=4

;-------------------bouwens07 i_850 filter pour z~5 et 6 (1600A) --------------
;; jeje_histo,mags.UV_1500A_Gab04,mini,maxi,nbins,h
;; jeje_histo,magsne.UV_1500A_Gab04,mini,maxi,nbins,hne
;; plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1600!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-18,-24],/xs,yr=[1.e-5,1.e-2],/ys,thick=3,title='Bouz5'
;; plot_obs_data,'bouwens07_1600A_z5',run.hubble  ; Change z
;; oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;---------iwata07 1570A filter pour z~5 ;  filtre pris arbitrairement --------
jeje_histo,mags.iwata07_1570A,mini,maxi,nbins,h
jeje_histo,magsne.iwata07_1570A,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1570!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='Iwata_z5'
plot_obs_data,'iwata_z5',run.hubble,symb=2  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3
plot_obs_data,'mclure_z5',run.hubble,symb=4  ; Change z
plot_obs_data,'bouwens07_1600A_z5',run.hubble,symb=5  ; Change z

;laisser decommenté ci dessous -> z=5.5

;---Mclure 1500A filter pour z~5 ; filtre pris arbitrairement (= Gabasch 1500 one --------
;; jeje_histo,mags.UV_1500A_Gab04,mini,maxi,nbins,h
;; jeje_histo,magsne.UV_1500A_Gab04,mini,maxi,nbins,hne
;; plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1500!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-18,-24],/xs,yr=[1.e-5,1.e-2],/ys,thick=3,title='Mclure_z5'
;; plot_obs_data,'mclure_z5',run.hubble  ; Change z
;; oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------------reddy08 R_G_steidel filter pour z~2 (1700A) --------------
;; jeje_histo,mags.,mini,maxi,nbins,h
;; plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1700!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-18,-24],/xs,yr=[1.e-5,1.e-2],/ys,thick=3
;; plot_obs_data,'eddy08_1700A_z2',run.hubble  

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; z = 5.5
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
z        = 5.25
timestep = z2ts(z,run.snapzdir)
if not keyword_set(nosubdir) then begin 
   for i = run.dirmin,run.dirmax do begin 
      outputdir= run.outputdir + strtrim(i,2) + '/'
      if i eq run.dirmin then begin
         read_gal_mags,outputdir,timestep,filters,mags,/total
         read_gal_mags,outputdir,timestep,filters,magsne,/total,/noext
         read_gal_results,outputdir,timestep,galres,grinfo
      endif else begin 
         read_gal_mags,outputdir,timestep,filters,magst,/total
         mags = [mags,magst]
         read_gal_mags,outputdir,timestep,filters,magstne,/total,/noext
         magsne = [magsne,magstne]
         read_gal_results,outputdir,timestep,galrest,grinfo
         galres = [galres,galrest]
      endelse
   endfor
endif else begin 
   read_gal_mags,run.outputdir,timestep,filters,mags,/total
   read_gal_mags,run.outputdir,timestep,filters,magsne,/total,/noext
   read_gal_results,run.outputdir,timestep,galres,grinfo
endelse
mini = -26
maxi = -15
nbins = 100

;-------------Gab04 2800 filter------------------------------------------------

jeje_histo,mags.UV_2800A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_2800A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d2800!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-18,-24],/xs,yr=[1.e-6,1.e-1],/ys,thick=3,title='gab04_2800_z5_6'
plot_obs_data,'gab04_2800_z5_6',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------Gab04 1500 filter------------------------------------------------

jeje_histo,mags.UV_1500A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_1500A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1500!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.6,xr=[-18,-24],/xs,yr=[1.e-5,1.e-2],/ys,thick=6,title='UV LF',charthick=6
plot_obs_data,'gab04_1500_z5_6',run.hubble  ; Change z
;oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=5

;--------------- z = 5 or 6 observations mais bon... ---------------------
plot_obs_data,'bouwens07_1600A_z6',run.hubble,symb=4
plot_obs_data,'mclure_z6',run.hubble ,symb=6
;plot_obs_data,'mclure_z5',run.hubble,symb=5  ; Change z
plot_obs_data,'bouwens07_1600A_z5',run.hubble,symb=5  ; Change z
;plot_obs_data,'iwata_z5',run.hubble,symb=2  ; Change z

XYOUTS,-22.,7d-4," z !9A!3 5.5",charsize=1.5,charthick=6.

legend,['Gabasch et al. (2004) 5<z<6','Bouwens et al. (2007) z !9A!3 5','Bouwens et al. (2007) z !9A!3 6','Mc Lure (2009) z !9A!3 6'],textcolors=[2,254,254,3],charsize=1.,charthick=5,psym=[7,5,4,6],colors=[2,254,254,3],/right;,/left,/bottom

;; oplot,[-21.2,-21.2],[0.007,0.007],psym=7,color=2,thick=5
;; XYOUTS,-21.4,0.0065,'Gabasch et al. (2004) 5<z<6',color=2,charthick=4.
;; oplot,[-21.2,-21.2],[0.005,0.005],psym=5,color=254,thick=5
;; XYOUTS,-21.4,0.0046,'Bouwens et al. (2007) z !9A!3 5',color=254,charthick=4.
;; oplot,[-21.2,-21.2],[0.0033,0.0033],psym=4,color=254,thick=5
;; XYOUTS,-21.4,0.0031,'Bouwens et al. (2007) z !9A!3 6',color=254,charthick=4.
;; oplot,[-21.2,-21.2],[0.0023,0.0023],psym=6,color=3,thick=5
;; XYOUTS,-21.4,0.00215,'Mc Lure (2009) z !9A!3 6',charthick=4.,color=3
;; RECTANGLE,-21.,0.0015,-3.2,0.0095;,color=3

;; oplot,[-22.1,-22.7],[0.0007,0.0007],thick=4
;; XYOUTS,-22.9,0.00065,'Model with extinction',charthick=4.,charsize=0.9
;; oplot,[-22.1,-22.7],[0.0004,0.0004],linestyle=2,thick=4
;; XYOUTS,-22.9,0.00037,'Model without extinction',charthick=4.,charsize=0.9

;-------------Gab04 u filter------------------------------------------------

jeje_histo,mags.SDSS_u,mini,maxi,nbins,h
jeje_histo,magsne.SDSS_u,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='U_sdss',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_u_z5_6'
plot_obs_data,'gab04_u_z5_6',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------Gab04 Johnson B filter------------------------------------------------

jeje_histo,mags.JOHNSON_BAB,mini,maxi,nbins,h
jeje_histo,magsne.JOHNSON_BAB,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='Johnson B',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_B_z5_6'
plot_obs_data,'gab04_B_z5_6',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;-------------Gab04 g sdss filter------------------------------------------------

jeje_histo,mags.SDSS_g,mini,maxi,nbins,h
jeje_histo,magsne.SDSS_g,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='g_sdss',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='gab04_g_z5_6'
plot_obs_data,'gab04_g_z5_6',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3


;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; z = 6
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
z        = 6
timestep = z2ts(z,run.snapzdir)
if not keyword_set(nosubdir) then begin 
   for i = run.dirmin,run.dirmax do begin 
      outputdir= run.outputdir + strtrim(i,2) + '/'
      if i eq run.dirmin then begin
         read_gal_mags,outputdir,timestep,filters,mags,/total
         read_gal_mags,outputdir,timestep,filters,magsne,/total,/noext
         read_gal_results,outputdir,timestep,galres,grinfo
      endif else begin 
         read_gal_mags,outputdir,timestep,filters,magst,/total
         mags = [mags,magst]
         read_gal_mags,outputdir,timestep,filters,magstne,/total,/noext
         magsne = [magsne,magstne]
         read_gal_results,outputdir,timestep,galrest,grinfo
         galres = [galres,galrest]
      endelse
   endfor
endif else begin 
   read_gal_mags,run.outputdir,timestep,filters,mags,/total
   read_gal_mags,run.outputdir,timestep,filters,magsne,/total,/noext
   read_gal_results,run.outputdir,timestep,galres,grinfo
endelse
mini = -26
maxi = -15
nbins = 100

;---Mclure 1500A filter pour z~6 ; filtre pris arbitrairement (= Gabasch 1500 one --------
jeje_histo,mags.UV_1500A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_1500A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1500!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='Mclure_z6'
plot_obs_data,'mclure_z6',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3

;----------------- Bouwens07 1600 A (changer filtre: mettre UV_1600Ang )-------------------------------------
jeje_histo,mags.UV_1500A_Gab04,mini,maxi,nbins,h
jeje_histo,magsne.UV_1500A_Gab04,mini,maxi,nbins,hne
plot_io,h.x,h.dn/h.dx/run.volumempc,psym=10,xtitle='M!d1600!n !6!sA!r!u!9 %!6!n',ytitle='!4U!3 [Mpc!u-3!n mag!u-1!n]',charsize=1.5,xr=[-16,-26],/xs,yr=[1.e-6,1.e-2],/ys,thick=3,title='bouwens07_1600A_z6'
plot_obs_data,'bouwens07_1600A_z6',run.hubble  ; Change z
oplot,hne.x,hne.dn/hne.dx/run.volumempc,psym=10,linestyle=2,thick=3


;++++++++++++++++++++++++++++++++++++++++++++++++++++++
;           Infrared galaxies  High redshift univer
;                       z = 3.0
;++++++++++++++++++++++++++++++++++++++++++++++++++++++
; define z=3.00 snapshot
;; z        = 3.0
;; timestep = z2ts(z,run.snapzdir)
;; ; read useful information 
;; if not keyword_set(nosubdir) then begin 
;;    for i = run.dirmin,run.dirmax do begin 
;;       outputdir= run.outputdir + strtrim(i,2) + '/'
;;       if i eq run.dirmin then begin
;;          read_bol_lum,outputdir,timestep,bol
;;          read_gal_mags,outputdir,timestep,filters,mags,/total
;;       endif else begin 
;;          read_bol_lum,outputdir,timestep,bolt
;;          bol = [bol,bolt]
;;          read_gal_mags,outputdir,timestep,filters,magst,/total
;;          mags = [mags,magst]
;;       endelse
;;    endfor
;; endif else begin 
;;    read_bol_lum,run.outputdir,timestep,bol
;;    read_gal_mags,run.outputdir,timestep,filters,mags,/total   
;; endelse


; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; IV.1   IR bolometric luminosities 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;; irbol = bol.ir*(run.hubble^2)
;; lirbol = alog10(irbol)
;; jeje_histo,  lirbol, 7.0, 14.0, 100, hirbol
;; plot_io,hirbol.x,hirbol.dn/hirbol.dx/(run.volumempc*run.hubble^3.0),psym=10,xtitle='Log!d10!n L!dir!n [h!u-2!n L!d!9n!n!3]', $  
;;         ytitle='dN / dLog!d10!nL!dir!n / Mpc!u3!nh!u-3!n',charsize=1.5,xr=[8,12.5],/xs,yr=[1.e-6,1.0],/ys

 device,/close

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; IV.2 Set of infrared filters
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ir_lfs,mags,run.volumempc*(run.hubble^3),run.hubble,z, no_obs_data=no_obs_data





end

