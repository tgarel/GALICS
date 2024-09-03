pro toulouse,psfile=psfile

if keyword_set(psfile) then begin 
   set_plot,'ps'
   device,filename=psfile,/color,bits_per_pixel=8
endif else begin 
   window,retain=2,xsize=720,ysize=720
   device,decomposed=0
endelse
loadct,39


toulouse_set_runs,'BestFit',runw1,runw3
toulouse_make_plots,runw1,runw3,3.0,titre='z = 3'
toulouse_make_plots,runw1,runw3,4.0,titre='z = 4'
toulouse_make_plots,runw1,runw3,5.0,titre='z = 5'
toulouse_make_plots,runw1,runw3,6.0,titre='z = 6'

if keyword_set(psfile) then begin 
   device,/close
   set_plot,'x'
endif


end

pro toulouse_make_plots,runw1,runw3,z,titre=titre

  toulouse_read_run,runw1,z,w1_gbk
  toulouse_read_run,runw3,z,w3_gbk
  w1dt = w1_gbk(0).deltat * 1.d9
  w3dt = w3_gbk(0).deltat * 1.d9

  ; plot accretion rate vs. Mh
;;   plot_oo,w1_gbk.mfof,w1_gbk.hmacc/w1dt*0.18,psym=3,yr=[1.,1.d4],ytitle='dMacc/dt [Msun / yr]',$
;;           xtitle='Mh [Msun]',title=titre
;;   oplot,w3_gbk.mfof,w3_gbk.hmacc/w3dt*0.18,psym=3,col=250

  ; plot distribution of accretion rates
  jeje_histo,alog10(w1_gbk.hmacc/w1dt*0.18),0.,3.5,30,hw1
  jeje_histo,alog10(w3_gbk.hmacc/w3dt*0.18),0.,3.5,30,hw3
  plot_io,hw1.x,hw1.dn/hw1.dx/runw1.volumempc,psym=10,thick=3,xtitle='Log dMacc/dt (Msun/yr)',$
          ytitle='dN/dLog(dMacc/dt) / Mpc!u3!n',charsize=1.5,xr=[1,3.5],yr=[1.d-6,1.d-1],/ys,$
          title=titre
  oplot,hw3.x,hw3.dn/hw3.dx/runw3.volumempc,psym=10,thick=3,col=250
  oplot,[0,4],[1.,1.]/runw3.volumempc*10.,linestyle=2,col=250
  oplot,[0,4],[1.,1.]/runw1.volumempc*10.,linestyle=2

  ii = where(hw3.x ge 2.0,ni) 
  if ni ne 0 then ni3 = total(hw3.dn(ii))
  ii = where(hw1.x ge 2.0,ni) 
  if ni ne 0 then ni1 = total(hw1.dn(ii))
  if ni1 ne 0 then print,'n1/n3 = ',float(ni1)/float(ni3),$ 
                         string(float(ni3)/runw3.volumempc,format='(e14.6)')

;  jeje_histo,alog10(w1_gbk.maccstream/w1dt),0.,3.5,30,hw1
;  jeje_histo,alog10(w3_gbk.maccstream/w3dt),0.,3.5,30,hw3
;  oplot,hw3.x,hw3.dn/hw3.dx/runw3.volumempc,psym=10,thick=3,col=250,linestyle=2
;  oplot,hw1.x,hw1.dn/hw1.dx/runw1.volumempc,psym=10,thick=3,linestyle=2

  ; plot distribution of masses
  jeje_histo,alog10(w1_gbk.mfof),10.,14,30,hw1
  jeje_histo,alog10(w3_gbk.mfof),10.,14,30,hw3
  plot_io,hw1.x,hw1.dn/hw1.dx/runw1.volumempc,psym=10,thick=3,xtitle='Log Mh (Msun)',$
          ytitle='dN/dLog(Mh) / Mpc!u3!n',charsize=1.5,xr=[10,14],/xs,yr=[1.d-6,1.d-1],/ys,title=titre
  oplot,hw3.x,hw3.dn/hw3.dx/runw3.volumempc,psym=10,thick=3,col=250
;  oplot,hw3.x+0.7,hw3.dn/hw3.dx/runw3.volumempc,psym=10,thick=2,col=250,linestyle=1
  oplot,[10,14],[1.,1.]/runw3.volumempc*10.,linestyle=2,col=250
  oplot,[10,14],[1.,1.]/runw1.volumempc*10.,linestyle=2
  
  ii = where(hw3.x ge 12.0,ni) 
  if ni ne 0 then ni3 = total(hw3.dn(ii))/runw3.volumempc else ni3 = 0.0
  ii = where(hw1.x ge 12.0,ni) 
  if ni ne 0 then ni1 = total(hw1.dn(ii))/runw1.volumempc else ni1 = 0.0
  print,'n1 =',string(float(ni1),format='(e14.6)')
  print,'n3 =',string(float(ni3),format='(e14.6)')


  ; plot star formation rate vs. gas accretion rate (from streams)
;;   plot_oo,w1_gbk.maccstream/w1dt,w1_gbk.mstarform/w1dt,psym=3,xtitle='Acc. rate (Msun/yr)',$
;;           ytitle='SFR (Msun/yr)',xr=[1.,1.d3],/xs,yr=[1.,1.d3],/ys
;;   oplot,w3_gbk.maccstream/w3dt,w3_gbk.mstarform/w3dt,psym=2,col=250,symsize=0.2
;;   oplot,[1.d-3,1.d3],[1.d-3,1.d3],col=200,thick=3

  ; plot star formation rate vs. gas accretion rate (from streams + wind)
;;   plot_oo,(w1_gbk.maccstream+w1_gbk.maccfountain)/w1dt,w1_gbk.mstarform/w1dt,psym=3,$
;;           xtitle='Acc. rate (Msun/yr)',ytitle='SFR (Msun/yr)',xr=[1.,1.d3],/xs,yr=[1.,1.d3],/ys
;;   oplot,(w3_gbk.maccstream+w3_gbk.maccfountain)/w3dt,w3_gbk.mstarform/w3dt,psym=2,col=250,symsize=0.2
;;   oplot,[1.d-3,1.d3],[1.d-3,1.d3],col=200,thick=3


end


pro toulouse_read_run,run,z,gbk
  
  timestep = z2ts(z,run.snapzdir)
  for i = run.dirmin, run.dirmax do begin 
     outputdir = run.outputdir + strtrim(i,2) + '/'
     if i eq run.dirmin then begin
        read_gbk,outputdir,timestep, gbk
     endif else begin 
        read_gbk,outputdir,timestep, gbkt
        gbk = [gbk,gbkt]
     endelse
  endfor
  
end ; pro toulouse_read_run


pro toulouse_set_runs,name,runw1,runw3
  
  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  
  runw1 = run
  runw1.outputdir = '/Users/blaizot/Documents/ASTRO/data/GM-RUNS/512-100h-1Mpc/'+name+'/'
  runw1.snapzdir=runw1.outputdir +'/1/'
  runw1.dirmin=1
  runw1.dirmax=15
  runw1.name=name+'W1'
  runw1.volumempc = 142.857^3. 
  runw1.hubble=0.7

  runw3 = run
  runw3.outputdir = '/Users/blaizot/Documents/ASTRO/data/GM-RUNS/512-100h-1Mpc-W3/'+name+'/'
  runw3.snapzdir=runw3.outputdir +'/1/'
  runw3.dirmin=1
  runw3.dirmax=15
  runw3.name=name+'W3'
  runw3.volumempc = 136.9863^3. 
  runw3.hubble=0.73

end ; pro toulouse_set_runs

