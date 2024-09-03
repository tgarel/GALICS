pro plot_haloprops,run,timestep,nosubdir=nosubdir,titre=titre

  if not keyword_set(nosubdir) then begin 
     for i = run.dirmin,run.dirmax do begin 
        if i eq run.dirmin then begin 
           read_halo_baryondec,run.outputdir+strtrim(i,2),timestep,h
           read_centralprops,run.outputdir+strtrim(i,2),timestep,g
        endif else begin 
           read_halo_baryondec,run.outputdir+strtrim(i,2),timestep,ht
           read_centralprops,run.outputdir+strtrim(i,2),timestep,gtt
           h = [h,ht]
           g = [g,gtt]
        endelse
     endfor
  endif else begin 
     read_halo_baryondec,run.outputdir,timestep,h
     read_centralprops,run.outputdir,timestep,g
  endelse
  
  jeje_histo,alog10(h.mfof)+11.,10.,13.,100,histo
  plot_io,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,xtitle='Log (Mfof)',ytitle='N (/Mpc!u3!n)/dLogM',xr=[10,13],/xs,yr=[1.d-5,0.5],/ys,thick=3
  
  jeje_histo,alog10(h.mvir)+11.,10.,13.,100,histo
  oplot,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,linestyle=2
  
  
  jeje_histo,alog10(h.mcoldgaz)+11.,10.,13,100,histo
  oplot,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,linestyle=1

  jeje_histo,alog10(h.mstar)+11.,10.,13,100,histo
  oplot,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,linestyle=3,thick=3,col=250

  jeje_histo,alog10(h.mism)+11.,10.,13,100,histo
  oplot,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,linestyle=1,thick=3
  
  jeje_histo,alog10(g.mstar)+11.,10.,13.,100,histo
  oplot,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,linestyle=2,thick=3,col=150


  ; weigh by number of galaxies 
  jeje_histo,alog10(h.mfof)+11.,10.,13.,100,histo
  plot_io,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,xtitle='Log (Mfof)',ytitle='N (/Mpc!u3!n)/dLogM',xr=[10,13],/xs,yr=[1.d-5,0.5],/ys
  jeje_histo,alog10(h.mfof)+11.0d0,10.,13.,100,histo,weight=h.nbgal
  oplot,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,thick=3

  ; plot hot/cold distribution
  coldfrac = h.macc / (h.mhot+h.macc)
  plot_oi,h.mfof*1.d11,coldfrac,xr=[1.d9,1.d13],/xs,yr=[0,1],/ys,psym=2,symsize=0.3,xtitle='Mfof',ytitle='Mstream/(Mhot+Mstream)',title=titre

  ; plot mass of winds 
  jeje_histo,alog10(h.mwind)+11.,6.,12.,100,histo
  plot_io,histo.x,histo.dn/histo.dx/run.volumempc,psym=10,linestyle=2,xr=[6,12],/xs,xtitle='Mwind',ytitle='#'


end

