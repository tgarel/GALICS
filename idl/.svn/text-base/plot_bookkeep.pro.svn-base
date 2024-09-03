pro plot_bookkeep, run, timestep, nosubdir=nosubdir

  if not keyword_set(nosubdir) then begin 
     for i = run.dirmin,run.dirmax do begin 
        outputdir= run.outputdir + strtrim(i,2) + '/'
        if i eq run.dirmin then begin
           read_gbk,outputdir,timestep, gbk
        endif else begin 
           read_gbk,outputdir,timestep, gbk
           gbk = [gbk,gbk]
        endelse
     endfor
  endif else begin 
     read_gbk,run.outputdir,timestep,gbk
  endelse
     
  dt = gbk(0).deltat * 1.d9



  ; histograms of macc stream
  mmin = 4.
  mmax  = 14.
  nbins = 100
  maccstream = alog10(gbk.maccstream)
  jeje_histo, maccstream, mmin, mmax, nbins, hstream
  plot_io,hstream.x,hstream.dn/hstream.dx/run.volumempc,psym=10,xtitle='Log!d10!n M!dstream!n [M!d!9n!n!3]', $  
        ytitle='dN / dLog!d10!nM!dstream!n / Mpc!u3!n',charsize=1.5,xr=[4,12.5],/xs,yr=[1.e-6,0.1],/ys

  ; macc stream as a fraction of halo mass
  frac_maccstream = gbk.maccstream/gbk.mfof
  mmin = 0.0
  mmax = 0.2
  nbins = 100
  jeje_histo, frac_maccstream, mmin, mmax, nbins, hstream
  plot_io,hstream.x,hstream.dn/hstream.dx/run.volumempc,psym=10,xtitle='M!dstream!n/M!dFOF!n', $  
        ytitle='dN / d(M!dstream!n/M!dFOF!n) / Mpc!u3!n',charsize=1.5,xr=[0.0,0.2],/xs,yr=[1.e-6,1.0],/ys

  ; macc stream as a fraction of mass in stars
  frac_maccstream = alog10((gbk.mstar + 1.0)/(gbk.maccstream + 1.0))
  mmin = -4.0
  mmax = 4.0
  nbins = 100
  jeje_histo, frac_maccstream, mmin, mmax, nbins, hstream
  plot_io,hstream.x,hstream.dn/hstream.dx/run.volumempc,psym=10,xtitle='Log!d10!nM!dstar!n/M!dstream!n', $  
        ytitle='dN / dLog!d10!n(M!dstar!n/M!dstream!n) / Mpc!u3!n',charsize=1.5,xr=[-4.0,4.0],/xs,yr=[1.e-6,1.0],/ys


  ; mout  as a fraction of mass in streams
  frac_mejcold = alog10((gbk.mejcold + 1.0)/(gbk.maccstream + 1.0))
  mmin = -4.0
  mmax = 4.0
  nbins = 100
  jeje_histo, frac_mejcold, mmin, mmax, nbins, hstream
  plot_io,hstream.x,hstream.dn/hstream.dx/run.volumempc,psym=10,xtitle='Log!d10!nM!dej,cold!n/M!dstream!n', $  
        ytitle='dN / dLog!d10!n(M!dej,cold!n/M!dstream!n) / Mpc!u3!n',charsize=1.5,xr=[-4.0,4.0],/xs,yr=[1.e-6,1.0],/ys


  ; mout  as a fraction of mass in stars
  frac_mejcold = alog10((gbk.mejcold + 1.0)/(gbk.mstar + 1.0))
  mmin = -4.0
  mmax = 4.0
  nbins = 100
  jeje_histo, frac_mejcold, mmin, mmax, nbins, hstream
  plot_io,hstream.x,hstream.dn/hstream.dx/run.volumempc,psym=10,xtitle='Log!d10!nM!dej,cold!n/M!dstar!n', $  
        ytitle='dN / dLog!d10!n(M!dej,cold!n/M!dstar!n) / Mpc!u3!n',charsize=1.5,xr=[-4.0,4.0],/xs,yr=[1.e-6,1.0],/ys


  ; fraction of hot gas as a function of halo mass
  baryon_fraction = 0.044/0.24
  frac_hot = (gbk.mhotgaz+ 1.0)/(gbk.mhotgaz + gbk.mcoldstream)
  plot,alog10(gbk.mfof), frac_hot,psym=4,xtitle='Log!d10!nM!dFOF!n', $  
        ytitle='M!dhot!n/(M!dhot!n+M!dcold!n)',charsize=1.5,xr=[9.0, 13.0],/xs,yr=[0.0,1.0],/ys


  ; mass bins (in log) for halo mass or stellar mass
  mmin  = 7.
  mmax  = 14.
  nbins = 40
  dm    = (mmax - mmin) / float(nbins)
  mstarmeds = replicate({m:0.0,hmacc:0.0,maccstream:0.0,maccfountain:0.0,mejcold:0.0,mejhot:0.0,mejout:0.0,mstarform:0.0},nbins)
  mfofmeds  = mstarmeds
  lmstar = alog10(gbk.mstar)
  lmfof  = alog10(gbk.mfof)
  ; compute median of everything 
  for imass = 0,nbins-1 do begin 
     m1 = mmin + imass * dm
     m2 = m1 + dm
     m  = 0.5*(m1+m2)
     mstarmeds(imass).m = m
     mfofmeds(imass).m  = m
     ; stellar mass
     ii = where(lmstar lt m2 and lmstar ge m1,ni)
     if ni ne 0 then begin 
        mstarmeds(imass).hmacc        = jeje_median(gbk(ii).hmacc)
        mstarmeds(imass).maccstream   = jeje_median(gbk(ii).maccstream)
        mstarmeds(imass).maccfountain = jeje_median(gbk(ii).maccfountain)
        mstarmeds(imass).mejcold      = jeje_median(gbk(ii).mejcold)
        mstarmeds(imass).mejhot       = jeje_median(gbk(ii).mejhot)
        mstarmeds(imass).mejout       = jeje_median(gbk(ii).mejout)
        mstarmeds(imass).mstarform    = jeje_median(gbk(ii).mstarform)
     endif 
     ; halo mass
     ii = where(lmfof lt m2 and lmfof ge m1,ni)
     if ni ne 0 then begin 
        mfofmeds(imass).hmacc        = jeje_median(gbk(ii).hmacc)
        mfofmeds(imass).maccstream   = jeje_median(gbk(ii).maccstream)
        mfofmeds(imass).maccfountain = jeje_median(gbk(ii).maccfountain)
        mfofmeds(imass).mejcold      = jeje_median(gbk(ii).mejcold)
        mfofmeds(imass).mejhot       = jeje_median(gbk(ii).mejhot)
        mfofmeds(imass).mejout       = jeje_median(gbk(ii).mejout)
        mfofmeds(imass).mstarform    = jeje_median(gbk(ii).mstarform)
     endif 
  endfor

  ; plot everything :) 
  loadct,39
  plot_io,mstarmeds.m,mstarmeds.hmacc,xtitle='Log!d10!n Mstar',ytitle='M',thick=3,yr=[1.d7,1.d12],/ys,xr=[7,11.5],/xs
  oplot,mstarmeds.m,mstarmeds.maccstream,thick=3,col=20
  oplot,mstarmeds.m,mstarmeds.maccfountain,thick=3,col=100
  oplot,mstarmeds.m,mstarmeds.mejcold,thick=3,col=150
  oplot,mstarmeds.m,mstarmeds.mejhot,thick=3,col=200
  oplot,mstarmeds.m,mstarmeds.mejout,thick=3,col=220
  oplot,mstarmeds.m,mstarmeds.mstarform,thick=3,linestyle=2
  oplot,mstarmeds.m,mstarmeds.maccstream+mstarmeds.maccfountain,thick=3,linestyle=1

  plot_io,mstarmeds.m,mstarmeds.hmacc/dt,xtitle='Log!d10!n Mstar',ytitle='[Msun/yr]',thick=3,yr=[0.01,1.d3],/ys,xr=[7,11.5],/xs
  oplot,mstarmeds.m,mstarmeds.maccstream/dt,thick=3,col=20
  oplot,mstarmeds.m,mstarmeds.maccfountain/dt,thick=3,col=100
  oplot,mstarmeds.m,mstarmeds.mejcold/dt,thick=3,col=150
  oplot,mstarmeds.m,mstarmeds.mejhot/dt,thick=3,col=200
  oplot,mstarmeds.m,mstarmeds.mejout/dt,thick=3,col=220
  oplot,mstarmeds.m,mstarmeds.mstarform/dt,thick=3,linestyle=2
  oplot,mstarmeds.m,(mstarmeds.maccstream+mstarmeds.maccfountain)/dt,thick=3,linestyle=1


  plot_io,mfofmeds.m,mfofmeds.hmacc,xtitle='Log!d10!n Mfof',ytitle='M',thick=3,yr=[1.d7,1.d12],/ys,xr=[10,13.5],/xs
  oplot,mfofmeds.m,mfofmeds.maccstream,thick=3,col=20
  oplot,mfofmeds.m,mfofmeds.maccfountain,thick=3,col=100
  oplot,mfofmeds.m,mfofmeds.mejcold,thick=3,col=150
  oplot,mfofmeds.m,mfofmeds.mejhot,thick=3,col=200
  oplot,mfofmeds.m,mfofmeds.mejout,thick=3,col=220
  oplot,mfofmeds.m,mfofmeds.mstarform,thick=3,linestyle=2
  oplot,mfofmeds.m,mfofmeds.maccstream+mfofmeds.maccfountain,thick=3,linestyle=1
  
  plot_io,mfofmeds.m,mfofmeds.hmacc/dt,xtitle='Log!d10!n Mfof',ytitle='[Msun/yr]',thick=3,yr=[0.1,1.d3],/ys,xr=[10,13.5],/xs
  oplot,mfofmeds.m,mfofmeds.maccstream/dt,thick=3,col=20
  oplot,mfofmeds.m,mfofmeds.maccfountain/dt,thick=3,col=100
  oplot,mfofmeds.m,mfofmeds.mejcold/dt,thick=3,col=150
  oplot,mfofmeds.m,mfofmeds.mejhot/dt,thick=3,col=200
  oplot,mfofmeds.m,mfofmeds.mejout/dt,thick=3,col=220
  oplot,mfofmeds.m,mfofmeds.mstarform/dt,thick=3,linestyle=2
  oplot,mfofmeds.m,(mfofmeds.maccstream+mfofmeds.maccfountain)/dt,thick=3,linestyle=1

  ; more plots
  plot_io,lmstar,gbk.mstarform/dt,psym=3,xtitle='Log!d10!n Mstar',ytitle='[Msun/yr]',yr=[0.01,1.d3],/ys,xr=[7,11.5],/xs
  oplot,mstarmeds.m,mstarmeds.mstarform/dt,thick=3,col=250


end

