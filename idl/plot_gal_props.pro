pro plot_gal_props,run,nosubdir=nosubdir


z        = 3.
timestep = z2ts(z,run.snapzdir)
if not keyword_set(nosubdir) then begin 
   for i = run.dirmin,run.dirmax do begin 
      outputdir= run.outputdir + strtrim(i,2) + '/'
      if i eq run.dirmin then begin
         read_gal_results,outputdir,timestep,galres,grinfo
         read_gal_sfr,outputdir,timestep,sfr
      endif else begin 
         read_gal_results,outputdir,timestep,galrest,grinfo
         galres = [galres,galrest]
         read_gal_sfr,outputdir,timestep,sfrt
         sfr    = [sfr,sfrt]
      endelse
   endfor
endif else begin 
   read_gal_results,run.outputdir,timestep,galres,grinfo   
   read_gal_sfr,outputdir,timestep,sfr
endelse


plot_oo,galres.disc_rgal*1000.,galres.burst_rgal*1000.,psym=3,xtitle='rdisc',ytitle='rburst'

plot_oo,sfr.disc_sfr100,sfr.burst_sfr100,psym=3,xtitle='SFR100 disc',ytitle='SFR100 burst'
oplot,[0.01,1000.],[0.01,1000.]



end

