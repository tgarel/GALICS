pro read_sfr_tabs, outputdir, timestep, time, deltat, sfrtabs, timetab,mettab,comp, dual_imf=dual_imf

  if (comp eq 'disc') then begin ; read discs 
     if not keyword_set(dual_imf) then $
        filename = outputdir + '/disc_sfrtabs.' + string(timestep,format='(i3.3)') $ 
     else $ 
        filename = outputdir + '/disc_sfrtabs_imf2.' + string(timestep,format='(i3.3)')
  endif
  if (comp eq 'bulge') then begin ; read bulges
     if not keyword_set(dual_imf) then $ 
        filename = outputdir + '/bulge_sfrtabs.' + string(timestep,format='(i3.3)') $ 
     else $
        filename = outputdir + '/bulge_sfrtabs_imf2.' + string(timestep,format='(i3.3)')
  endif
  if (comp eq 'burst') then begin ; read bursts
     if not keyword_set(dual_imf) then $ 
        filename = outputdir + '/burst_sfrtabs.' + string(timestep,format='(i3.3)') $ 
     else $ 
        filename = outputdir + '/burst_sfrtabs_imf2.' + string(timestep,format='(i3.3)')
  endif
  
  ; open file and read header
  openr,11,filename
  r = (ng = (nage = (nmet = 0L)))
  readu,11,r,ng,nage,nmet,r
  timetab = dblarr(nage)
  mettab  = dblarr(nmet)
  readu,11,r,timetab,r
  readu,11,r,mettab,r

  ; compute delta-time over which SFR is cumulated in each bin
  deltat = shift(timetab,-1) - timetab ; wrong for last bin, which is 20Gyr anyway so doesn't matter... 

  ; allocate comp with oversize arrays
  sfrst = {sfrtot:0.0d0, sfrtab:dblarr(nage,nmet)}

  ; read the data
  sfrtabs = replicate(sfrst,ng)
  for ig = 0L, ng - 1L do begin 
     sfrtot = 0.0d0
     n1 = (n2 = 0L)
     readu,11,r,sfrtot,n1,n2,r
     if sfrtot gt 0.0d0 then begin 
        tab = dblarr(n1,n2)
        readu,11,r,tab,r
        sfrtabs(ig).sfrtot = sfrtot
        sfrtabs(ig).sfrtab(0:n1-1,0:n2-1) = tab
     endif
  endfor

  close,11

end
