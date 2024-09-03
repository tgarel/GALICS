pro read_gal_nphot, outputdir, timestep, gal_N, gal_N_info

  file = outputdir + 'gal_Nion_phot.'+string(timestep,format='(i3.3)')
  openr,11,file

;;   disc_Nphot_tab  = dblarr(99999)
;;   bulge_Nphot_tab = dblarr(99999)
;;   burst_Nphot_tab = dblarr(99999)

;;   disc_Nphot_tab(*)  = -888.0d0
;;   bulge_Nphot_tab(*) = -888.0d0
;;   burst_Nphot_tab(*) = -888.0d0

  disc_Nphot  = 0.0d0
  bulge_Nphot = 0.0d0
  burst_Nphot = 0.0d0

;;   i = 0LL
;;   timestep = 0L
  ngal = 0LL

  readf,11,format='(i3,1x,i6)',st,ngal

  if ngal gt 0 then begin
  
     gal_N_info  = {st:st,ngal:ngal}
     gal_N       = replicate({disc_N:0.0d0,bulge_N:0.0d0, burst_N:0.0d0},ngal)
     
     for igal = 0LL, ngal - 1LL do begin 
        readf,11,format='(e14.6,1x,e14.6,1x,e14.6)', disc_Nphot,bulge_Nphot,burst_Nphot
        
        gal_N(igal).disc_N = disc_Nphot
        gal_N(igal).bulge_N = bulge_Nphot
        gal_N(igal).burst_N = burst_Nphot
     endfor
  endif else begin
     print,'no galaxy'
     gal_N = -1
  endelse

     ;; while not eof(11) do begin
;;      readf,11,disc_Nphot,bulge_Nphot,burst_Nphot
;;      disc_Nphot_tab(i) = disc_Nphot
;;      bulge_Nphot_tab(i) = bulge_Nphot
;;      burst_Nphot_tab(i) = burst_Nphot
;;      i = i + 1
;;   endwhile
     
;;   close,11
     
;;   print,disc_Nphot_tab
;;   w = where(disc_Nphot_tab ne -888.0d0,nw)
;;   print,nw
     
;;   disc_Nphot_tab  = disc_Nphot_tab(w)
;;   bulge_Nphot_tab = bulge_Nphot_tab(w)
;;   burst_Nphot_tab = burst_Nphot_tab(w)
     
                                ;readf,11,gal_N
     close,11
     
  end
  
