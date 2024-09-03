pro read_gal_results, outputdir, timestep, galres, galres_info

; open appropriate gal_results file
file = outputdir + 'gal_results.'+string(timestep,format='(i3.3)')
openr,11,file

; read header 
st = (n_total = 0L) & aexp = ( age_univ = 0.0) 
readf,11,format='(i3,1x,i6,1x,e14.6,1x,e14.6)',st,n_total,aexp,age_univ
galres_info = {st:st,ng:n_total,aexp:aexp,age_univ:age_univ}

if n_total gt 0 then begin 
; define and allocate galres structure
   galres = define_structure('gal_results')
   galres = replicate(galres,n_total)
   
; read the rest of the file 
   ON_IOERROR, pure_disc
   for igal = 0L, n_total - 1L do begin 
      readf,11,format='(11(e14.6,1x),2(e14.6,1x),i7,1x,2(e14.6,1x),1x,i6,1x,4(e14.6,1x),18(e14.6,1x))',$
            mwind_gal,mwind_gal_z,disc_mgal,disc_mcold,disc_sfr,disc_minstar,disc_mcoldz,disc_rgal,disc_tdyn,disc_speed,disc_transp, $
            gal_disturb,gal_rorbit,gal_hno,gal_cvel,gal_incl,nb_merg,tbirth,mvir,rvir,mfof, $
            bulge_mgal,bulge_mcold,bulge_sfr,bulge_minstar,bulge_mcoldz,bulge_rgal,bulge_tdyn,bulge_speed,bulge_transp, $
            burst_mgal,burst_mcold,burst_sfr,burst_minstar,burst_mcoldz,burst_rgal,burst_tdyn,burst_speed,burst_transp
      galres(igal).bulge_mgal       = bulge_mgal   
      galres(igal).bulge_mcold      = bulge_mcold  
      galres(igal).bulge_sfr        = bulge_sfr    
      galres(igal).bulge_minstar    = bulge_minstar
      galres(igal).bulge_mcoldz     = bulge_mcoldz 
      galres(igal).bulge_rgal       = bulge_rgal   
      galres(igal).bulge_tdyn       = bulge_tdyn   
      galres(igal).bulge_speed      = bulge_speed  
      galres(igal).bulge_transp     = bulge_transp 
      galres(igal).burst_mgal       = burst_mgal   
      galres(igal).burst_mcold      = burst_mcold  
      galres(igal).burst_sfr        = burst_sfr    
      galres(igal).burst_minstar    = burst_minstar
      galres(igal).burst_mcoldz     = burst_mcoldz 
      galres(igal).burst_rgal       = burst_rgal   
      galres(igal).burst_tdyn       = burst_tdyn   
      galres(igal).burst_speed      = burst_speed  
      galres(igal).burst_transp     = burst_transp 
      pure_disc: 
      galres(igal).disturb          = gal_disturb      
      galres(igal).r_orbit          = gal_rorbit      
      galres(igal).hno              = gal_hno          
      galres(igal).cvel             = gal_cvel         
      galres(igal).inclination      = gal_incl
      galres(igal).nb_merg          = nb_merg      
      galres(igal).tbirth           = tbirth
      galres(igal).mvir             = mvir
      galres(igal).rvir             = rvir
      galres(igal).mfof             = mfof
      galres(igal).mwind_gal        = mwind_gal
      galres(igal).mwind_gal_z      = mwind_gal_z
      galres(igal).disc_mgal        = disc_mgal    
      galres(igal).disc_mcold       = disc_mcold   
      galres(igal).disc_sfr         = disc_sfr     
      galres(igal).disc_minstar     = disc_minstar 
      galres(igal).disc_mcoldz      = disc_mcoldz  
      galres(igal).disc_rgal        = disc_rgal    
      galres(igal).disc_tdyn        = disc_tdyn    
      galres(igal).disc_speed       = disc_speed   
      galres(igal).disc_transp      = disc_transp 
   endfor
endif else begin 
   print,'No galaxy found ... '
   galres = -1
endelse
; close and return
close,11

end
