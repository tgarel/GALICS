pro plot_obs_data,thing,hubble,datadir=datadir,symb=symb

; hubble is H/100
  
  if not keyword_set(symb) then symb = 7

 
  if not keyword_set(datadir) then begin 
     spawn,'echo $GALICS_PATH',path
     datadir = strtrim(path,2)+'/idl/ObsData/'
  endif
  
  if thing eq 'SDSS_u_LF' then begin 
     read_blanton_lf,'u',lf
     oplot,lf.absm,lf.phi,psym=symb,symsize=0.3 ; m-5log(h) vs. N/[Mpc/h]^3/mag
  endif
  if thing eq 'SDSS_g_LF' then begin 
     read_blanton_lf,'g',lf
     oplot,lf.absm,lf.phi,psym=symb,symsize=0.3 ; m-5log(h) vs. N/[Mpc/h]^3/mag
  endif
  if thing eq 'SDSS_r_LF' then begin 
     read_blanton_lf,'r',lf
     oplot,lf.absm,lf.phi,psym=symb,symsize=0.3 ; m-5log(h) vs. N/[Mpc/h]^3/mag
  endif
  if thing eq 'SDSS_i_LF' then begin 
     read_blanton_lf,'i',lf
     oplot,lf.absm,lf.phi,psym=symb,symsize=0.3 ; m-5log(h) vs. N/[Mpc/h]^3/mag
  endif
  if thing eq 'SDSS_z_LF' then begin 
     read_blanton_lf,'z',lf
     oplot,lf.absm,lf.phi,psym=symb,symsize=0.3 ; m-5log(h) vs. N/[Mpc/h]^3/mag
  endif
  
  if thing eq 'BDP_DR3MF' then begin 
     restore,datadir+'DR3_MF_data_BdP.sav'
     oplot,mass_x,freq_y,psym=symb,thick=4 ; LogMstar vs N/Mpc^3/LogMstar
  endif
  
  if thing eq 'Babbedge_MIPS_24mic_z0' then begin
     bab = read_babbedge('MIPS_24mic', 0)
     oplot, bab.loglum, 10^bab.logLF, psym=symb,thick=2 ; log10(\nu L_{\nu}/(h^-2 L_sun)) vs (dN/d(log10(L)) / h^3 Mpc^-3)
     errplot, bab.loglum, 10^bab.logLFlow, 10^bab.logLFhigh, thick=2
  endif
  
  if thing eq 'Babbedge_MIPS_24mic_z1' then begin
     bab = read_babbedge('MIPS_24mic', 1)
     oplot, bab.loglum, 10^bab.logLF, psym=symb,thick=2 ; log10(\nu L_{\nu}/(h^-2 L_sun)) vs (dN/d(log10(L)) / h^3 Mpc^-3)
     errplot, bab.loglum, 10^bab.logLFlow, 10^bab.logLFhigh, thick=2
  endif
  
  if thing eq 'Babbedge_MIPS_24mic_z2' then begin
     bab = read_babbedge('MIPS_24mic', 2)
     oplot, bab.loglum, 10^bab.logLF, psym=symb,thick=2 ; log10(\nu L_{\nu}/(h^-2 L_sun)) vs (dN/d(log10(L)) / h^3 Mpc^-3)
     errplot, bab.loglum, 10^bab.logLFlow, 10^bab.logLFhigh, thick=2
  endif

  if thing eq 'Babbedge_MIPS_24mic_z3' then begin
     bab = read_babbedge('MIPS_24mic', 3)
     oplot, bab.loglum, 10^bab.logLF, psym=symb,thick=2 ; log10(\nu L_{\nu}/(h^-2 L_sun)) vs (dN/d(log10(L)) / h^3 Mpc^-3)
     errplot, bab.loglum, 10^bab.logLFlow, 10^bab.logLFhigh, thick=2
  endif
  
  if thing eq 'Caputi_IRbol_z1' then begin
     cap = read_caputi('IRbol', 1)
     oplot, cap.loglum, 10^cap.logLF, psym=symb,thick=2 ; log10(\nu L_{bol}/(h^-2 L_sun)) vs (dN/d(log10(L)) / h^3 Mpc^-3)
     errplot, cap.loglum, 10^cap.logLFlow, 10^cap.logLFhigh, thick=2
  endif


  if thing eq 'Caputi_IRbol_z2' then begin
     cap = read_caputi('IRbol', 2)
     oplot, cap.loglum, 10^cap.logLF, psym=symb,thick=2 ; log10(\nu L_{bol}/(h^-2 L_sun)) vs (dN/d(log10(L)) / h^3 Mpc^-3)
     errplot, cap.loglum, 10^cap.logLFlow, 10^cap.logLFhigh, thick=2
  endif

  if thing eq 'Drory_Mstar_z2' then begin
     drd = read_drory('FDF', 2)
     oplot, drd.loglum, 10^drd.logLF, psym=symb,thick=2 ; log10(M_{star}/(M_sun)) vs (dN/d(log10(L)) / Mpc^-3)
     errplot, drd.loglum, 10^drd.logLFlow, 10^drd.logLFhigh, thick=2

     drdd = read_drory('GOODS', 2)
     oplot, drdd.loglum, 10^drdd.logLF, psym=symb+1,thick=2 ; log10(M_{star}/(M_sun)) vs (dN/d(log10(L)) / Mpc^-3)
     errplot, drdd.loglum, 10^drdd.logLFlow, 10^drdd.logLFhigh, thick=2
  endif

  if thing eq 'Drory_Mstar_z3' then begin
     drd = read_drory('FDF', 3)
     oplot, drd.loglum, 10^drd.logLF, psym=symb,thick=2 ; log10(M_{star}/(M_sun)) vs (dN/d(log10(L)) / Mpc^-3)
     errplot, drd.loglum, 10^drd.logLFlow, 10^drd.logLFhigh, thick=2

     drdd = read_drory('GOODS', 3)
     oplot, drdd.loglum, 10^drdd.logLF, psym=symb+1,thick=2 ; log10(M_{star}/(M_sun)) vs (dN/d(log10(L)) / Mpc^-3)
     errplot, drdd.loglum, 10^drdd.logLFlow, 10^drdd.logLFhigh, thick=2
  endif

  if thing eq 'Drory_Mstar_z4' then begin
     drd = read_drory('FDF', 4)
     oplot, drd.loglum, 10^drd.logLF, psym=symb,thick=2 ; log10(M_{star}/(M_sun)) vs (dN/d(log10(L)) / Mpc^-3)
     errplot, drd.loglum, 10^drd.logLFlow, 10^drd.logLFhigh, thick=2

     drdd = read_drory('GOODS', 4)
     oplot, drdd.loglum, 10^drdd.logLF, psym=symb+1,thick=2 ; log10(M_{star}/(M_sun)) vs (dN/d(log10(L)) / Mpc^-3)
     errplot, drdd.loglum, 10^drdd.logLFlow, 10^drdd.logLFhigh, thick=2
  endif



  if thing eq 'HI_MF' then begin 
     read_zwaan,zwa
     ; convert Msun/h_75^2 to Msun and Mpc^3/h_75^3 to Mpc^3
     zwa.logmgas  = zwa.logmgas - 2.*alog10(hubble * 100./75.)
     zwa.logN     = zwa.logN + 3. * alog10(hubble * 100./75.)
     zwa.logNlow  = zwa.logNlow + 3. * alog10(hubble * 100./75.)
     zwa.logNhigh = zwa.logNhigh + 3. * alog10(hubble * 100./75.)
     oplot,zwa.logmgas,zwa.logN,psym=symb,thick=2
     errplot,zwa.logmgas,zwa.logNlow,zwa.logNhigh,thick=2 
  endif
  
  if (thing eq 'gab04_2800_z3') or (thing eq 'gab04_2800_z1') or (thing eq 'gab04_2800_z2') or (thing eq 'gab04_2800_z2_3') or (thing eq 'gab04_2800_z4') or (thing eq 'gab04_2800_z5_6') or $
     (thing eq 'gab04_1500_z3') or (thing eq 'gab04_1500_z1') or (thing eq 'gab04_1500_z2') or (thing eq 'gab04_1500_z2_3') or (thing eq 'gab04_1500_z4') or (thing eq 'gab04_1500_z5_6') or $
     (thing eq 'gab04_g_z3') or (thing eq 'gab04_g_z1') or (thing eq 'gab04_g_z2') or (thing eq 'gab04_g_z2_3') or (thing eq 'gab04_g_z4') or (thing eq 'gab04_g_z5_6') or $
     (thing eq 'gab04_u_z3') or (thing eq 'gab04_u_z1') or (thing eq 'gab04_u_z2') or (thing eq 'gab04_u_z2_3') or (thing eq 'gab04_u_z4') or (thing eq 'gab04_u_z5_6') or $
     (thing eq 'gab04_B_z3') or (thing eq 'gab04_B_z1') or (thing eq 'gab04_B_z2') or (thing eq 'gab04_B_z2_3') or (thing eq 'gab04_B_z4') or (thing eq 'gab04_B_z5_6') $
  then begin 
     read_gabash,thing,gab,datadir=datadir
     oplot,gab.m,gab.n,psym=symb
     errplot,gab.m,gab.n+gab.l,gab.n+gab.u
  endif
  
  if (thing eq 'bouwens07_1600A_z4') or (thing eq 'bouwens07_1600A_z5') or (thing eq 'bouwens07_1600A_z6') then begin

     read_bouwens07,thing,bouw,datadir=datadir
     oplot,bouw.m,10.^bouw.n,psym=symb
     errplot,bouw.m,10.^(bouw.n+bouw.l),10.^(bouw.n+bouw.u)

  endif

  if (thing eq 'reddy08_1700A_z2') or (thing eq 'reddy08_1700A_z3') then begin

     read_reddy08,thing,reddy,datadir=datadir
     oplot,reddy.m,reddy.n,psym=symb
     errplot,reddy.m,reddy.n+reddy.l,reddy.n+reddy.u
     ;errplot,reddy.m+reddy.left,reddy.m+reddy.right,reddy.n
     ; there are horizontal error bars
  endif
     
  if (thing eq 'mclure_z5') or (thing eq 'mclure_z6') then begin

     read_mclure,thing,mc,datadir=datadir
     oplot,mc.m,10.^mc.n,psym=symb
     errplot,mc.m,10.^(mc.n+mc.l),10.^(mc.n+mc.u)
     
  endif

  if (thing eq 'iwata_z5') then begin

     read_iwata,thing,iwa,datadir=datadir
     oplot,iwa.m,10.^iwa.n,psym=symb
     errplot,iwa.m,10.^(iwa.n+iwa.l),10.^(iwa.n+iwa.u)
     
  endif
     

end
