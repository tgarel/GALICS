pro jj

snapzdir = '/simu1/data2/garel/GM-RUNS/1024-100h-1Mpc-W3/kenn_sfr/sf20_fd025_wind_corr/1/'
  ;snapzdir = '/data/garel/GM-RUNS/512-100h-1Mpc-W3/kenn_sfr/sf20_feedB002_best3/1/'
  timestep = z2ts(z,snapzdir)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for i=1L,111L do begin

     outputdir='/simu1/data2/garel/GM-RUNS/1024-100h-1Mpc-W3/kenn_sfr/sf20_fd025_wind_corr/'+strtrim(i,2)+'/'     

     if i eq 1 then begin
        read_gal_results,outputdir,timestep,S
        read_gal_sfr,outputdir,timestep,a
        read_gal_nphot,outputdir,timestep,n
        read_gal_uv1200,outputdir,timestep,uv
        read_gal_uv1200_ext,outputdir,timestep,e
     endif else begin
        read_gal_results,outputdir,timestep,S2
        read_gal_sfr,outputdir,timestep,a2
        read_gal_nphot,outputdir,timestep,n2
        read_gal_uv1200,outputdir,timestep,uv2
        read_gal_uv1200_ext,outputdir,timestep,e2
        S=[S,S2]
        a=[a,a2]
        n   = [n,n2]
        uv  = [uv,uv2]
        e   = [e,e2]

     endelse
     checkbox++
  endfor








end
