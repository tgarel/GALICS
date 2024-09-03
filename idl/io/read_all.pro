pro read_all, outputdir, timestep

; read gal_results.xxx files
read_gal_results, outputdir, timestep, galres

; read bolometric luminosities
;read_bol_lum, outputdir, timestep, bol

; read DM cores props 
read_cores, outputdir, timestep, cores

; read galaxy IDs and tree information 
;read_galaxyids, outputdir, timestep, galid

; read instantaneous SFR 
;read_gal_sfr, outputdir, timestep, sfr

; read halo props (DM only)
;read_halo_dm_results, outputdir, timestep, hdmres

; read halo gas results 
;read_halo_gas_results, outputdir, timestep, hgres,/renumbering

; read halo IDs 
read_haloids, outputdir, timestep, haloid

; read halo reservoir content 
read_halo_reservoir, outputdir, timestep, hreserv


; READ MAGNITUDES ... 

; read the magnitudes (rest-frame)
read_gal_mags, outputdir, timestep, filters, mags, /total

; read the magnitudes (rest-frame with no extinction)
read_gal_mags, outputdir, timestep, filters, mags_noext, /noext,/total


end
