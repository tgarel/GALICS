function define_structure,struct_name,renumbering=renumbering

; structure of the gal_result files
; ---------------------------------
if struct_name eq 'gal_results' then begin 
   struct = { $
            disturb          : 0.0, $
            r_orbit          : 0.0, $
            hno              : 0L,  $
            cvel             : 0.0, $
            inclination      : 0.0, $
            nb_merg          : 0L,  $
            tbirth           : 0.0, $
            mvir             : 0.0, $
            rvir             : 0.0, $
            mfof             : 0.0, $
            mwind_gal        : 0.0, $
            mwind_gal_z      : 0.0, $
            disc_mgal        : 0.0, $
            disc_mcold       : 0.0, $
            disc_sfr         : 0.0, $
            disc_minstar     : 0.0, $
            disc_mcoldz      : 0.0, $
            disc_rgal        : 0.0, $         
            disc_tdyn        : 0.0, $
            disc_speed       : 0.0, $
            disc_transp      : 0.0, $
            bulge_mgal       : 0.0, $
            bulge_mcold      : 0.0, $
            bulge_sfr        : 0.0, $
            bulge_minstar    : 0.0, $
            bulge_mcoldz     : 0.0, $
            bulge_rgal       : 0.0, $         
            bulge_tdyn       : 0.0, $
            bulge_speed      : 0.0, $
            bulge_transp     : 0.0, $
            burst_mgal       : 0.0, $
            burst_mcold      : 0.0, $
            burst_sfr        : 0.0, $
            burst_minstar    : 0.0, $
            burst_mcoldz     : 0.0, $
            burst_rgal       : 0.0, $         
            burst_tdyn       : 0.0, $
            burst_speed      : 0.0, $
            burst_transp     : 0.0 $
            }
endif   


if struct_name eq 'bol_lum' then begin 
   struct = {lum:0.0, irlum:0.0}
endif


if struct_name eq 'cores' then begin 
   struct = {r:0.0,mass:0.0,cooling_frac:0.0}
endif


if struct_name eq 'galaxyids' then begin 
   struct = {BushID            : lon64arr(1), $
             TreeID            : lon64arr(1), $
             GalaxyID          : lon64arr(1), $
             DescendantID      : lon64arr(1), $
             FirstProgenitorID : lon64arr(1), $
             NextProgenitorID  : lon64arr(1), $
             LastProgenitorID  : lon64arr(1), $
             HaloID            : lon64arr(1)}
endif


if struct_name eq 'haloids' then begin 
   struct = {BushID            : lon64arr(1), $
             TreeID            : lon64arr(1), $
             HaloID            : lon64arr(1), $
             DescendantID      : lon64arr(1), $
             FirstProgenitorID : lon64arr(1), $
             NextProgenitorID  : lon64arr(1), $
             LastProgenitorID  : lon64arr(1)}
endif


if struct_name eq 'gal_sfr' then begin 
   struct = {disc_sfr1:0.0,bulge_sfr1:0.0,burst_sfr1:0.0,$
             disc_sfr10:0.0,bulge_sfr10:0.0,burst_sfr10:0.0,$
             disc_sfr100:0.0,bulge_sfr100:0.0,burst_sfr100:0.0}
endif


if struct_name eq 'halo_dm_results' then begin 
   struct = {x:0.0 ,y:0.0 ,z:0.0,  $
             vx:0.0,vy:0.0,vz:0.0, $
             Lx:0.0,Ly:0.0,Lz:0.0, $
             ek:0.0,ep:0.0,et:0.0,spin:0.0}
endif


if struct_name eq 'halo_gas_results' then begin 
   if keyword_set(renumbering) then begin 
      struct = {my_number   :0L, $  
                mfof        : 0.0, $
                mvir        : 0.0, $
                mcoldgaz    : 0.0, $    
                mtemp       : 0.0, $
                mhotgaz     : 0.0, $
                mhotz       : 0.0, $
                mcoldz      : 0.0, $
                rfof        : 0.0, $
                rvir        : 0.0, $
                tvir        : 0.0, $
                fthtemp     : 0.0, $
                nbgal       : 0L,  $
                r_c         : 0.0, $
                rho_0       : 0.0, $
                mstar       : 0.0, $
                metsinstars : 0.0}
   endif else begin 
      struct = {mfof        : 0.0, $
                mvir        : 0.0, $
                mcoldgaz    : 0.0, $    
                mtemp       : 0.0, $
                mhotgaz     : 0.0, $
                mhotz       : 0.0, $
                mcoldz      : 0.0, $
                rfof        : 0.0, $
                rvir        : 0.0, $
                tvir        : 0.0, $
                fthtemp     : 0.0, $
                nbgal       : 0L,  $
                r_c         : 0.0, $
                rho_0       : 0.0, $
                mstar       : 0.0, $
                metsinstars : 0.0}
   endelse
endif

if struct_name eq 'halo_reservoir' then begin 
   struct = {accretion:0.0,mhotgaz:0.0,mhotmet:0.0}
endif

return,struct

end
