pro read_gabash,thing,gab,datadir=datadir

  if not keyword_set(datadir) then begin 
     spawn,'echo $GALICS_PATH',path
     datadir = strtrim(path,2)+'/idl/ObsData/'
  end

;---------- 2800A band----------------------------

  if thing eq 'gab04_2800_z3' then begin 
     gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
     fichier = datadir + 'gabasch04/2800_band/gab04_2800_z3.txt'
     openr,11,fichier
     readf,11,gab
     close,11
  endif
  
  if thing eq 'gab04_2800_z1' then begin 
     gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
     fichier = datadir + 'gabasch04/2800_band/gab04_2800_z1.txt'
     openr,11,fichier
     readf,11,gab
     close,11
  endif

 if thing eq 'gab04_2800_z2' then begin 
     gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},5)
     fichier = datadir + 'gabasch04/2800_band/gab04_2800_z2.txt'
     openr,11,fichier
     readf,11,gab
     close,11
  endif

 if thing eq 'gab04_2800_z2_3' then begin 
     gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},5)
     fichier = datadir + 'gabasch04/2800_band/gab04_2800_z2_3.txt'
     openr,11,fichier
     readf,11,gab
     close,11
  endif

 if thing eq 'gab04_2800_z4' then begin 
     gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
     fichier = datadir + 'gabasch04/2800_band/gab04_2800_z4.txt'
     openr,11,fichier
     readf,11,gab
     close,11
  endif

 if thing eq 'gab04_2800_z5_6' then begin 
     gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
     fichier = datadir + 'gabasch04/2800_band/gab04_2800_z5_6.txt'
     openr,11,fichier
     readf,11,gab
     close,11
  endif

;-----------------1500A band-----------------------------------

 if thing eq 'gab04_1500_z3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/1500_band/gab04_1500_z3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif

 if thing eq 'gab04_1500_z2' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/1500_band/gab04_1500_z2.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 
 if thing eq 'gab04_1500_z2_3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},5)
    fichier = datadir + 'gabasch04/1500_band/gab04_1500_z2_3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
  endif
 
 if thing eq 'gab04_1500_z4' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},3)
    fichier = datadir + 'gabasch04/1500_band/gab04_1500_z4.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_1500_z5_6' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/1500_band/gab04_1500_z5_6.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
;---------------------- u band -----------------------------------
 
 if thing eq 'gab04_u_z1' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
    fichier = datadir + 'gabasch04/u_band_sdss/gab04_u_z1.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_u_z2' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},5)
    fichier = datadir + 'gabasch04/u_band_sdss/gab04_u_z2.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_u_z2_3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
    fichier = datadir + 'gabasch04/u_band_sdss/gab04_u_z2_3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_u_z3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/u_band_sdss/gab04_u_z3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_u_z4' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},5)
    fichier = datadir + 'gabasch04/u_band_sdss/gab04_u_z4.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_u_z5_6' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},3)
    fichier = datadir + 'gabasch04/u_band_sdss/gab04_u_z5_6.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
;---------------------- B Johnson band -----------------------------------
 
 if thing eq 'gab04_B_z1' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
    fichier = datadir + 'gabasch04/B_Johnson/gab04_B_z1.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_B_z2' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
    fichier = datadir + 'gabasch04/B_Johnson/gab04_B_z2.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_B_z2_3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
    fichier = datadir + 'gabasch04/B_Johnson/gab04_B_z2_3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_B_z3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/B_Johnson/gab04_B_z3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_B_z4' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},5)
    fichier = datadir + 'gabasch04/B_Johnson/gab04_B_z4.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_B_z5_6' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/B_Johnson/gab04_B_z5_6.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
;---------------------- g band -----------------------------------

 if thing eq 'gab04_g_z1' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
    fichier = datadir + 'gabasch04/g_band_sdss/gab04_g_z1.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_g_z2' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
    fichier = datadir + 'gabasch04/g_band_sdss/gab04_g_z2.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_g_z2_3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},5)
    fichier = datadir + 'gabasch04/g_band_sdss/gab04_g_z2_3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_g_z3' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/g_band_sdss/gab04_g_z3.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_g_z4' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/g_band_sdss/gab04_g_z4.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 if thing eq 'gab04_g_z5_6' then begin 
    gab = replicate({m:0.0,n:0.0,l:0.0,u:0.0},4)
    fichier = datadir + 'gabasch04/g_band_sdss/gab04_g_z5_6.txt'
    openr,11,fichier
    readf,11,gab
    close,11
 endif
 
 
 
end
