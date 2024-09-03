pro read_bouwens07,thing,bouw,datadir=datadir

  if not keyword_set(datadir) then begin 
     spawn,'echo $GALICS_PATH',path
     datadir = strtrim(path,2)+'/idl/ObsData/'
  end

  if thing eq 'bouwens07_1600A_z4' then begin 
     bouw = replicate({m:0.0,n:0.0,l:0.0,u:0.0},14)
     fichier = datadir + 'bouwens07/bouwens07_1600A_z4.txt'
     openr,11,fichier
     readf,11,bouw
     close,11
  endif
  
  if thing eq 'bouwens07_1600A_z5' then begin 
     bouw = replicate({m:0.0,n:0.0,l:0.0,u:0.0},12)
     fichier = datadir + 'bouwens07/bouwens07_1600A_z5.txt'
     openr,11,fichier
     readf,11,bouw
     close,11
  endif

  if thing eq 'bouwens07_1600A_z6' then begin 
     bouw = replicate({m:0.0,n:0.0,l:0.0,u:0.0},8)
     fichier = datadir + 'bouwens07/bouwens07_1600A_z6.txt'
     openr,11,fichier
     readf,11,bouw
     close,11
  endif


end
