pro read_mclure,thing,mc,datadir=datadir

  if not keyword_set(datadir) then begin 
     spawn,'echo $GALICS_PATH',path
     datadir = strtrim(path,2)+'/idl/ObsData/'
  end


  if thing eq 'mclure_z5' then begin 
     mc = replicate({m:0.0,n:0.0,l:0.0,u:0.0},8)
     fichier = datadir + 'mclure08/mclure_z5.txt'
     openr,11,fichier
     readf,11,mc
     close,11
  endif

if thing eq 'mclure_z6' then begin 
     mc = replicate({m:0.0,n:0.0,l:0.0,u:0.0},6)
     fichier = datadir + 'mclure08/mclure_z6.txt'
     openr,11,fichier
     readf,11,mc
     close,11
  endif

end
