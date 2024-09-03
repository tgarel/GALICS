pro read_iwata,thing,iwa,datadir=datadir

  if not keyword_set(datadir) then begin 
     spawn,'echo $GALICS_PATH',path
     datadir = strtrim(path,2)+'/idl/ObsData/'
  end

  if thing eq 'iwata_z5' then begin 
     iwa = replicate({m:0.0,n:0.0,l:0.0,u:0.0},7)
     fichier = datadir + 'iwata07/iwata_z5.txt'
     openr,11,fichier
     readf,11,iwa
     close,11
  endif

end
