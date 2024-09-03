pro read_reddy08,thing,reddy,datadir=datadir

  if not keyword_set(datadir) then begin 
     spawn,'echo $GALICS_PATH',path
     datadir = strtrim(path,2)+'/idl/ObsData/'
  end

  if thing eq 'reddy08_1700A_z2' then begin 
     reddy = replicate({m:0.0,n:0.0,left:0.0,right:0.0,l:0.0,u:0.0},9)
     fichier = datadir + 'reddy08/reddy08_1700A_z2.txt'
     openr,11,fichier
     readf,11,reddy
     close,11
  endif
  
  if thing eq 'reddy08_1700A_z3' then begin 
     reddy = replicate({m:0.0,n:0.0,left:0.0,right:0.0,l:0.0,u:0.0},10)
     fichier = datadir + 'reddy08/reddy08_1700A_z3.txt'
     openr,11,fichier
     readf,11,reddy
     close,11
  endif

end
