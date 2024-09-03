function read_tremonti

  spawn,'echo $GALICS_PATH',path
  dir = strtrim(path,2)+'/idl/ObsData/Tremonti04/'
  filename = dir + 'tremonti.txt'
  openr,11,filename
  pipo='' & readf,11,pipo
  tremonti = fltarr(6,28)
  readf,11,tremonti
  close,11
  mzr = {mstar:reform(tremonti(0,*)),$
         P2: reform(tremonti(1,*)),$
         p16:reform(tremonti(2,*)),$
         median:reform(tremonti(3,*)),$
         p84:reform(tremonti(4,*)),$
         p97:reform(tremonti(5,*))$
        }
  return,mzr

end
