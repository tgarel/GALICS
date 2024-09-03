pro read_zwaan,zwa

  spawn,'echo $GALICS_PATH',path
  dir = strtrim(path,2)+'/idl/ObsData/Zwaan05/'
  filename = dir+'HI_Zwaan.cnt'
  openr,11,filename
  readf,11,n1,n2
  a=fltarr(n1,n2)
  readf,11,a
  close,11

  zwa = {logmgas : reform(a(0,*)),$  ; Log(Msun / h_75^2)
        logN     : reform(a(1,*)),$
        logNlow  : reform(a(2,*)),$
        logNhigh : reform(a(3,*))}
  
end
