pro read_centralprops,outputdir,ts,cp

  filename=outputdir + '/centralprops'+string(ts,format='(i3.3)')+'.dat'
  openr,11,filename
  readf,11,ng
  cp = replicate({rvir:0.0,mvir:0.0,cvel:0.0,rfof:0.0,mfof:0.0,sfr100:0.0,sfr10:0.0,mstar:0.0,mgas:0.0,macc100:0.0},ng)
  readf,11,cp
  close,11
end

