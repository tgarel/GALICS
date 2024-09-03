pro read_halo_baryondec,outputdir,ts,hbd

  filename=outputdir + '/halo_baryondec'+string(ts,format='(i3.3)')+'.dat'
  spawn,'cat '+filename+' | wc -l',nh
  hbd = replicate({nbgal:0L,coldgas:0.0,mhot:0.0,mcoldgaz:0.0,mstar:0.0,mism:0.0,mwind:0.0,macc:0.0,mfof:0.0,mvir:0.0},nh)
  openr,11,filename
  readf,11,hbd
  close,11
end


