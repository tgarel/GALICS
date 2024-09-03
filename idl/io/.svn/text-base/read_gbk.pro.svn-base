pro read_gbk, outputdir, timestep, gbk

  fichier = outputdir + '/gbk.'+string(timestep,'(i3.3)')
  spawn,'cat '+fichier + ' | wc -l',n
  n = long(n)
  gbk = replicate({hmacc:0.0d0,mfof:0.0d0,mcoldstream:0.0d0,mhotgaz:0.0d0,mhotz:0.0d0,maccstream:0.0d0,maccfountain:0.0d0,$
                   mejcold:0.0d0,mejhot:0.0d0,mejout:0.0d0,mstarform:0.0d0,$
                   deltat:0.0d0,mstar:0.0d0},n)
  openr,11,fichier
  readf,11,gbk
  close,11
end 
