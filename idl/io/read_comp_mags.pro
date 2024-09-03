pro read_comp_mags, outputdir, timestep, comp, compmags, noext=noext

if not keyword_set(noext) then $
  file = outputdir + 'rf_mags_' + comp + '.' + string(timestep,format='(i3.3)') $
else $
  file = outputdir + 'rf_noext_mags_' + comp + '.' + string(timestep,format='(i3.3)')

nfilters = (ngalaxies = 1L)
openr,11,file
readf,11,nfilters,ngalaxies
if ngalaxies ne 0 then begin
compmags = fltarr(nfilters,ngalaxies,/nozero)
readf,11,compmags
endif else begin
compmags = 888
endelse

close,11

end
