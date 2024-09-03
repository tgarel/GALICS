pro read_gal_mags, outputdir, timestep, filters, mags, noext=noext, total=total,face_on=face_on

; KEYWORDS: 
; - noext to true to read mags without dust extinction
; - total to include total (sum over components) mag of each gal in the
;   structure mags ... 
; - face_on to read face-on disc mags


; read filter list
file = outputdir + 'filters.dat'
openr,11,file
readf,11,nfilters
filters = strarr(nfilters)
readf,11,filters
close,11
filters = strtrim(filters,2)

; read disc magnitudes
if keyword_set(face_on) then $
   read_comp_mags, outputdir, timestep, 'face', discmags, noext=noext $
else $
   read_comp_mags, outputdir, timestep, 'disc', discmags, noext=noext

; read bulge/burst magnitudes
read_comp_mags, outputdir, timestep, 'bulge', blgmags, noext=noext
read_comp_mags, outputdir, timestep, 'burst', bstmags, noext=noext

dsize  = size(discmags,/dimensions)
bgsize = size(blgmags,/dimensions)
btsize = size(bstmags,/dimensions)

if (dsize(0) ne 0 and bgsize(0) ne 0 and btsize(0) ne 0) then begin
mags_to_struct, filters, discmags, blgmags, bstmags, mags, total=total
endif

end
