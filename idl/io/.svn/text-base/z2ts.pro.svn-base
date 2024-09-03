function z2ts,z,dir,silent=silent

; reads the file snaps_redshifts.dat which should be in the directory dir 
; uses it to find closest timestep to redshift z passed as input. 

file = dir + 'snaps_redshifts.dat'
z1 = (z2 = 0.)
ts = 0
openr,11,file
while not eof(11) do begin 
   ts = ts + 1
   readf,11,z2
   if z2 le z then begin 
      ; ts is either z2 or z1 
      d1 = abs(z1 - z)
      d2 = abs(z2 - z)
      if d1 lt d2 then begin 
         ts = ts - 1
         if not keyword_set(silent) then $
            print,'z2ts returns timestep '+string(ts)+$
                  ' at redshift '+string(Z1),' when asked for '+string(z)
      endif else begin 
         if not keyword_set(silent) then $
            print,'z2ts returns timestep '+string(ts)+$
                  ' at redshift '+string(Z2),' when asked for '+string(z)
      endelse
      close,11
      return,ts
   endif
   z1 = z2
endwhile

close,11
stop,'No redshift found ... in z2ts'

end 


      
