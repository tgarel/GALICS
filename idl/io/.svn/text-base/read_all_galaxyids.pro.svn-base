pro read_all_galaxyids,outputdir,nsteps,gtree

; count the galaxies 
  ntot = 0L
  for ts = 1,nsteps do begin
     filename = outputdir + '/GalaxyIDs.'+string(ts,format='(i3.3)')
     openr,11,filename
     ng = 0L
     readf,11,ng
     close,11
     ntot = ntot + ng
  endfor

  treestruct = {BushID            : lon64arr(1), $
                TreeID            : lon64arr(1), $
                GalaxyID          : lon64arr(1), $
                DescendantID      : lon64arr(1), $
                FirstProgenitorID : lon64arr(1), $
                NextProgenitorID  : lon64arr(1), $
                LastProgenitorID  : lon64arr(1), $
                HaloID            : lon64arr(1)}
  treestruct2 = {BushID            : lon64arr(1), $
                TreeID            : lon64arr(1), $
                GalaxyID          : lon64arr(1), $
                DescendantID      : lon64arr(1), $
                FirstProgenitorID : lon64arr(1), $
                NextProgenitorID  : lon64arr(1), $
                LastProgenitorID  : lon64arr(1), $
                HaloID            : lon64arr(1), $
                ts                : 0}

  gtree = replicate(treestruct2,ntot)
  
  np = 0L
  for ts = 1,nsteps do begin
     filename = outputdir + '/GalaxyIDs.'+string(ts,format='(i3.3)')
     openr,11,filename
     ng = 0L
     readf,11,ng
     if ng gt 0 then begin 
        gtemp = replicate(treestruct,ng)
        readf,11,gtemp
        for i = 0,7 do begin 
           gtree(np:np+ng-1L).(i) = gtemp.(i)
        endfor
        gtree(np:np+ng-1L).ts  = ts
     endif
     close,11
     np = np + ng
  endfor  

end
