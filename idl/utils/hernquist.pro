function hernquist, a, rmin, rmax, nbins, petrosian=petrosian

; returns the cumulative surface brightness S(R) (Hernquist 1990, Eq. 37),
; normalized to total luminosity
; 
; note that the Hernquist profile reaches 90% of total luminosity at r~15*a,
; which is far out! This gives a concentration c=R90/R50=8, which is larger
; than a typical de Vaucouleur's law (where c ~ 3.3). 

; INPUTS : 
; - a : scale-length of the profile
; - rmin,rmax,nbins : define the sampling of the profile  
;
; KEYWORD : 
; - petrosian : set this to renormalise the flux to 80% of the total
;   (i.e. divide profile by 0.8). This mimicks petrosian bias in the
;   estimate of hernquist-type profiles.

; sample radii from r=0 to 10 scale-lengths
if rmin lt 0. then stop,'rmin should be >= 0 in funtion hernquist'
r  = (dindgen(nbins)+0.5) / float(nbins) * (rmax - rmin) + rmin ; avoid zero ... 
s  = r/a
s2 = s^2

x = dblarr(nbins)

ii = where(s lt 1.0d0 and s gt 0.0d0,ni)
if ni ne 0 then x(ii) = alog((1.d0+sqrt(1.d0-s2(ii)))/s(ii)) / sqrt(1.d0 -s2(ii))

ii = where(s gt 1.0d0,ni) 
if ni ne 0 then x(ii) = acos(1/s(ii)) / sqrt(s2(ii)-1.d0)

ii = where(s eq 1.0d0, ni) 
if ni ne 0 then x(ii) = 1.d0

result = dblarr(nbins)
if ni ne 0 then begin 
   jj = where(s ne 1.0d0)
   result(jj) = s(jj)^2 * (x(jj) - 1.d0) / (1.d0 - s(jj)^2)
   result(ii) = 1.d0/3.d0
endif else begin 
   result = s^2 * (x - 1.d0)/ (1.d0 - s^2)
endelse

if keyword_set(petrosian) then $
   result = result / 0.8d0

; build struct to return
hern = {r:r,s:result}
return, hern

end
