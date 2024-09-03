function expdisc, rd, rmin, rmax, nbins 

; return the cumulative surface brightness of an exponential disc
; S(R) = 2!pi * \int_0^R \Sigma_0 exp(-r/rd) * r dr
;      = 2 pi rd^2 Sigma_0 * [1 -exp(-R/rd)*(1+R/rd)]
; NB: i normalise to the total luminosity, which is 2*pi*rd^2*Sigma_0

if rmin lt 0. then stop,'rmin should be >= 0 in funtion expdisc'
r = (dindgen(nbins)+0.5) / float(nbins) * (rmax - rmin) + rmin

result = (1.d0 - exp(-r/rd)*(1.d0+r/rd))

expdisc = {r:r, s:result}
return,expdisc

end
