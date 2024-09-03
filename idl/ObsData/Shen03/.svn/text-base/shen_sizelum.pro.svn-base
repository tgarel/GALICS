function shen_sizelum, mag, rmin, rmax, rbins

; returns Shen's fits to distributions of galaxy sizes in different
; magnitude bins (from Shen et al, 2003)
; NB: this is a fit of fits ... Fig. 4 show that there is a lot of scatter
; (i.e. the expressions used below are not very accurate in every bin). 
; -> use with caution.

;
; INPUTS : 
; - mag : central mag of mag bin of width 0.5
; - rmin,rmax,rbins : min max radii to sample with rbins

; early-type galaxies average size : 
etg_logr = -0.4 * 0.6 * mag -4.63
etg_lnr  = etg_logr * 2.30259

; late-type galaxies average size :
ltg_logr = -0.4*0.21*mag + (0.53-0.21)*alog10(1.+10.^(-0.4*(mag+20.52))) - 1.31
ltg_lnr  = ltg_logr * 2.30259

; dispersion (for both early and late types):
sigma = 0.25 + (0.48 - 0.25) / (1. + 10.^(-0.8*(mag+20.52)))  ; NB: this is sigma_ln(r)

; sample the size distribution 
lx   = (findgen(rbins)+0.5d0) / float(rbins) * (alog(rmax)-alog(rmin)) + alog(rmin) ; this is ln(r)...
x    = exp(lx) 
dlx  = lx(2)-lx(1)

if (mag le -18.75 and mag ge -23.75) then $
   netg = shen_lognormal(lx,etg_lnr,sigma) $
else $ 
   netg = x * 0.

if (mag lt -16.75 and mag ge -23.75) then $
   nltg = shen_lognormal(lx,ltg_lnr,sigma) $
else $
   nltg = x * 0.

; normalise (so that area in Log10 space is one) 
totnetg = total(netg,/double)
if totnetg gt 0. then netg = netg / totnetg * 2.3 / dlx
totnltg = total(nltg,/double)
if totnltg gt 0. then nltg = nltg / totnltg * 2.3 / dlx

sizedist={rmin:rmin,rmax:rmax,rbins:rbins,x:x,netg:netg,nltg:nltg}
return,sizedist

end
