function baldry_coldist, mag, cmin,cmax,cbins

; A. Cattaneo & J. Blaizot (2008)
;
; returns Baldry's fit for galaxies with r-band magnitude between
; magmin and magmax. (NB: really these should be the bins that Baldry used,
; and in particular have magmax-magmin=0.5)
;
; INPUTS : 
; - mag : central mag of mag bin (of width 0.5)
; - cmin, cmac, cbins : min max of colors to sample with cbins bins.

; distribution of red galaxies : 
cred  = baldry_tfunc(mag,2.279,-0.037,-0.108,-19.81,0.96) ; mean color of reds
sred  = baldry_tfunc(mag,0.152,0.008,0.044,-19.91,0.94)   ; sdev of re

; distribution of blue galaxies :
cblue = baldry_tfunc(mag,1.79,-0.053,-0.363,-20.75,1.12)  ; mean color of blues
sblue = baldry_tfunc(mag,0.298,0.014,-0.067,-19.9,0.58)   ; sdev of blues

; color sampling and number densities : 
x     = findgen(cbins)/float(cbins-1L)*(cmax-cmin) + cmin
nred  = 0.921*0.00225*exp(-0.921*(mag+21.49)*(1-0.83))*exp(-exp(-0.921*(mag+21.49))) $ ; red galaxies r-band LF
        * exp(-0.5*(x-cred)^2/sred^2)/sred/sqrt(6.28)                                  ; Gaussian colour distribution
nblue =(0.921*0.00282*exp(-0.921*(mag+20.60)*(1+0.26))*exp(-exp(-0.921*(mag+20.60))) $
        +0.921*0.00235*exp(-0.921*(mag+20.60)*(1-1.35))*exp(-exp(-0.921*(mag+20.60))))$ ; blue galaxies r-band LF
       * exp(-0.5*(x-cblue)^2/sblue^2)/sblue/sqrt(6.28)  ; Gaussian colour distribution 
; NB: 
; 0.4*ln(10) = 0.921
; 0.00225 is phi* for early type
; -21.49  is M*   for early type
; -0.83 is the Schechter index for early type
; same things for late types, where a double Schechter function has been used

coldist = {cmin:cmin,cmax:cmax,nbins:cbins,col:x,nred:nred,nblue:nblue}
return,coldist


end

