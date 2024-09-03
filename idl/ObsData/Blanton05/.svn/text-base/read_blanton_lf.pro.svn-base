pro read_blanton_lf, band, lf

; band is one of u, g, r, i, z
; lf is a struct with 
; - lf.absmag      : the absolute mag (M-5log(h)), K+E corrected to z=0
; - lf.phi         : the number density in (Mpc/h)^-3 mag^-1 
; - lf.lg10phi_err : the error on phi, in log10.
; NB: although given as a function of h, Blanton uses h=1 ... 

spawn,'echo $GALICS_PATH',path
dir = strtrim(path,2)+'/idl/ObsData/Blanton05/'
filename = dir+'eep_lowz_'+strtrim(band,2)+'.drtwo14.sbweight.fits'
lf = mrdfits(filename,1)

end
