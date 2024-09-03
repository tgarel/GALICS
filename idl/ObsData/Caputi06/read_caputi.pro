function read_caputi, filter, z
; J.E. Forero-Romero (2009)
;
; loads the bolometric IR luminosity fromCaputi'h 06 paper
;
; INPUTS
; -filter: the name of the filter
; -z: redshift
;

 spawn,'echo $GALICS_PATH',path
 dir = strtrim(path,2)+'/idl/ObsData/Caputi06/'
 filename = dir+filter+'_Caputi_z' + strtrim(string(FIX(z)),2)+".lf"
 openr,11,filename
 readf,11,n1,n2
 a=fltarr(n2,n1)
 readf,11,a
 close,11
 
 lf_data = {loglum : reform(a(0,*)),$    ; log10(\nu L_bol}/(h^-2 L_sun))
            logLF     : reform(a(1,*)),$ ; log10(dN/d(log10(L)) / h^3 Mpc^-3)
            logLFlow  : reform(a(2,*)),$
            logLFhigh : reform(a(3,*))} 
 return, lf_data
end
