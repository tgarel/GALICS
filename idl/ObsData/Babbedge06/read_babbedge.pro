function read_babbedge, filter, z
; J.E. Forero-Romero (2009)
;
; loads the luminosity function data for redshifts 0, 1, 2 and 3
; for the filters IRAC_3_6mic,  IRAC_4_5mic, IRAC_5_8mic, IRAC_8_0mi and
; MIPS_24mic
;
; INPUTS
; -filter: the name of the filter
; -z: redshift
;

 spawn,'echo $GALICS_PATH',path
 dir = strtrim(path,2)+'/idl/ObsData/Babbedge06/'
 filename = dir+filter+'_Babbedge_z' + strtrim(string(FIX(z)),2)+".lf"
 openr,11,filename
 readf,11,n1,n2
 a=fltarr(n2,n1)
 readf,11,a
 close,11
 
 lf_data = {loglum : reform(a(0,*)),$    ; log10(\nu L_{\nu}/(h^-2 L_sun))
            logLF     : reform(a(1,*)),$ ; log10(dN/d(log10(flux)) / h^3 Mpc^-3)
            logLFlow  : reform(a(2,*)),$
            logLFhigh : reform(a(3,*))} 
 return, lf_data
end
