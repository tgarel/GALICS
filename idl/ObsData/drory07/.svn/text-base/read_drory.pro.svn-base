function read_drory, survey, z
; J.E. Forero-Romero (2009)
;
; loads the m_star functions from Drory et al 07 paper
;
; INPUTS
; -survey: the name of the survey (GOODS or FDF)
; -z: redshift
;

 spawn,'echo $GALICS_PATH',path
 dir = strtrim(path,2)+'/idl/ObsData/drory07/' 
 filename = dir+'drory_'+survey+'_mstar_z'+strtrim(string(FIX(z)),2)+"-"+strtrim(string(FIX(z+1)),2)+".txt"
 openr,11,filename
 readf,11,n1,n2
 a=fltarr(n2,n1)
 readf,11,a
 close,11
 
 lf_data = {loglum : reform(a(0,*)),$    ; log10(\nu L_bol}/(M_sun))
            logLF     : reform(a(1,*)),$ ; log10(dN/d(log10(L)) / Mpc^-3)
            logLFlow  : reform(a(2,*)),$
            logLFhigh : reform(a(3,*))} 
 return, lf_data
end
