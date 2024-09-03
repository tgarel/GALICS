pro read_shell_params, outputdir, timestep,shell_prop

  file = outputdir + 'shell_props.'+string(timestep,format='(i3.3)')
  openr,11,file

  r_sh    = 0.0d0
  v_sh    = 0.0d0
  age_sh  = 0.0d0
  nh_sh   = 0.0d0
  tau_sh  = 0.0d0
  M_sh    = 0.0d0
  M_shz   = 0.0d0
  cf_sh   = 0.0d0
  zmet_sh = 0.0d0
  onoff      = 'A'

;;   i = 0LL
;;   timestep = 0L
  ngal     = 0LL
  st       = 0
  aexp     = 0.0d0
  age_univ = 0.0d0

  readf,11,format='(i3,1x,i6,1x,e14.6,1x,e14.6)',st,ngal,aexp,age_univ
  shell_prop       = replicate({r_shell:0.0d0,v_shell:0.0d0,age_shell:0.0d0,nh_shell:0.0d0,tau_shell:0.0d0,M_shell:0.0d0,M_shellz:0.0d0,cf_shell:0.0d0,zmet_shell:0.0d0,onoff_shell:'A'},ngal)
 print,'NGAL',ngal
  for igal = 0LL, ngal - 1LL do begin ;ngal - 1LL
     readf,11,format='(9e14.6,A1)', r_sh,v_sh,age_sh,nh_sh,tau_sh,M_sh,M_shz,cf_sh,zmet_sh,onoff

     shell_prop(igal).r_shell     = r_sh
     shell_prop(igal).v_shell     = v_sh
     shell_prop(igal).age_shell   = age_sh
     shell_prop(igal).nh_shell    = nh_sh
     shell_prop(igal).tau_shell   = tau_sh
     shell_prop(igal).M_shell     = M_sh
     shell_prop(igal).M_shellz    = M_shz
     shell_prop(igal).cf_shell    = cf_sh
     shell_prop(igal).zmet_shell  = zmet_sh
     shell_prop(igal).onoff_shell = onoff

  endfor

  close,11

end
