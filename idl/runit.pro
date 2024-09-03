pro runit,runs
  
  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
;  run.outputdir = '/data/blaizot/GM-RUNS/512-100h-1Mpc/BestFit/'
  run.outputdir = '/Users/blaizot/Documents/ASTRO/data/GM-RUNS/512-100h-1Mpc-W3/BestFit/'
  run.snapzdir=run.outputdir +'/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='bestFitW3'
  run.volumempc = 136.9863^3. 
;  run.volumempc = 142.9863^3. 
  run.hubble=0.73
  
  set_plot,'ps'
  device,/color,bits_per_pixel=8,filename=run.name+'.ps'
  loadct,39
;;   plot_gal_props,run
;;  ts = z2ts(3.0,run.snapzdir)
;;  plot_haloprops,run,ts,titre='z = 3'
;;   plot_haloprops,run,6;,/nosubdir
;;   plot_haloprops,run,8;,/nosubdir
;;   plot_haloprops,run,10;,/nosubdir
;;   plot_haloprops,run,12;,/nosubdir
;;   plot_haloprops,run,14;,/nosubdir
;;   plot_haloprops,run,16;,/nosubdir
;;   plot_haloprops,run,18;,/nosubdir
;;   plot_haloprops,run,20;,/nosubdir
  plot_high_z_lfs,run;,/nosubdir
;;  plot_bookkeep,run,16;,/nosubdir
  device,/close
stop







for i = 0,n_elements(runs)-1 do begin 
   irun = runs(i)
   run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
   run.outputdir='/data/forero/GM-RUNS/'+strtrim(irun,2)+'/'
   ;run.outputdir = '/data/blaizot/TEST/newKennicutt/fullRun1/'
   run.snapzdir=run.outputdir;+'1/'
   run.dirmin=1
   run.dirmax=15
   run.name='run'+strtrim(irun,2)
   run.volumempc = 136.9863^3. / 15.0d0
   run.hubble=0.73

   set_plot,'ps'
   device,/color,bits_per_pixel=8
   loadct,39
   plot_haloprops,run,6,/nosubdir
   plot_haloprops,run,8,/nosubdir
   plot_haloprops,run,10,/nosubdir
   plot_haloprops,run,12,/nosubdir
   plot_haloprops,run,14,/nosubdir
   plot_haloprops,run,16,/nosubdir
   plot_haloprops,run,18,/nosubdir
   plot_haloprops,run,20,/nosubdir
   device,/close
;;    set_plot,'ps'
;;    device,filename='run'+strtrim(irun,2)+'.ps'
;;    plot_high_z_lfs,run,/nosubdir
;;    plot_bookkeep,run,16,/nosubdir
;;    device,/close
   
endfor

;  calibration_plots,run,'calib_runw3.ps'
;  ir_calibration_plots,run,'ir_calib_runw3.ps'
;  mass_sfr,run,'mass_sfrw3.ps'
stop

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/run105/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/run105/1/'
  run.dirmin=1
  run.dirmax=10
  run.name='TEST'
  run.volumempc = 142.85714^3. / 60. * 10.
  run.hubble=0.7
  calibration_plots,run,'calib_run105.ps'
  ir_calibration_plots,run,'ir_calib_run105.ps'
  mass_sfr,run,'mass_sfr_105.ps'

stop

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/LowFeedback/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/LowFeedback/1/'
  run.dirmin=1
  run.dirmax=10
  run.name='TEST'
  run.volumempc = 142.85714^3. / 60. * 10.
  run.hubble=0.7
  calibration_plots,run,'calib_lowfeedback_hr.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_SN_FEEDBACK/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_SN_FEEDBACK/1/'
  run.dirmin=1
  run.dirmax=10
  run.name='TEST'
  run.volumempc = 142.85714^3. / 60. * 10.
  run.hubble=0.7
  calibration_plots,run,'calib_nofeedback_hr.ps'


  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/256-100h-1Mpc/RUN_TEST/'
  run.snapzdir='/data/blaizot/GM-RUNS/256-100h-1Mpc/RUN_TEST/1/'
  run.dirmin=1
  run.dirmax=24
  run.name='TEST'
  run.volumempc = 142.85714^3. /24. * 24.
  run.hubble=0.7
  calibration_plots,run,'calib_test_lr.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/RUN_TEST/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/RUN_TEST/1/'
  run.dirmin=1
  run.dirmax=60
  run.name='TEST'
  run.volumempc = 142.85714^3. /24. * 24.
  run.hubble=0.7
  calibration_plots,run,'calib_test_hr.ps'

stop

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NEW_lowall/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NEW_lowall/1/'
  run.dirmin=1
  run.dirmax=60
  run.name='normal run'
  run.volumempc = 142.85714^3. ;/ 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_new_lowall.ps'
stop
  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_SN_FEEDBACK/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_SN_FEEDBACK/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_no_sn_feedbacl.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_DISC_INSTAB/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_DISC_INSTAB/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_no_disc_instab.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_MERGERS/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NO_MERGERS/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='no disc instab'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_no_mergers.ps'
 
  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_highSN/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_highSN/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run (high SN)'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_normal_highSN.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_lowSN/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_lowSN/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run (low SN)'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_normal_lowSN.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_lowSF/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_lowSF/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run (low SF)'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_normal_lowSF.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_highSF/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_highSF/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run (high SF)'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_normal_highSF.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_HighMcrit/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_HighMcrit/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run (high SF)'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_normal_HighMcrit.ps'

  run = {outputdir:' ',dirmin:0,dirmax:0,snapzdir:' ',name:' ',volumempc:0.0, hubble:0.0}
  run.outputdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_LowMcrit/'
  run.snapzdir='/data/blaizot/GM-RUNS/512-100h-1Mpc/NORMAL_LowMcrit/1/'
  run.dirmin=1
  run.dirmax=15
  run.name='normal run (high SF)'
  run.volumempc = 142.85714^3. / 60. * 15.
  run.hubble=0.7
  calibration_plots,run,'calib_hr_normal_LowMcrit.ps'


end

