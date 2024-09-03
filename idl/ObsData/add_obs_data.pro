pro add_obs_data,thing

;;;;; Data from Baldry et al. Apj, 600, 681, 2004

check = strpos(thing,'baldry')
if(check ne -1) then begin
    readcol, '/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/baldry-bimodality-data.txt', $
      mag, col, num, err, $
             format = 'f, f, f, f', /silent
    magbinpos = strpos(thing, '-')
    magbin = strmid(thing,magbinpos,6)
    magbin = float(magbin)
    ind = where(mag eq magbin)
    if(ind(0) ne -1) then begin
        symbols, 2, 1
        normnum = num(ind)/float(total(num(ind)))
        normerr = err(ind)/float(total(num(ind)))
        oploterror, col(ind), normnum, normerr, psym=8, errsym=8, $
                    errthick = 4
        print, total(normnum)
        ;plots, 0.2, 0.11, psym = 8
        ;xyouts, 0.3, 0.11, 'Baldry et al. (2004)', color = !lgreen
    endif
endif

;;;;; Fig. 2 from Kauffmann et al. 2003, 341, 54

check = strpos(thing,'dstat')
if (check ne -1) then begin 
    readcol, '/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/FromGuinevere/'+thing, $
      xdata, ydata, format='f,f', /silent
    symbols, 2, 1
    oplot, xdata, ydata, psym=8, color = !black, thick=4
endif

if(thing eq 'tremonti') then begin
    readcol, '/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/tremonti.txt', $
      xx, yy2_5, yy16, yy50, yy84, yy97_5, FORMAT='F,F,F,F,F,F'    
    symbols, 2, 0.8
    plots, xx, yy50, psym=8,    color=!dgreen
    plots, xx, yy2_5,  color=!dgreen, thick=4
    plots, xx, yy16,   color=!dgreen, thick=4
    plots, xx, yy50,   color=!dgreen, thick=4
    plots, xx, yy84,   color=!dgreen, thick=4
    plots, xx, yy97_5, color=!dgreen, thick=4
    xyouts, 10., 8.2, 'Tremonti et al. (2004)', charthick=4, color=!dgreen
endif

if(thing eq 'gallazzi') then begin
    readcol, '/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/gallazzi.txt', $
      xx, yy50, yy16, yy84, tt50, tt16, tt84, FORMAT='F,F,F,F,F,F,F'
   symbols, 2, 0.8
   plots, xx, yy50, psym=8, color=!dgreen
   plots, xx, yy16, color=!dgreen, thick=4
   plots, xx, yy84, color=!dgreen, thick=4
   xyouts, 10., 8.2, 'Gallazzi et al. (2005)', charthick=4, color=!dgreen
endif

if (thing eq 'SDSS_u') then begin 
    openr,11,'/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/sdss_u_blanton.dat'
    pipo='' & readf,11,pipo
    a = fltarr(3,620) & readf,11,a
    close,11
    oploterr,a(0,*),a(1,*)*0.73^3,a(2,*)*0.73^3
    oplot,a(0,*),a(1,*)*0.73^3,thick=4,color=!dgreen
    xyouts, -20, 3e-6, 'Blanton',color=!dgreen
endif
if (thing eq 'SDSS_g') then begin 
    openr,11,'/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/sdss_g_blanton.dat'
    pipo='' & readf,11,pipo
    a = fltarr(3,740) & readf,11,a
    close,11
    oploterr,a(0,*),a(1,*)*0.73^3,a(2,*)*0.73^3
    oplot,a(0,*),a(1,*)*0.73^3,thick=4,color=!dgreen
    xyouts, -20, 3e-6, 'Blanton',color=!dgreen
endif
if (thing eq 'SDSS_r') then begin 
    openr,11,'/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/sdss_r_blanton.dat'
    pipo='' & readf,11,pipo
    a = fltarr(3,620) & readf,11,a
    close,11
    oploterr,a(0,*),a(1,*)*0.73^3,a(2,*)*0.73^3
    oplot,a(0,*),a(1,*)*0.73^3,thick=4,color=!dgreen
    xyouts, -20, 3e-6, 'Blanton',color=!dgreen
endif
if (thing eq 'SDSS_i') then begin 
    openr,11,'/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/sdss_i_blanton.dat'
    pipo='' & readf,11,pipo
    a = fltarr(3,830) & readf,11,a
    close,11
    oploterr,a(0,*),a(1,*)*0.73^3,a(2,*)*0.73^3
    oplot,a(0,*),a(1,*)*0.73^3,thick=4,color=!dgreen
    xyouts, -20, 3e-6, 'Blanton',color=!dgreen
endif
if (thing eq 'SDSS_z') then begin 
    openr,11,'/afs/mpa/sim/Millennium/SA-lucia/Extinction/ObsData/sdss_z_blanton.dat'
    pipo='' & readf,11,pipo
    a = fltarr(3,830) & readf,11,a
    close,11
    oploterr,a(0,*),a(1,*)*0.73^3,a(2,*)*0.73^3
    oplot,a(0,*),a(1,*)*0.73^3,thick=4,color=!dgreen
    xyouts, -20, 3e-6, 'Blanton',color=!dgreen
endif
end



