; code provided by A. Cattaneo (2008)

function Tfunc,mag,p0,p1,q0,q1,q2
; Tfunc is a function of r-band magnitude 
; it contains five free parameters
; Baldry et al have successfully used this function to fit :
; the mean colour of early type galaxies
; the mean colour of late type galaxies
; the standard deviation from mean colour in late type
; the standard deviation from mean colour in early type
;;Tfunc=p0+p1*(-23.25+0.5*magin+20.)+q0*tanh((-23.25+0.5*magin-q1)/q2)
Tfunc=p0+p1*(mag+20.)+q0*tanh((mag-q1)/q2)
return, Tfunc
end

pro SDSSbimodality_plot

; magnitudes at which the colour distribution is to be computed
nmags=10
mag=-23.25+0.5*findgen(nmags)
cet=Tfunc(mag,2.279,-0.037,-0.108,-19.81,0.96) ; mean colour early type
clt=Tfunc(mag,1.790,-0.053,-0.363,-20.75,1.12) ; mean colour late type
slt=Tfunc(mag,0.298,0.014,-0.067,-19.9,0.58)   ; standard deviation from mean colour in late type
set=Tfunc(mag,0.152,0.008,0.044,-19.91,0.94)   ; standard deviation from mean colour in early type

; colours at which the colour distribution is to be computed
ncols=350
col=0.5+0.01*findgen(ncols)

; joint colour-magnitude distribution
colmagdist=fltarr(ncols,nmags)
for i=0,nmags-1 do begin 
   for j=0,ncols-1 do begin
       red_colmagdist = 0.921*0.00225*exp(-0.921*(mag(i)+21.49)*(1-0.83))*exp(-exp(-0.921*(mag(i)+21.49)))$ ; red galaxies r-band LF
                        * exp(-0.5*(col(j)-cet(i))^2/set(i)^2)/set(i)/sqrt(6.28)                            ; Gaussian colour distribution
       blu_colmagdist =(0.921*0.00282*exp(-0.921*(mag(i)+20.60)*(1+0.26))*exp(-exp(-0.921*(mag(i)+20.60)))$
                       +0.921*0.00235*exp(-0.921*(mag(i)+20.60)*(1-1.35))*exp(-exp(-0.921*(mag(i)+20.60))))$; blue galaxies r-band LF
                        * exp(-0.5*(col(j)-clt(i))^2/slt(i)^2)/slt(i)/sqrt(6.28)                            ; Gaussian colour distribution 
       colmagdist(j,i)= red_colmagdist + blu_colmagdist  
       ; 0.4*ln(10) = 0.921
       ; 0.00225 is phi* for early type
       ; -21.49  is M*   for early type
       ; -0.83 is the Schechter index for early type
       ; same things for late types, where a double Schechter function has
       ; been used
   endfor
endfor  

;plotting routine
!p.multi=[0,2,5,0,0]
!p.thick=3. 
!x.thick=3. 
!y.thick=3.
!p.charthick=3.
orig_device=!d.name
set_plot,'PS'
device,filename='colour_histograms.ps',bits=8,/encapsulated
for i=0,9 do begin
   plot,col,colmagdist(*,i)
endfor
device,/close
set_plot,'x'
end
