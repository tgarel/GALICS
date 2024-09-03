function shen_lognormal,lx,mean,dev

; lx = ln(r)

a = 1.d0/(2.5066283d0*dev)   ; 2.5066283 is sqrt(2*pi)
b = (lx - mean)^2 / 2.0d0 / dev^2
y = a * exp(-b)

return, y

end
