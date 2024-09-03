pro titi

  sTemplate = ASCII_TEMPLATE(gabO4_2800_z3)
  
  rdfloat, 'gabO4_2800_z3.txt',FIELD1,FIELD2,FIELD3,FIELD4,NUMLINE =4 ; z=6

  plot,[0,0],thick=3.5,/NODATA,xr=[-23.,-18.],yr=[0.00001,0.1],xtitle='M',ytitle='log(N / Mpc^-3 / mag)',charsize=1.5,charthick=2.,TITLE='UV LF distribution function',/ylog

  oplot,FIELD1,FIELD2,psym=2,thick=4
  errplot,FIELD1,FIELD2+FIELD3,FIELD2+FIELD4

end
